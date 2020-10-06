"""
Hotspot Finder identifies mutational hotspots and generates basic annotations
"""

# Import modules
from collections import defaultdict
import gzip
import itertools
import logging
import os
import sys

from bgconfig import BGConfig
import bgdata
import bgreference as bgref
from bgparsers import readers
import click
import daiquiri
from intervaltree import IntervalTree
import tabix


# Global variables
VERSION = '0.1.0'
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
LOGS = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}


def load_configuration(config_file, override=None):
    """
    Load the configuration file and checks the format.

    Args:
        config_file: configuration file path
        override: override values

    Returns:
        :class:`bgconfig.BGConfig`: configuration as a :obj:`dict`

    """
    config_template = f'{config_file}.template'

    try:
        return BGConfig(
            config_template, config_file=config_file, use_env_vars=True, override_values=override, unrepr=False)
    except ValueError as e:
        logger.error(e)
        sys.exit(-1)


class HotspotFinder:
    """Class to identify and annotate hotspots with basic information"""

    def __init__(self,
                 input_file,
                 mappable_regions,
                 blacklisted_regions,
                 population_variants,
                 genomic_elements,
                 output_file_results,
                 output_file_warning,
                 output_format,
                 is_gzip,
                 hotspot_mutations,
                 split_alternates,
                 remove_unknown_nucleotides,
                 remove_nonannotated_hotspots,
                 genome,
                 group_by,
                 cores
                 ):
        """
        Initialize HotspotFinder class

        Args:
            input_file (str): path to input mutations data
            mappable_regions (str): file with regions of high mappability
            blacklisted_regions (str): file with artifact regions of low mappability
            population_variants (str): file with polimorphisms data to insersect
            genomic_elements (str): file with genomic elements data to intersect
            output_file_results (str): path to output file
            output_file_warning (str): path to genomic positions where warning is found
            output_format (str): format of output file (TSV or VCF)
            is_gzip (bool): GZIP output files
            hotspot_mutations (int): cutoff of mutations to define a hotspot
            split_alternates (bool): compute hotspots per nucleotide alternate independently
            remove_unknown_nucleotides (bool): remove hotspots with N nucleotides in their reference context
            remove_nonannotated_hotspots (bool): remove hotspots not overlapping genomic elements
            genome (str): reference genome
            group_by (str): name of the column to group hotspots identification
            cores (int): number of cpu

        Returns:
            None

        """

        # Input output files
        self.input_file = input_file
        self.input_file_name = input_file.split('/')[-1].split('.')[0]
        self.output_file_hotspots = output_file_results
        self.output_file_warning = output_file_warning

        # Mappability data
        self.mappable_regions_file = mappable_regions
        self.blacklisted_regions_file = blacklisted_regions

        # Variation data
        self.variation_directory = population_variants

        # Genomic elements data
        self.genomic_elements = genomic_elements
        self.regions_tree = None
        self.genomic_elements_list = [
            'cds',
            '5utr',
            '3utr',
            'proximal_promoters',
            'distal_promoters',
            'introns'
        ]
        self.genomic_elements_priority = {e: self.genomic_elements_list.index(e) for e in self.genomic_elements_list}

        # Params
        self.output_format = output_format
        self.open_function = gzip.open if is_gzip else open
        self.write_mode = 'wt' if is_gzip else 'w'
        self.hotspot_mutations = hotspot_mutations
        self.split_alternates = split_alternates
        self.remove_unknown_nucleotides = remove_unknown_nucleotides
        self.remove_nonannotated_hotspots = remove_nonannotated_hotspots
        self.genome = genome
        self.group_by = group_by
        self.cores = cores

        # Initialize variables to load data
        self.mutations_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.original_reference = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
        self.cohort_total_mutations = defaultdict(lambda: defaultdict(int))
        self.cohort_to_sample = defaultdict(set)
        self.cohort_to_mutation_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        self.cohort_to_mutation_alts = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.hotspots = defaultdict(lambda: defaultdict(dict))
        self.hotspots_samples = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
        self.warning_chr_position = set()

    @staticmethod
    def load_genomic_elements(files):
        """
        Load genomic element regions into intervaltree
        Args:
            files (list): list of tuples containing element_type and path(s) to file(s)

        Returns:
            tree (dict): tree with regions, keys are chromosomes, data are regions
        """

        tree = defaultdict(IntervalTree)
        for genomic_element, file in files:
            with gzip.open(file, 'rt') as fd:
                next(fd)
                for line in fd:
                    chrom, start, end, strand, gene_id, transcript_id, symbol = line.strip().split('\t')
                    tree[chrom].addi(int(start), int(end) + 1, f'{symbol}::{gene_id}::{transcript_id}::{genomic_element}')  # +1
        return tree

    @staticmethod
    def load_mappability(file, chr_format='chrN'):
        """
        Load mappability regions into intervaltree
        Args:
            file (path): path to file containing mappability data
            chr_format (str): chromosome format used in file. HotspotFinder works with format '1' instead of 'chr1'

        Returns:
            tree (dict): tree with regions, keys are chromosomes, data are regions
        """

        tree = defaultdict(IntervalTree)
        trim = 3 if chr_format == 'chrN' else None
        with gzip.open(file, 'rt') as fd:
            for line in fd:
                chrom, start, end = line.strip().split('\t')
                tree[chrom[trim:]].addi(int(start), int(end) + 1, '')  # +1 interval
        return tree

    def parse_mutations(self):
        """
        Load mutations file into dictionaries containing number of mutations and alternates

        Returns:
            None

        """

        chromosomes = list(map(str, range(1, 23))) + ['X', 'Y']
        mutations_ref_nomatch = 0

        # Read mutations
        for row in readers.variants(
                file=self.input_file,
                required=['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE'],
                extra=['GROUP', 'GROUP_BY', 'COHORT', 'CANCER_TYPE', 'PLATFORM']
        ):
            chromosome = row['CHROMOSOME']
            position = str(row['POSITION'])
            ref = row['REF']
            alt = row['ALT']
            sample = row['SAMPLE']
            chr_position = '{}_{}'.format(chromosome, position)
            # Identify group
            # If no group, hotspots are computed using the whole input file
            if not self.group_by:
                cohort = self.input_file_name
            else:
                cohort = row[self.group_by]
            # Keep track of samples from each group
            self.cohort_to_sample[cohort].add(sample)
            # Keep track of how many mutations each group has
            self.cohort_total_mutations['total'][cohort] += 1

            # Read mutations in autosomal + sexual chromosomes
            if chromosome in set(chromosomes):
                if ref != alt:
                    # Read substitutions of any length
                    if ref != '-' and alt != '-':
                        if len(alt) == 1:
                            # Check reference
                            if ref == bgref.refseq(self.genome, chromosome, int(position), 1):
                                self.mutations_dict['snv'][sample][chr_position].append(alt)
                                self.cohort_total_mutations['snv'][cohort] += 1
                                self.hotspots_samples[cohort]['snv'][chr_position].add(sample)
                            else:
                                mutations_ref_nomatch += 1
                                self.cohort_total_mutations['total'][cohort] -= 1
                        else:
                            self.mutations_dict['mnv'][sample][chr_position].append(alt)
                            self.cohort_total_mutations['mnv'][cohort] += 1
                            self.hotspots_samples[cohort]['mnv'][chr_position].add(sample)

                    # Read indels of any length
                    else:
                        # Insertions
                        if ref == '-':
                            self.mutations_dict['ins'][sample][chr_position].append(alt)
                            self.cohort_total_mutations['ins'][cohort] += 1
                            self.hotspots_samples[cohort]['ins'][chr_position].add(sample)
                        # Deletion
                        elif alt == '-':
                            self.mutations_dict['del'][sample][chr_position].append(alt)
                            self.cohort_total_mutations['del'][cohort] += 1
                            self.hotspots_samples[cohort]['del'][chr_position].add(sample)
                            self.original_reference['del'][sample][chr_position].add(ref)

        if mutations_ref_nomatch > 0:
            logger.warning(f'A total of {mutations_ref_nomatch} SNVs REF nucleotides do not match the reference. '
                           f'These mutations are discarded from the analysis')

        for muttype, data in self.cohort_total_mutations.items():
            for cohort, nmuts in data.items():
                logger.info(f'Input {muttype} mutations in {cohort} = {nmuts}')

        # Check mutations per sample and add to cohort_to_mutation dicts
        # Write warning positions
        header = ['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE', 'WARNING', 'SKIP']
        with self.open_function(self.output_file_warning, self.write_mode) as ofd:
            ofd.write('{}\n'.format('\t'.join(header)))
            warning_samples_n = 0
            for cohort, set_of_samples in self.cohort_to_sample.items():
                for sample in set_of_samples:
                    total_muts_in_sample = defaultdict(list)
                    for muttype in ('snv', 'mnv', 'ins', 'del'):
                        for chr_pos, alts in self.mutations_dict[muttype][sample].items():
                            total_muts_in_sample[chr_pos] += [(a, muttype) for a in alts]
                    # Warning mutations, when there's more than one alternate in the same position
                    warning_mutations = [(k, v) for k, v in total_muts_in_sample.items() if len(v) > 1]

                    logger.debug(sample)
                    logger.debug(total_muts_in_sample)
                    logger.debug(warning_mutations)

                    if warning_mutations:
                        warning_samples_n += 1
                        # Load warning positions
                        for chr_position, alts in warning_mutations:
                            chromosome, position = chr_position.split('_')
                            alts_unique = set(alts)
                            alts_simplified = [a for a, muttype in alts]
                            self.warning_chr_position.add(chr_position)
                            if len(alts_unique) == 1:
                                logger.debug(
                                    f'Sample "{sample}" position chr{chr_position} has 2 alternates: {alts_simplified}')
                                alt, muttype = list(alts_unique)[0]
                                self.cohort_to_mutation_alts[cohort][muttype][f'{chr_position}_{muttype}'] += [alt]
                                self.cohort_to_mutation_counts[cohort][muttype][f'{chr_position}_{muttype}'] += 1
                                reference_n = bgref.refseq(self.genome, chromosome, int(position), 1)
                                ofd.write('{}\n'.format('\t'.join([
                                    chromosome, position, reference_n, ','.join(alts_simplified), sample, 'warning_1', 'False'
                                ])))
                            elif len(alts_unique) == 2:
                                logger.debug(
                                    f'Sample "{sample}" position chr{chr_position} has 2 or 3 alternates: {alts_simplified}. '
                                    f'Number of mutations and number of mutated samples will not match.')
                                reference_n = bgref.refseq(self.genome, chromosome, int(position), 1)
                                ofd.write('{}\n'.format('\t'.join([
                                    chromosome, position, reference_n, ','.join(alts_simplified), sample, 'warning_2', 'False'
                                ])))
                                for mutation in alts_unique:
                                    alt, muttype = mutation
                                    self.cohort_to_mutation_alts[cohort][muttype][f'{chr_position}_{muttype}'] += [alt]
                                    self.cohort_to_mutation_counts[cohort][muttype][f'{chr_position}_{muttype}'] += 1
                            else:
                                logger.debug(
                                    f'Sample "{sample}" position chr{chr_position} has 3 or more different alternates: '
                                    f'{alts_simplified}. Mutations are skipped from analysis')
                                reference_n = bgref.refseq(self.genome, chromosome, int(position), 1)
                                ofd.write('{}\n'.format('\t'.join([
                                    chromosome, position, reference_n, ','.join(alts_simplified), sample, 'warning_3', 'True'
                                ])))
                                self.hotspots_samples[cohort][muttype][chr_position].discard(sample)

                    # Load non-warning positions
                    for chr_position, mutation in total_muts_in_sample.items():
                        if not chr_position in self.warning_chr_position:
                            muttype = mutation[0][1]
                            self.cohort_to_mutation_alts[cohort][muttype][f'{chr_position}_{muttype}'] += [mutation[0][0]]
                            self.cohort_to_mutation_counts[cohort][muttype][f'{chr_position}_{muttype}'] += 1

        if warning_samples_n > 0:
            logger.warning(f'A total of {warning_samples_n} samples '
                           f'contain mutations flagged with warnings. Please check '
                           f'{self.output_file_warning}')

    def find_hotspots(self):
        """
        Identifies hotspots and generates raw hotspots file

        Returns:
            None
        """

        if self.split_alternates is False:
            for cohort, data in self.cohort_to_mutation_counts.items():
                for muttype, mutated_positions in data.items():
                    self.hotspots[cohort][muttype] = {k: v for k, v in mutated_positions.items() if
                                                      v >= self.hotspot_mutations}
        else:
            for cohort, data in self.cohort_to_mutation_counts.items():
                for muttype, mutated_positions in data.items():
                    if muttype == 'snv':
                        for hotspot_id, n_samples in mutated_positions.items():
                            list_of_alternates = self.cohort_to_mutation_alts[cohort][muttype][hotspot_id]
                            for alternate in set(list_of_alternates):
                                count = list_of_alternates.count(alternate)
                                if count >= self.hotspot_mutations:
                                    self.hotspots[cohort][muttype][f'{hotspot_id}>{alternate}'] = count
                                    self.cohort_to_mutation_alts[cohort][muttype][f'{hotspot_id}>{alternate}'] = [alternate] * count
                    else:
                        self.hotspots[cohort][muttype] = {k: v for k, v in mutated_positions.items() if
                                                          v >= self.hotspot_mutations}

    def write_hotspots(self):
        """
        Writes output file for HotspotFinder

        Returns:
            None
        """
        nucleotides = ('A', 'C', 'G', 'T')
        header = [
            'CHROMOSOME',
            'POSITION',
            'CHR_POS',
            'HOTSPOT_ID',
            'MUT_TYPE',
            'COHORT',
            'N_MUTATIONS',
            'N_MUTATED_SAMPLES',
            'FRAC_MUTATED_SAMPLES',
            'REF',
            'ALT',
            'ALT_COUNTS',
            'FRAC_ALT',
            'CONTEXT_3',
            'CONTEXT_5',
            'GENOMIC_ELEMENT',
            'SYMBOL',
            'GENE_ID',
            'TRANSCRIPT_ID',
            'GENOMIC_REGION',
            'GENOMIC_REGION_PRIORITY',
            'CODING_NONCODING',
            'OVERLAP_WARNING_POSITION',
            'MAPPABILITY',
            'MAPPABILITY_BLACKLIST',
            'VARIATION_AF',
            'HOTSPOTFINDER_FILTERS',
            'N_COHORT_SAMPLES',
            'N_COHORT_MUTATIONS_TOTAL',
            'N_COHORT_MUTATIONS_SNV',
            'N_COHORT_MUTATIONS_MNV',
            'N_COHORT_MUTATIONS_INS',
            'N_COHORT_MUTATIONS_DEL',
            'MUTATED_SAMPLES',
            'MUTATED_SAMPLES_ALTS'
        ]

        variation_tb_dict = {}
        for chromosome in list(map(str, range(1, 23))) + ['X', 'Y']:
            variation_tb_dict[chromosome] = tabix.open(f'{self.variation_directory}/gnomad.genomes.r3.0.sites.chr{chromosome}.af_0.01.tsv.gz')
        hotspot_count = 0

        with self.open_function(self.output_file_hotspots, self.write_mode) as ofd:
            ofd.write('{}\n'.format('\t'.join(header)))
            for cohort, data in self.hotspots.items():
                n_cohort_samples = len(self.cohort_to_sample[cohort])
                n_cohort_mut_total = self.cohort_total_mutations['total'][cohort]
                n_cohort_mut_snv = self.cohort_total_mutations['snv'][cohort]
                n_cohort_mut_mnv = self.cohort_total_mutations['mnv'][cohort]
                n_cohort_mut_ins = self.cohort_total_mutations['ins'][cohort]
                n_cohort_mut_del = self.cohort_total_mutations['del'][cohort]
                for muttype, hotspots in data.items():
                    for hotspot_id_type, n_mutations in sorted(hotspots.items(), key=lambda item: item[1], reverse=True):
                        chromosome, position, _ = hotspot_id_type.split('>')[0].split('_')    # remove SNV alternate
                        chr_pos = f'{chromosome}_{position}'

                        # Mutated samples
                        mut_samples = set(list(self.hotspots_samples[cohort][muttype][chr_pos]))
                        # Alternates
                        mut_samples_to_alt = {}
                        for sample in mut_samples:
                            mut_samples_to_alt[sample] = self.mutations_dict[muttype][sample][chr_pos]
                        list_of_alternates = self.cohort_to_mutation_alts[cohort][muttype][hotspot_id_type]

                        alternates_counts = []
                        alternates_fractions = []
                        if muttype == 'snv':
                            if self.split_alternates:
                                alternate = hotspot_id_type.split('>')[1]
                                mut_samples = set()
                                for sample, sample_alternates in mut_samples_to_alt.items():
                                    for unique_alternate in sample_alternates:
                                        if unique_alternate == alternate:
                                            mut_samples.add(sample)
                                mut_samples_to_alt = {k: v for k, v in mut_samples_to_alt.items() if k in mut_samples}
                                list_of_alternates = list(itertools.chain.from_iterable(list(mut_samples_to_alt.values())))
                            for nucleotide in nucleotides:
                                count = list_of_alternates.count(nucleotide)
                                alternates_counts.append(f'{nucleotide}={count}')
                                alternates_fractions.append(f'{nucleotide}={count/n_mutations}')
                        else:
                            alternates_counts = [f'{a}={list_of_alternates.count(a)}' for a in set(list_of_alternates)]
                            alternates_fractions = [f'{a}={list_of_alternates.count(a) / n_mutations}' for a in set(list_of_alternates)]

                        alternates = ','.join(sorted(list(set(list_of_alternates))))
                        alternates_counts = ','.join(alternates_counts)
                        alternates_fractions = ','.join(alternates_fractions)
                        n_mut_samples = len(mut_samples)
                        frac_mut_samples = str(len(mut_samples) / n_cohort_samples)
                        mut_samples_to_alt = ';'.join([f'{k}::{",".join(v)}' for k, v in mut_samples_to_alt.items()])
                        mut_samples = ';'.join(list(mut_samples))

                        # Sequence info
                        pentamer_sequence = bgref.refseq(self.genome, chromosome, int(position) - 2, 5)
                        trimer_sequence = pentamer_sequence[1:4]
                        if muttype == 'ins':
                            ref = '-'
                        elif muttype == 'del':
                            ref = set()
                            for sample in mut_samples.split(';'):
                                ref.update(self.original_reference['del'][sample][chr_pos])
                            ref = ','.join(list(ref))
                        else:
                            ref = pentamer_sequence[2]

                        if self.remove_unknown_nucleotides and 'N' in trimer_sequence or 'N' in pentamer_sequence:
                            break

                        # Warning flag: at least one sample in the dataset contains a warning in this position
                        warning_flag = 'True' if chr_pos in self.warning_chr_position else 'False'

                        # HotspotFinder filters
                        hotspotfinder_filters = 2

                        # Mappability
                        map_data = 'low'
                        for intersect in self.mappable_regions_tree[chromosome][int(position)]:
                            if intersect:
                                map_data = 'high'
                                hotspotfinder_filters += 1
                                break

                        # Mappability blacklisted regions
                        blacklist_data = 'PASS'
                        for intersect in self.blacklisted_regions_tree[chromosome][int(position)]:
                            if intersect:
                                blacklist_data = 'FAIL'
                                hotspotfinder_filters -= 1
                                break

                        # Variation
                        variation_tb = variation_tb_dict[chromosome]
                        var_data = 'PASS'
                        try:
                            for info in variation_tb.query(f'chr{chromosome}', int(position), int(position)):
                                var_data = 'FAIL'
                                hotspotfinder_filters -= 1
                                break
                        except tabix.TabixError:
                            var_data = 'PASS'

                        hotspotfinder_filters = 'PASS' if hotspotfinder_filters == 3 else 'FAIL'

                        # Genomic elements
                        genomic_elements_full = []
                        symbol = set()
                        geneid = set()
                        transcriptid = set()
                        genomic_elements_type = set()
                        for intersect in self.regions_tree[chromosome][int(position)]:
                            if intersect:
                                genomic_elements_full += [intersect.data]
                                isymbol, igeneid, itranscriptid, igenomicelement = intersect.data.split('::')
                                symbol.add(isymbol)
                                geneid.add(igeneid)
                                transcriptid.add(itranscriptid)
                                genomic_elements_type.add(igenomicelement)

                        if genomic_elements_full:
                            genomic_elements_full = ';'.join(genomic_elements_full)
                            symbol = ';'.join(sorted(list(symbol)))
                            geneid = ';'.join(sorted(list(geneid)))
                            transcriptid = ';'.join(sorted(list(transcriptid)))
                            genomic_elements_type = ';'.join(sorted(
                                list(genomic_elements_type), key=lambda x: self.genomic_elements_priority.get(x, 999)))
                            genomic_elements_type_priority = genomic_elements_type.split(';')[0]
                        else:
                            genomic_elements_full = 'None'
                            symbol = 'None'
                            geneid = 'None'
                            transcriptid = 'None'
                            genomic_elements_type = 'None'
                            genomic_elements_type_priority = 'None'

                        if os.path.isfile(self.genomic_elements):
                            coding_noncoding = 'Unknown'
                        else:
                            coding_noncoding = 'CODING' if 'cds' in genomic_elements_type else 'NONCODING'

                        # Merge data
                        data_to_write = list(map(str,
                                                 [
                                                     chromosome,
                                                     position,
                                                     chr_pos,
                                                     hotspot_id_type,
                                                     muttype,
                                                     cohort,
                                                     n_mutations,
                                                     n_mut_samples,
                                                     frac_mut_samples,
                                                     ref,
                                                     alternates,
                                                     alternates_counts,
                                                     alternates_fractions,
                                                     trimer_sequence,
                                                     pentamer_sequence,
                                                     genomic_elements_full,
                                                     symbol,
                                                     geneid,
                                                     transcriptid,
                                                     genomic_elements_type,
                                                     genomic_elements_type_priority,
                                                     coding_noncoding,
                                                     warning_flag,
                                                     map_data,
                                                     blacklist_data,
                                                     var_data,
                                                     hotspotfinder_filters,
                                                     n_cohort_samples,
                                                     n_cohort_mut_total,
                                                     n_cohort_mut_snv,
                                                     n_cohort_mut_mnv,
                                                     n_cohort_mut_ins,
                                                     n_cohort_mut_del,
                                                     mut_samples,
                                                     mut_samples_to_alt,
                                                  ]))

                        # Write
                        if self.remove_nonannotated_hotspots is True:
                            if genomic_elements_full != 'None':
                                ofd.write('{}\n'.format('\t'.join(data_to_write)))
                                hotspot_count += 1
                        else:
                            ofd.write('{}\n'.format('\t'.join(data_to_write)))
                            hotspot_count += 1

        logger.info(f'Total hotspots identified: {hotspot_count}')

    def run(self):
        """
        Run HotspotFinder analysis

        Returns:
            None

        """

        # Parse mutations and apply warning checks for samples with more than one mutation in the same chr:position
        self.parse_mutations()
        logger.info('Mutations parsed')

        # Load mappability data into IntervalTree
        if self.mappable_regions_file == 'bgdata':
            pass
            #FIXME
        else:
            self.mappable_regions_tree = self.load_mappability(self.mappable_regions_file)
        if self.blacklisted_regions_file == 'bgdata':
            pass
            # FIXME
        else:
            self.blacklisted_regions_tree = self.load_mappability(self.blacklisted_regions_file)
        logger.info('Mappability data loaded')

        # Load genomic elements into IntervalTree
        regions_data = []
        if os.path.isfile(self.genomic_elements):
            regions_data += [('regions', self.genomic_elements)]
        elif self.genomic_elements == 'all':
            for genomic_element in self.genomic_elements_list:
                regions_data += [(genomic_element, bgdata.get(f'genomicregions/{self.genome}/{genomic_element}'))]
        else:
            regions_data += [(self.genomic_elements,
                              bgdata.get(f'genomicregions/{self.genome}/{self.genomic_elements}'))]
        self.regions_tree = self.load_genomic_elements(regions_data)
        logger.info('Genomic regions loaded')

        # Identify hotspots based on the threshold of mutations
        self.find_hotspots()
        logger.info('Hotspots identified')

        # Write info
        # TODO implement VCF output
        self.write_hotspots()
        logger.info('Hotspots saved to file')


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input-mutations', default=None, required=True, type=click.Path(exists=True),
              help='User input file containing somatic mutations in TSV format')
@click.option('-o', '--output-directory', default=None, required=True, help='Output directory')
@click.option('-g', '--genome', default=None, type=click.Choice(['hg38']), help='Genome to use')
@click.option('-cutoff', '--mutations-cutoff', type=click.IntRange(min=2, max=None, clamp=False), default=None,
              help='Cutoff of mutations to define a hotspot. Default is 3')
@click.option('-group', '--group-by', default=None,
              type=click.Choice(['GROUP', 'GROUP_BY', 'COHORT', 'CANCER_TYPE', 'PLATFORM']),
              help='Header of the column to group hotspots identification')
@click.option('-c', '--cores', type=click.IntRange(min=1, max=os.cpu_count(), clamp=False), default=None,
              help='Number of cores to use in the computation. By default it uses all available cores.')
@click.option('-conf', '--configuration-file', default='./hotspotfinder_v0.1.0.conf',
              required=False, type=click.Path(exists=True), help='User input configuration file')
@click.option('--log-level', default='info', help='Verbosity of the logger',
              type=click.Choice(['debug', 'info', 'warning', 'error', 'critical']))
def main(
        input_mutations,
        output_directory,
        mutations_cutoff,
        genome,
        group_by,
        cores,
        configuration_file,
        log_level
):
    """
    HotspotFinder looks for hotspots of somatic mutations across the genome

    Args:
        input_mutations (str): user path to input mutation data, required
        output_directory (str): path to output directory, required
        mutations_cutoff (int): cutoff of mutations to define a hotspot
        genome (str): reference genome
        group_by (str): name of the column to group hotspots identification
        cores (int): number of cpu
        configuration_file (str): path to configuration file
        log_level (str): logger verbosity

    Returns:
        None
    """
    global logger

    # File name
    file_name = input_mutations.split('/')[-1].split('.')[0]

    # Output directories
    if not os.path.exists(output_directory):
        os.makedirs(output_directory, exist_ok=True)

    # Log
    daiquiri.setup(level=LOGS[log_level], outputs=(
        daiquiri.output.STDERR,
        daiquiri.output.File(
            filename=f'{file_name}.log',
            directory=output_directory
        )
    ))
    logger = daiquiri.getLogger()
    logger.info('HotspotFinder')
    logger.info('Initializing HotspotFinder...')

    # Read configuration file
    config = load_configuration(config_file=configuration_file)

    # Use configuration file parameters when they are not provided by the command line
    if not genome:
        genome = config['genome']['build']
    if not mutations_cutoff:
        mutations_cutoff = config['hotspot_mutations']['cutoff']
    if not group_by:
        group_by = config['group']['groupby']
    if not cores:
        cores = config['settings']['cores']
        if cores > os.cpu_count():
            cores = os.cpu_count()
            logger.warning('Input cores greater than maximum CPU cores. Using maximum CPU cores instead')

    # Mappability
    if os.path.isfile(config['mappability']['mappable_regions']):
        mappable_regions = config['mappability']['mappable_regions']
    else:
        if config['mappability']['mappable_regions'] == 'bgdata':
            # TODO add bgdata
            logger.error(f"bgdata mappable regions file is not available")
            sys.exit(-1)
        else:
            logger.error(f"Mappable regions file does not exist: {config['mappability']['mappable_regions']}")
            sys.exit(-1)
    if os.path.isfile(config['mappability']['blacklisted_regions']):
        blacklisted_regions = config['mappability']['blacklisted_regions']
    else:
        if config['mappability']['blacklisted_regions'] == 'bgdata':
            # TODO add bgdata
            logger.error(f"bgdata blacklisted regions file is not available")
            sys.exit(-1)
        else:
            logger.error(f"Blacklisted regions file does not exist: {config['mappability']['blacklisted_regions']}")
            sys.exit(-1)

    # Population variants
    if os.path.isdir(config['polymorphisms']['population_variants']):
        population_variants = config['polymorphisms']['population_variants']
    else:
        if config['polymorphisms']['population_variants'] == 'bgdata':
            # TODO add bgdata
            logger.error(f"bgdata population variants info is not available")
            sys.exit(-1)
        else:
            logger.error(f"Population variants file does not exist: {config['polymorphisms']['population_variants']}")
            sys.exit(-1)

    # Genomic elements
    if os.path.isfile(config['genomic_regions']['genomic_elements']):
        genomic_elements = config['genomic_regions']['genomic_elements']
    else:
        if config['genomic_regions']['genomic_elements'] in {'all', 'cds', '5utr', '3utr', 'proximal_promoters',
                                                             'distal_promoters', 'introns'}:
            genomic_elements = config['genomic_regions']['genomic_elements']
        else:
            logger.error(f"Genomic regions file does not exist: {config['genomic_regions']['genomic_elements']}")
            sys.exit(-1)

    # Extra parameters
    remove_unknown_nucleotides = config['reference_nucleotides']['remove_unknowns']
    remove_nonannotated_hotspots = config['genomic_regions']['remove_nonannotated_hotspots']
    split_alternates = config['alternates']['split']
    output_format = config['settings']['output_format']
    is_gzip = config['settings']['gzip']

    # Output file names
    compression = '.gz' if is_gzip else ''
    output_file_results = os.path.join(output_directory, f'{file_name}.results.txt{compression}')
    output_file_warning = os.path.join(output_directory, f'{file_name}.warningpositions.txt{compression}')

    logger.info('Analysis parameters loaded')

    logger.info('\n'.join([
        '\n'
        f'* Input mutations file: {input_mutations}',
        f'* Mapable regions file: {mappable_regions}',
        f'* Blacklisted regions file: {blacklisted_regions}',
        f'* Population variants directory: {population_variants}',
        f'* Genomic elements file: {genomic_elements}',
        f'* Output results directory: {output_directory}',
        f'* Output format: {output_format}',
        f'* GZIP compression: {is_gzip}',
        f'* Cutoff hotspots mutations: {mutations_cutoff}',
        f'* Hotspots split alternates: {split_alternates}',
        f'* Remove unknown nucleotides: {remove_unknown_nucleotides}',
        f'* Remove nonannotated hotspots: {remove_nonannotated_hotspots}',
        f'* Genome: {genome}',
        f'* Group analysis by: {group_by}',
        f'* Cores: {cores}',
    ]))

    # Suppress log messages from some libraries
    daiquiri.getLogger('bgdata').setLevel(logging.WARNING)

    # Generate stage 1 hotspots (no annotations)
    experiment = HotspotFinder(
        input_mutations,
        mappable_regions,
        blacklisted_regions,
        population_variants,
        genomic_elements,
        output_file_results,
        output_file_warning,
        output_format,
        is_gzip,
        mutations_cutoff,
        split_alternates,
        remove_unknown_nucleotides,
        remove_nonannotated_hotspots,
        genome,
        '' if group_by == 'none' else group_by,
        cores
    )

    experiment.run()
    logger.info('HotspotFinder analysis finished!')


if __name__ == '__main__':
    main()
