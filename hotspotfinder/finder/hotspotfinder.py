"""
Hotspot Finder identifies mutational hotspots and generates basic annotations
"""

# Import modules
from collections import defaultdict
import gzip
import itertools
import logging

import bgreference as bgref
from bgparsers import readers
from intervaltree import IntervalTree

from hotspotfinder import __logger_name__
from hotspotfinder.finder.counter import MutationCounter
from hotspotfinder.utils import file_open

logger = logging.getLogger(__logger_name__)


class HotspotFinder:
    """Class to identify and annotate hotspots with basic information"""

    GENOMIC_ELEMENTS = [
            'cds',
            'splice_sites',
            '5utr',
            '3utr',
            'proximal_promoters',
            'distal_promoters',
            'introns',
            'lncrna_exons',
            'lncrna_splice_sites',
            'lncrna_proximal_promoters',
            'lncrna_distal_promoters',
            'lncrna_introns',
            'mirna',
            'misc_rna',
            'mt_rrna',
            'mt_trna',
            'rrna',
            'scrna',
            'snorna',
            'snrna',
            'pseudogenes'
        ]
    CHROMOSOMES = list(map(str, range(1, 23))) + ['X', 'Y']

    def __init__(self,
                 input_file,
                 output_file_results,
                 output_file_warning,
                 config
                 ):
        """
        Initialize HotspotFinder class

        Args:
            input_file (str): path to input mutations data
            output_file_results (str): path to output file
            output_file_warning (str): path to genomic positions where warning is found
            config (dict): configuration parameters

        Returns:
            None

        """

        # Input output files
        self.input_file = input_file
        self.output_file_hotspots = output_file_results
        self.output_file_warning = output_file_warning

        # Parameters
        self.genome = config['genome']
        self.cores = config['cores']
        self.hotspot_mutated_samples = config['finder']['samples_cutoff']
        self.hotspot_mutations = config['finder']['mutations_cutoff']
        self.remove_nonannotated_hotspots = config['finder']['remove_nonannotated_hotspots']
        self.group_by = config['finder']['groupby']
        self.split_alternates = config['finder']['split_alternates']
        self.annotate_hotspots = config['finder']['annotate']

        # Initialize variables to load data
        self.mutation_counts = MutationCounter(self.genome)
        self.cohort_to_mutation_alts = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.hotspots = defaultdict(lambda: defaultdict(dict))
        self.warning_chr_position = set()

        # Load genomic annotations if required
        if config['finder']['annotate']:

            logger.info('Loading genomic annotations...')

            # Mappability data
            self.mappable_regions_tree = self.load_mappability(config['mappable_regions'])
            self.blacklisted_regions_tree = self.load_mappability(config['blacklisted_regions'])
            logger.info('Mappability data loaded')

            # Variation data
            self.variation_data_set = self.load_variation(config['population_variants'])
            logger.info('Population variants data loaded')

            # Repeats data
            self.repeats_tree = self.load_repeats(config['repeats'], 'N')
            logger.info('Repeats data loaded')

            # Immunoglobulin and T-cell receptor loci
            self.ig_tr_tree = self.load_ig_tr(config['ig_tr_regions'])
            logger.info('Ig-TR data loaded')

            # Genomic elements data
            self.regions_tree = self.load_genomic_elements(config['genomic_elements'])
            self.genomic_elements_priority = {e: i for i, e in enumerate(HotspotFinder.GENOMIC_ELEMENTS)}
            logger.info('Genomic regions loaded')

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
            for line in fd:    # no header
                chrom, start, end = line.strip().split('\t')
                tree[chrom[trim:]].addi(int(start), int(end) + 1)  # +1 interval
        return tree

    @staticmethod
    def load_variation(file, chr_format='chrN'):
        """
        Load population variants into set
        Args:
            file (path): path to file containing population variants data
            chr_format (str): chromosome format used in file. HotspotFinder works with format '1' instead of 'chr1'

        Returns:
            set_of_interest (set): set with variants annotated as 'N_position'
        """

        set_of_interest = set()
        trim = 3 if chr_format == 'chrN' else None
        with gzip.open(file, 'rt') as fd:
            for line in fd:
                chrom, position = line.strip().split('\t')[:2]
                set_of_interest.add(f'{chrom[trim:]}_{position}')
        return set_of_interest

    @staticmethod
    def load_repeats(file, chr_format='chrN'):
        """
        Load repeats regions into intervaltree
        Args:
            file (path): path to file containing repeats data
            chr_format (str): chromosome format used in file. HotspotFinder works with format '1' instead of 'chr1'

        Returns:
            tree (dict): tree with regions, keys are chromosomes, data are regions
        """

        tree = defaultdict(IntervalTree)
        trim = 3 if chr_format == 'chrN' else None
        with gzip.open(file, 'rt') as fd:
            # Skip header
            next(fd)
            for line in fd:
                chrom, start, end, element = line.strip().split('\t')
                tree[chrom[trim:]].addi(int(start), int(end) + 1, element)  # +1 interval
        return tree

    @staticmethod
    def load_ig_tr(file):
        """
        Load Ig and T-cell receptor regions into intervaltree

        Args:
            file (path): path to file containing data

        Returns:
            tree (dict): tree with regions, keys are chromosomes, data are regions
        """
        tree = defaultdict(IntervalTree)
        with gzip.open(file, 'rt') as fd:
            # Skip header
            next(fd)
            for line in fd:
                chrom, start, end, strand, _, transcript_id, _ = line.strip().split('\t')
                tree[chrom].addi(int(start), int(end) + 1,
                                 f'{transcript_id}::ig_tr')  # +1
        return tree

    @staticmethod
    def load_genomic_elements(file):
        """
        Load genomic element regions into intervaltree
        Args:
            file: file to genomic elements annotations

        Returns:
            tree (dict): tree with regions, keys are chromosomes, data are regions
        """
        tree = defaultdict(IntervalTree)
        with gzip.open(file, 'rt') as fd:
            for line in fd:
                chrom, start, end, strand, gene_id, transcript_id, symbol, genomic_element = line.strip().split('\t')
                tree[chrom].addi(
                    int(start), int(end) + 1, f'{symbol}::{gene_id}::{transcript_id}::{genomic_element}')  # +1
        return tree

    def parse_mutations(self):
        """
        Load mutations file into dictionaries containing number of mutations and alternates

        Returns:
            None

        """

        # Cohort name
        # By default, hotspots are computed using all mutations in the input file
        # This is overwritten when self.group_by is True
        cohort = self.input_file.split('/')[-1].split('.')[0]

        # Read mutations
        for row in readers.variants(
                file=self.input_file,
                required=['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE'],
                extra=['GROUP', 'GROUP_BY', 'COHORT', 'CANCER_TYPE', 'PLATFORM']
        ):
            chromosome = row['CHROMOSOME']
            position = row['POSITION']
            ref = row['REF']
            alt = row['ALT']
            alt_type = row['ALT_TYPE']
            sample = row['SAMPLE']

            # Identify group
            # If analysis is carried out per group, use the group value of the row
            if self.group_by:
                cohort = row[self.group_by]

            # Read mutations in autosomal + sexual chromosomes
            # Skip mutations whose reference nucleotide(s) is/are equal to the alternate nucleotide(s)
            if chromosome in set(HotspotFinder.CHROMOSOMES) and ref != alt:
                self.mutation_counts.add_mutation(chromosome, position, ref, alt, alt_type, sample, cohort)

        if self.mutation_counts.reference_mismatch > 0:
            logger.warning(f'A total of {self.mutation_counts.reference_mismatch} SNV/MNV/deletion REF nucleotides '
                           f'do not match the reference. '
                           f'These mutations are discarded from the analysis')
        if self.mutation_counts.unknown_nucleotides_in_mutation > 0:
            logger.warning(
                f'A total of {self.mutation_counts.unknown_nucleotides_in_mutation} mutations '
                f'contain unknown nucleotides. '
                f'These mutations are discarded from the analysis')
        if self.mutation_counts.unknown_nucleotides_in_context > 0:
            logger.warning(
                f'A total of {self.mutation_counts.unknown_nucleotides_in_context} mutations '
                f'contain unknown nucleotides in their pentamer nucleotide context. '
                f'These mutations are discarded from the analysis')

        # Just logging information
        for cohort in self.mutation_counts.get_cohorts():
            muts_cohort = 0
            for muttype in ('snv', 'mnv', 'ins', 'del'):
                nmuts = self.mutation_counts.n_mutations_cohort(cohort, muttype)
                logger.info(f'Input {muttype} mutations in {cohort} = {nmuts}')
                muts_cohort += nmuts
            logger.info(f'Input total mutations in {cohort} = {muts_cohort}')

    def load_mutations(self):
        """
        Check that there is one mutation per sample and add it to a dictionary to compute hotspots
        If a sample has more than one mutation in a given position, report it and filter mutations out if necessary

        This function writes the warning positions output file
        """

        header = ['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE', 'WARNING', 'SKIP']
        with file_open(self.output_file_warning) as ofd:
            ofd.write('{}\n'.format('\t'.join(header)))
            warning_samples_n = 0
            for cohort in self.mutation_counts.get_cohorts():
                for sample in self.mutation_counts.get_samples(cohort):

                    # Get mutations per sample
                    total_muts_in_sample = defaultdict(list)
                    for chr_pos, muttype, alts in self.mutation_counts.get_mutations_and_alternates(cohort, sample):
                        total_muts_in_sample[chr_pos] += [(a, muttype) for a in alts]

                    # Check number of mutations per position
                    for chr_position, alts_data in total_muts_in_sample.items():

                        # There is only one mutation per position (this is expected)
                        if len(alts_data) == 1:
                            alt, muttype = alts_data[0]
                            self.cohort_to_mutation_alts[cohort][muttype][f'{chr_position}_{muttype}'] += [alt]

                        # There is more than one mutation per position
                        elif len(alts_data) > 1:
                            warning_samples_n += 1
                            chromosome, position = chr_position.split('_')
                            alts_unique = set(alts_data)
                            alts_simplified = [a for a, muttype in alts_data]
                            self.warning_chr_position.add(chr_position)

                            # The mutations are the same
                            if len(alts_unique) == 1:
                                logger.debug(
                                    f'Sample "{sample}" position chr{chr_position} has repeated alternates: '
                                    f'{alts_simplified}. All alternates are kept, but '
                                    f'the number of mutations and the number of mutated samples will not match.'
                                )
                                alt, muttype = list(alts_unique)[0]
                                self.cohort_to_mutation_alts[cohort][muttype][
                                    f'{chr_position}_{muttype}'] += alts_simplified
                                reference_n = bgref.refseq(self.genome, chromosome, int(position), 1)
                                ofd.write('{}\n'.format('\t'.join([
                                    chromosome, position, reference_n, ','.join(alts_simplified),
                                    sample, 'warning_1', 'False'
                                ])))
                            # There are two different mutations
                            elif len(alts_unique) == 2:
                                logger.debug(
                                    f'Sample "{sample}" position chr{chr_position} has 2 different alternates: '
                                    f'{alts_simplified}. All alternates are kept, but '
                                    f'the number of mutations and the number of mutated samples will not match.'
                                )
                                reference_n = bgref.refseq(self.genome, chromosome, int(position), 1)
                                ofd.write('{}\n'.format('\t'.join([
                                    chromosome, position, reference_n, ','.join(alts_simplified),
                                    sample, 'warning_2', 'False'
                                ])))
                                for mutation in alts_data:
                                    alt, muttype = mutation
                                    self.cohort_to_mutation_alts[cohort][muttype][f'{chr_position}_{muttype}'] += [alt]
                            # There are 3 or more different mutations
                            else:
                                logger.debug(
                                    f'Sample "{sample}" position chr{chr_position} has 3 or more different alternates: '
                                    f'{alts_simplified}. Mutations are skipped from analysis')
                                reference_n = bgref.refseq(self.genome, chromosome, int(position), 1)
                                ofd.write('{}\n'.format('\t'.join([
                                    chromosome, position, reference_n, ','.join(alts_simplified),
                                    sample, 'warning_3', 'True'
                                ])))
                                for mutation in alts_data:
                                    alt, muttype = mutation
                                    self.mutation_counts.discard_mutation(cohort, sample, muttype, chr_position)

        if warning_samples_n > 0:
            logger.warning(f'A total of {warning_samples_n} sample{"s" if warning_samples_n > 1 else ""} '
                           f'contain{"" if warning_samples_n > 1 else "s"} mutations flagged with warnings. '
                           f'Please check: {self.output_file_warning}')

    def find_hotspots(self):
        """
        Function to identify hotspots: genomic positions with number of mutations >= hotspot_mutations
        Further filters are applied in write_hotspots() function
        """

        # Identify hotspots as genomic positions with number of mutations >= hotspot_mutations
        for cohort, data in self.cohort_to_mutation_alts.items():
            for muttype, mutated_positions in data.items():
                if self.split_alternates:
                    for hotspot_id, list_of_alternates in mutated_positions.items():
                        for alternate in set(list_of_alternates):
                            count = list_of_alternates.count(alternate)
                            if count >= self.hotspot_mutations:
                                self.hotspots[cohort][muttype][f'{hotspot_id}>{alternate}'] = [alternate] * count
                else:
                    self.hotspots[cohort][muttype] = {k: v for k, v in mutated_positions.items() if
                                                      len(v) >= self.hotspot_mutations}

    def write_hotspots(self):
        """
        Writes output file for HotspotFinder

        Returns:
            None
        """
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
            'N_COHORT_SAMPLES',
            'N_COHORT_MUTATIONS_TOTAL',
            'N_COHORT_MUTATIONS_SNV',
            'N_COHORT_MUTATIONS_MNV',
            'N_COHORT_MUTATIONS_INS',
            'N_COHORT_MUTATIONS_DEL',
            'MUTATED_SAMPLES',
            'MUTATED_SAMPLES_ALTS',
            'OVERLAP_WARNING_POSITION'
        ]
        if self.annotate_hotspots:
            header += [
                'GENOMIC_ELEMENT',
                'SYMBOL',
                'GENE_ID',
                'TRANSCRIPT_ID',
                'GENOMIC_REGION',
                'GENOMIC_REGION_PRIORITY',
                'CODING_NONCODING',
                'MAPPABILITY',
                'MAPPABILITY_BLACKLIST',
                'VARIATION_AF',
                'HOTSPOTFINDER_FILTERS',
                'REPEATS',
                'REPEATS_OVERLAP',
                'IG_TR',
                'IG_TR_OVERLAP',
            ]

        hotspot_count = 0
        with file_open(self.output_file_hotspots) as ofd:
            ofd.write('{}\n'.format('\t'.join(header)))

            # Cohort data
            for cohort, data in self.hotspots.items():
                n_cohort_samples = len(self.mutation_counts.get_samples(cohort))
                n_cohort_mut_snv = self.mutation_counts.n_mutations_cohort(cohort, 'snv')
                n_cohort_mut_mnv = self.mutation_counts.n_mutations_cohort(cohort, 'mnv')
                n_cohort_mut_ins = self.mutation_counts.n_mutations_cohort(cohort, 'ins')
                n_cohort_mut_del = self.mutation_counts.n_mutations_cohort(cohort, 'del')
                n_cohort_mut_total = n_cohort_mut_snv + n_cohort_mut_mnv + \
                                     n_cohort_mut_ins + n_cohort_mut_del
                # Read hotspots
                for muttype, hotspots in data.items():

                    for hotspot_id_type, list_of_alternates in sorted(
                            hotspots.items(), key=lambda item: len(item[1]), reverse=True):

                        data_to_write = []

                        # Chromosome and position
                        # Remove alternate if self.split_alternates option is True
                        # Works also in default option self.split_alternates False
                        chromosome, position, _ = hotspot_id_type.split('>')[0].split('_')
                        chr_pos = f'{chromosome}_{position}'

                        # Number of mutations
                        n_mutations = len(list_of_alternates)

                        # Mutated samples and alternates
                        # Get mutated samples and their alternates from counter module
                        mut_samples = self.mutation_counts.get_samples_per_mutation(chr_pos, muttype, cohort)
                        mut_samples_to_alt = self.mutation_counts.get_alternates_per_mutation(chr_pos, muttype, cohort)

                        # If hotspots are alternate specific (split alternates is True)
                        # check which samples have the alternate of the hotspot and update variables
                        if self.split_alternates:
                            alternate = hotspot_id_type.split('>')[1]
                            mut_samples = set()
                            for sample, sample_alternates in mut_samples_to_alt.items():
                                for unique_alternate in sample_alternates:
                                    if unique_alternate == alternate:
                                        mut_samples.add(sample)
                            mut_samples_to_alt = {k: v for k, v in mut_samples_to_alt.items() if k in mut_samples}
                            list_of_alternates = list(itertools.chain.from_iterable(list(mut_samples_to_alt.values())))
                        # Number of mutated samples
                        n_mut_samples = len(mut_samples)

                        # Skip hotspots below threshold of number of mutated samples
                        if n_mut_samples < self.hotspot_mutated_samples:
                            logger.debug('Number of mutated samples is below the hotspot definition threshold:\n'
                                         f'{hotspot_id_type} has {n_mut_samples}'
                                         f' mutated sample{"s" if n_mut_samples > 1 else ""}')
                            continue

                        # Alternates
                        # Reformat alternates and add other information to write in output file
                        if muttype != 'del':
                            alternates = ','.join(sorted(list(set(list_of_alternates))))
                        else:
                            alternates = '-'
                        alternates_counts = [f'{a}={list_of_alternates.count(a)}' for a in set(list_of_alternates)]
                        alternates_counts = ','.join(alternates_counts)
                        alternates_fractions = [f'{a}={list_of_alternates.count(a) / n_mutations}' for a in
                                               set(list_of_alternates)]
                        alternates_fractions = ','.join(alternates_fractions)
                        frac_mut_samples = str(n_mut_samples / n_cohort_samples)
                        mut_samples_to_alt = ';'.join([f'{k}::{",".join(v)}' for k, v in mut_samples_to_alt.items()])
                        mut_samples = ';'.join(list(mut_samples))

                        # Sequence context and reference nucleotide(s)
                        # Retrieve reference genome sequence contexts
                        pentamer_sequence = bgref.refseq(self.genome, chromosome, int(position) - 2, 5)
                        trimer_sequence = pentamer_sequence[1:4]
                        # Retrieve reference nucleotides
                        if muttype == 'snv':
                            ref = pentamer_sequence[2]
                        elif muttype == 'ins':
                            ref = '-'
                        else:  # MNVs and deletions
                            ref = set()
                            for sample in mut_samples.split(';'):
                                ref.update(self.mutation_counts.reference(chr_pos, muttype, sample))
                            ref = ','.join(list(ref))

                        logger.debug(list(map(str, [
                                                  hotspot_id_type,
                                                  muttype,
                                                  cohort,
                                                  n_mutations,
                                                  n_mut_samples,
                                                  ref,
                                                  alternates,
                                                  alternates_counts,
                                                  alternates_fractions,
                                                  mut_samples_to_alt
                                              ])))

                        # Warning flag: at least one sample in the dataset contains a warning in this position
                        warning_flag = 'True' if chr_pos in self.warning_chr_position else 'False'

                        # Append hotspot data to list
                        # When the method runs without annotations, this is the output per hotspot
                        data_to_write += list(map(str,
                                                 [   chromosome,
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
                                                     n_cohort_samples,
                                                     n_cohort_mut_total,
                                                     n_cohort_mut_snv,
                                                     n_cohort_mut_mnv,
                                                     n_cohort_mut_ins,
                                                     n_cohort_mut_del,
                                                     mut_samples,
                                                     mut_samples_to_alt,
                                                     warning_flag]))

                        # Genomic annotations
                        # If the method runs with annotations, compute them and add them to results
                        if self.annotate_hotspots:

                            # HotspotFinder flags
                            # The hotspot is checked against 3 main flags:
                            # mappability, mappability blacklisted regions and variation (SNP overlap)
                            hotspotfinder_filters = 2

                            # Mappability
                            map_data = 'low'
                            for _ in self.mappable_regions_tree[chromosome][int(position)]:
                                map_data = 'high'
                                hotspotfinder_filters += 1
                                break

                            # Mappability blacklisted regions
                            blacklist_data = 'PASS'
                            for _ in self.blacklisted_regions_tree[chromosome][int(position)]:
                                blacklist_data = 'FAIL'
                                hotspotfinder_filters -= 1
                                break

                            # Variation
                            var_data = 'PASS'
                            if {chr_pos}.intersection(self.variation_data_set):
                                var_data = 'FAIL'
                                hotspotfinder_filters -= 1

                            hotspotfinder_filters = 'PASS' if hotspotfinder_filters == 3 else 'FAIL'

                            # Additional annotations
                            # Repeats
                            repeats = set()
                            for intersect in self.repeats_tree[chromosome][int(position)]:
                                if intersect:
                                    repeats.add(intersect.data)
                            in_repeat = 'True' if repeats else 'False'
                            repeats = ','.join(sorted(list(repeats))) if repeats else 'None'

                            # Ig and TR
                            ig_tr_loci = set()
                            for intersect in self.ig_tr_tree[chromosome][int(position)]:
                                if intersect:
                                    ig_tr_loci.add(intersect.data)
                            in_ig_tr_loci = 'True' if ig_tr_loci else 'False'
                            ig_tr_loci = ','.join(sorted(list(ig_tr_loci))) if ig_tr_loci else 'None'

                            # Genomic elements
                            genomic_elements_full = []
                            symbol = set()
                            geneid = set()
                            transcriptid = set()
                            genomic_elements_type = set()
                            for intersect in self.regions_tree[chromosome][int(position)]:
                                if intersect:
                                    genomic_elements_full.append(intersect.data)
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

                            # If specified (remove_nonannotated_hotspots is True)
                            # filter out hotspots that do not overlap any known region
                            if self.remove_nonannotated_hotspots is True and genomic_elements_full == 'None':
                                continue

                            # Coding or noncoding flag
                            coding_noncoding = 'CODING' if 'cds' in genomic_elements_type else 'NONCODING'

                            data_to_write += list(map(str,
                                                 [   genomic_elements_full,
                                                     symbol,
                                                     geneid,
                                                     transcriptid,
                                                     genomic_elements_type,
                                                     genomic_elements_type_priority,
                                                     coding_noncoding,
                                                     map_data,
                                                     blacklist_data,
                                                     var_data,
                                                     hotspotfinder_filters,
                                                     repeats,
                                                     in_repeat,
                                                     ig_tr_loci,
                                                     in_ig_tr_loci]))

                        # Write
                        ofd.write('{}\n'.format('\t'.join(data_to_write)))
                        hotspot_count += 1

        logger.info(f'Total hotspots identified: {hotspot_count}')

    def run(self):
        """
        Run HotspotFinder analysis

        Returns:
            None

        """

        logger.info('Parsing mutations...')
        # Parse mutations and apply warning checks for samples with more than one mutation in the same chr:position
        self.parse_mutations()
        self.load_mutations()
        logger.info('Mutations parsed')

        # Identify hotspots based on the threshold of mutations
        logger.info('Identifying hotspots...')
        self.find_hotspots()
        logger.info('Hotspots identified')

        # Write info
        logger.info('Writing output file...')
        self.write_hotspots()
        logger.info('Hotspots saved to file')
