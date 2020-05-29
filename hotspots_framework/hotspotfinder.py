"""
Hotspot Finder identifies mutational hotspots and generates basic annotations
"""

# Import modules
from collections import defaultdict
import gzip
import json
import logging
import os
import shutil

import bgreference as bgref
from bgparsers import readers
import click
import daiquiri
from intervaltree import IntervalTree
import tabix

# Global variables
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
LOGS = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}


class HotspotFinder:
    """Class to identify and annotate hotspots with basic information"""

    def __init__(self, input_file, output_file_results, output_file_warning, output_path_tmp, hotspot_mutations, genome, group_by, cores):
        """
        Initialize HotspotFinder class
        Args:
            input_file (str): path to input data
            output_file_results (str): path to output file
            output_file_warning (str): path to genomic positions where warning is found
            output_path_tmp (str): path to tmp directory
            hotspot_mutations (int): cutoff of mutations to define a hotspot
            genome (str): reference genome
            group_by (str): name of the column to group hotspots identification
            cores (int): number of cpu

        Returns:
            None
        """

        self.input_file = input_file
        self.output_file_hotspots = output_file_results
        self.output_file_warning = output_file_warning
        # self.output_file_tmp_hotspots = os.path.join(output_path_tmp, 'hotspots_v1.txt')
        # self.output_file_tmp_warning = os.path.join(output_path_tmp, 'warning_mutations.txt')
        self.hotspot_mutations = hotspot_mutations
        self.genome = genome
        self.group_by = group_by
        self.cores = cores

        # TODO add command line
        merge_mutations_types = False
        if not merge_mutations_types:
            self.muttype_dict = {
                's': 'snv',
                'm': 'mnv',
                'i': 'ins',
                'd': 'del'
            }
        else:
            self.muttype_dict = {
                's': 'mut',
                'm': 'mut',
                'i': 'mut',
                'd': 'mut'
            }

        self.cohort_to_sample = defaultdict(set)
        self.cohort_to_mutation_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        self.cohort_to_mutation_alts = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.hotspots = defaultdict(lambda: defaultdict(dict))

        #TODO refactor mutations as namedtuple so that we can check the REF in file with bgreference

    @staticmethod
    def json_convert_save(dictionary, output_file):
        """
        Convert dictionary to json format and save it

        Args:
            dictionary (dict): dictionary of hotspots
            output_file (str): path to output json file

        Returns:
            None
        """
        json_dict = json.dumps(dict(dictionary))
        with open(output_file, 'w') as ofd:
            ofd.write(json_dict)

    @staticmethod
    def tabix_convert_save(file):
        """
        Convert to tabix and save
        Args:
            file (str): path to file

        Returns:
            None
        """

        file_gz = '{}.gz'.format(file)
        os.system('cat {} | sort -k1V -k2n -k3n | bgzip > {}'.format(file, file_gz))
        os.system('tabix --sequence 1 --begin 2 --end 3 {}'.format(file_gz))
        os.system('rm {}'.format(file))

    @staticmethod
    def load_intervaltree(files):
        """
        Load regions into intervaltree
        Args:
            files (list): list of path(s) to file(s)

        Returns:
            tree (dict): tree with regions, keys are chromosomes, data are regions
        """

        tree = defaultdict(IntervalTree)
        for file in files:
            with gzip.open(file, 'rt') as fd:
                for line in fd:
                    chrom, start, end, strand, info = line.strip().split('\t')
                    tree[chrom].addi(int(start), int(end) + 1, info)  # +1 interval
        return tree

    def parse_mutations(self):
        """
        Load mutations file into dictionaries containing number of mutations and alternates

        Args:
            None

        Returns:
            None

        """

        chromosomes = list(map(str, range(1, 23))) + ['X', 'Y']
        cohort_nmuts = defaultdict(lambda: defaultdict(int))
        substitutions_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        insertions_dict = defaultdict(lambda: defaultdict(list))
        deletions_dict = defaultdict(lambda: defaultdict(list))

        # Read mutations
        for row in readers.variants(
                file=self.input_file,
                required=['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE'],
                extra=['GROUP', 'GROUP_BY', 'COHORT', 'CANCER_TYPE']
        ):
            chromosome = row['CHROMOSOME']
            chr_position = str(row['POSITION'])
            ref = row['REF']
            alt = row['ALT']
            sample = row['SAMPLE']
            chr_pos = '{}_{}'.format(chromosome, chr_position)
            # Identify group
            # If no group, hotspots are computed using the whole input file
            if not self.group_by:
                cohort = 'cohort'
            else:
                cohort = row[self.group_by]
            # Keep track of samples from each group
            self.cohort_to_sample[cohort].add(sample)
            # Keep track of how many mutations each group has
            cohort_nmuts['raw'][cohort] += 1

            # Read mutations in autosomal + sexual chromosomes
            if chromosome in set(chromosomes):
                if ref != alt:
                    # Read substitutions of any length
                    if ref != '-' and alt != '-':
                        if len(alt) == 1:
                            substitutions_dict['snvs'][sample][chr_pos].append(alt)
                        else:
                            substitutions_dict['mnvs'][sample][chr_pos].append(alt)
                    # Read indels of any length
                    else:
                        # Insertions
                        if ref == '-':
                            insertions_dict[sample][chr_pos].append(alt)
                        # Deletion
                        elif alt == '-':
                            deletions_dict[sample][chr_pos].append(alt)

        for cohort, nmuts in cohort_nmuts['raw'].items():
            logger.info(f'Input mutations in {cohort} = {nmuts}')

        # Check mutations per sample and add to cohort_to_mutation dicts
        # Write warning positions
        header = ['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE', 'WARNING', 'SKIP']
        with open(self.output_file_warning, 'w') as ofd:
            ofd.write('{}\n'.format('\t'.join(header)))
            for cohort, set_of_samples in self.cohort_to_sample.items():
                for sample in set_of_samples:
                    snvs_in_sample = substitutions_dict['snvs'][sample]
                    mnvs_in_sample = substitutions_dict['mnvs'][sample]
                    ins_in_sample = insertions_dict[sample]
                    dels_in_sample = deletions_dict[sample]

                    total_muts_in_sample = defaultdict(list)
                    for set_of_muts, muttype in [
                        (snvs_in_sample, self.muttype_dict['s']),
                        (mnvs_in_sample, self.muttype_dict['m']),
                        (ins_in_sample, self.muttype_dict['i']),
                        (dels_in_sample, self.muttype_dict['d'])
                    ]:
                        for chr_pos, alts in set_of_muts.items():
                            total_muts_in_sample[chr_pos] += [(a, muttype) for a in alts]
                    warning_mutations = [(k, v) for k, v in total_muts_in_sample.items() if len(v) > 1]

                    logger.debug(sample)
                    logger.debug(total_muts_in_sample)
                    logger.debug(warning_mutations)

                    self.warning_chr_position = set()
                    if warning_mutations:
                        # Load warning positions
                        for chr_position, alts in warning_mutations:
                            chromosome, position = chr_position.split('_')
                            alts_unique = set(alts)
                            alts_simplified = [a for a, muttype in alts]
                            self.warning_chr_position.add(chr_position)
                            if len(alts_unique) == 1:
                                logger.warning(
                                    f'Sample "{sample}" position chr{chr_position} has 2 alternates: {alts_simplified}')
                                alt, muttype = list(alts_unique)[0]
                                self.cohort_to_mutation_alts[cohort][muttype][chr_position] += [alt]
                                self.cohort_to_mutation_counts[cohort][muttype][chr_position] += 1
                                #FIXME add reference when mutations are namedtuple
                                ofd.write('{}\n'.format('\t'.join([
                                    chromosome, position, 'NA', ','.join(alts_simplified), sample, 'warning_1', 'False'
                                ])))
                            elif len(alts_unique) == 2:
                                logger.warning(
                                    f'Sample "{sample}" position chr{chr_position} has 2 or 3 alternates: {alts_simplified}. '
                                    f'Mutation counts and alternates might not match.')
                                for mutation in alts_unique:
                                    alt, muttype = mutation
                                    self.cohort_to_mutation_alts[cohort][muttype][chr_position] += [alt]
                                    self.cohort_to_mutation_counts[cohort][muttype][chr_position] += 1
                                    ofd.write('{}\n'.format('\t'.join([
                                        chromosome, position, 'NA', ','.join(alts_simplified), sample, 'warning_2', 'False'
                                    ])))
                            else:
                                logger.warning(
                                    f'Sample "{sample}" position chr{chr_position} has 3 or more different alternates: '
                                    f'{alts_simplified}. Mutations are skipped from analysis')
                                ofd.write('{}\n'.format('\t'.join([
                                    chromosome, position, 'NA', ','.join(alts_simplified), sample, 'warning_3', 'True'
                                ])))

                    # Load non-warning positions
                    for chr_position, mutation in total_muts_in_sample.items():
                        if not chr_position in self.warning_chr_position:
                            self.cohort_to_mutation_alts[cohort][mutation[0][1]][chr_position] += [mutation[0][0]]
                            self.cohort_to_mutation_counts[cohort][mutation[0][1]][chr_position] += 1

        # logger.debug(self.cohort_to_mutation_counts)
        logger.debug(self.cohort_to_mutation_alts)

    def find_hotspots(self):
        """
        Identifies hotspots and generates raw hotspots file

        Returns:
            None
        """
        for cohort, data in self.cohort_to_mutation_counts.items():
            for muttype, mutated_positions in data.items():
                self.hotspots[cohort][muttype] = {k: v for k, v in mutated_positions.items() if v >= self.hotspot_mutations}

        logger.debug(self.hotspots)

    def write_raw_hotspots(self):
        """
        Writes output file for HotspotFinder

        Returns:
            None
        """
        nucleotides = ('A', 'C', 'G', 'T')
        header = [
            'CHROMOSOME',
            'POSITION',
            'ID',
            'COHORT',
            'N_COHORT_SAMPLES',
            'N_MUT_SAMPLES',
            'FRAC_MUT_SAMPLES',
            'MUT_TYPE',
            'N_MUT_TYPE',
            'REF',
            'ALT',
            'FRAC_ALT',
            'CONTEXT_3',
            'CONTEXT_5',
            'WARNING_POSITION'
        ]

        # FIXME how does it work with merged?
        with open(self.output_file_hotspots, 'w') as ofd:
            ofd.write('{}\n'.format('\t'.join(header)))
            for cohort, data in self.hotspots.items():
                n_cohort_samples = len(self.cohort_to_sample[cohort])
                for muttype, hotspots in data.items():
                    for hotspot_id, n_mut_samples in sorted(hotspots.items(), key=lambda item: item[1], reverse=True):
                        chromosome, position = hotspot_id.split('_')
                        frac_mut_samples = str(n_mut_samples / n_cohort_samples)

                        # Sequence info
                        pentamer_sequence = bgref.refseq(self.genome, chromosome, int(position) - 2, 5)
                        trimer_sequence = pentamer_sequence[1:4]
                        ref = pentamer_sequence[2]

                        # Alternates
                        list_of_alternates = self.cohort_to_mutation_alts[cohort][muttype][hotspot_id]
                        alternates_counts = []
                        alternates_fractions = []
                        for nucleotide in nucleotides:
                            count = list_of_alternates.count(nucleotide)
                            alternates_counts.append(f'{nucleotide}={count}')
                            alternates_fractions.append(f'{nucleotide}={count / n_mut_samples}')
                        alternates_counts = ','.join(alternates_counts)
                        alternates_fractions = ','.join(alternates_fractions)

                        # Warning position
                        # FIXME rename
                        warning_position = 'True' if hotspot_id in self.warning_chr_position else 'False'

                        # Write
                        ofd.write('{}\n'.format('\t'.join(list(map(str,
                            [
                                chromosome,
                                position,
                                hotspot_id,
                                cohort,
                                n_cohort_samples,
                                n_mut_samples,
                                frac_mut_samples,
                                muttype,
                                '1',
                                ref,
                                alternates_counts,
                                alternates_fractions,
                                trimer_sequence,
                                pentamer_sequence,
                                warning_position
                             ]

                        )))))

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input-file', default=None, required=True, type=click.Path(exists=True),
              help='File containing somatic mutations in TSV format')
@click.option('-o', '--output-path', default=None, required=True, help='Output directory')
@click.option('-mut', '--hotspot-mutations', type=click.IntRange(min=2, max=None, clamp=False), default=2,
              help='Cutoff of mutations to define a hotspot. Default is 3')
@click.option('-g', '--genome', default='hg19',
              type=click.Choice(['hg38', 'hg19']),
              help='Genome to use')
@click.option('-group', '--group-by', default=None, type=click.Choice(['GROUP', 'GROUP_BY', 'COHORT', 'CANCER_TYPE']),
              help='Header of the column to group hotspots identification')
@click.option('-c', '--cores', type=click.IntRange(min=1, max=os.cpu_count(), clamp=False), default=os.cpu_count(),
              help='Number of cores to use in the computation. By default it uses all available cores.')
@click.option('--log-level', default='info', help='Verbosity of the logger',
              type=click.Choice(['debug', 'info', 'warning', 'error', 'critical']))
def main(input_file, output_path, hotspot_mutations, genome, group_by, cores, log_level):
    """
    Look for mutational hotspots across the genome
    Args:
        input_file (str): path to input data
        output_path (str): path to output directory
        hotspot_mutations (int): cutoff of mutations to define a hotspot
        genome (str): reference genome
        group_by (str): name of the column to group hotspots identification
        cores (int): number of cpu
        log_level (str): logger verbosity

    Returns:
        None
    """
    global logger

    # Output directories
    if not os.path.exists(output_path):
        os.makedirs(output_path, exist_ok=True)
    output_path_tmp = os.path.join(output_path, 'tmp')
    # if not os.path.exists(output_path_tmp):
    #     os.makedirs(output_path_tmp, exist_ok=True)

    # File name
    file_name = input_file.split('/')[-1].split('.')[0]

    # Output file name
    output_file_results = os.path.join(output_path, f'{file_name}.results.txt')
    output_file_warning = os.path.join(output_path, f'{file_name}.warningpositions.txt')


    # Log
    daiquiri.setup(level=LOGS[log_level], outputs=(
        daiquiri.output.STDERR,
        daiquiri.output.File(
            filename=f'{file_name}.log',
            directory=output_path
        )
    ))
    logger = daiquiri.getLogger()

    logger.info('HotspotFinder')
    logger.info('Initializing HotspotFinder...')
    logger.info('\n'.join([
        '\n'
        f'Input data file: {input_file}',
        f'Output results directory: {output_path}',
        f'Genome: {genome}',
        f'Cutoff hotspots mutations: {hotspot_mutations}'
    ]))

    # Generate stage 1 hotspots (no annotations)
    experiment = HotspotFinder(
        input_file,
        output_file_results,
        output_file_warning,
        output_path_tmp,
        hotspot_mutations,
        genome,
        group_by,
        cores
    )

    experiment.parse_mutations()
    experiment.find_hotspots()
    experiment.write_raw_hotspots()


if __name__ == '__main__':
    main()
