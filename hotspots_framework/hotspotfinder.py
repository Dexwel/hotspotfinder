"""
Hotspot Finder identifies hotspots and generates basic annotations
"""

# Import modules
from collections import defaultdict
import gzip
import json
import logging
import os
import shutil

from bgparsers import readers
import bgreference as bgref
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

    def __init__(self, input_file, output_file, output_path_tmp, hotspot_mutations, genome, group_by, cores):
        """
        Initialize HotspotFinder class
        Args:
            input_file (str): path to input data
            output_file (str): path to output file
            output_path_tmp (str): path to tmp directory
            hotspot_mutations (int): cutoff of mutations to define a hotspot
            genome (str): reference genome
            group_by (str): name of the column to group hotspots identification
            cores (int): number of cpu

        Returns:
            #TODO
        """

        self.input_file = input_file
        self.output_file = output_file
        self.output_path_tmp = output_path_tmp
        self.tmp_output_file_subs = os.path.join(output_path_tmp, 'subs.json')
        self.tmp_output_file_alts = os.path.join(output_path_tmp, 'alts.json')
        self.tmp_output_file_samples = os.path.join(output_path_tmp, 'samples.json')
        self.tmp_output_file_indels = os.path.join(output_path_tmp, 'indels.json')
        self.hotspot_mutations = hotspot_mutations
        self.genome = genome
        self.group_by = group_by
        self.cores = cores

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

    def preprocessing(self):
        """
        Load mutations file into hotspot dictionary with SNVs.

        Args:
            None

        Returns:
            None

        """

        substitutions_dict = defaultdict(lambda: defaultdict(int))
        mut_sample_alt_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
        samples_set = set()

        # Read mutations
        with open(self.tmp_output_file_samples, 'w') as samp_fd:
            with open(self.tmp_output_file_alts, 'w') as alts_fd:
                with open(self.tmp_output_file_indels, 'w') as indels_fd:
                    for row in readers.variants(
                            file=self.input_file,
                            required=['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE'],
                            extra=['GROUP', 'GROUP_BY', 'COHORT', 'CANCER_TYPE']
                    ):
                        chromosome = row['CHROMOSOME']
                        position = str(row['POSITION'])
                        ref = row['REF']
                        alt = row['ALT']
                        sample = row['SAMPLE']
                        samples_set.add(sample)
                        mutation = '{}_{}'.format(chromosome, position)
                        if not self.group_by:
                            cohort = 'cohort'
                        else:
                            cohort = row[self.group_by]

                        if ref != '-' and alt != '-':
                            if len(ref) == 1 and len(alt) == 1:  # 1 bp SNVs, DBS are skipped
                                if not mut_sample_alt_dict[cohort][mutation][sample]:  # First time sample is listed in mutation
                                    mutation_5 = '{}_{}'.format(chromosome, int(position) - 1)
                                    mutation_3 = '{}_{}'.format(chromosome, int(position) + 1)
                                    if not mut_sample_alt_dict[cohort][mutation_5][sample]:     # 5' SNV same patient
                                        if not mut_sample_alt_dict[cohort][mutation_3][sample]:    # 3' SNV same patient
                                            mut_sample_alt_dict[cohort][mutation][sample] = alt
                                            del mut_sample_alt_dict[cohort][mutation_5][sample]
                                            del mut_sample_alt_dict[cohort][mutation_3][sample]
                                            # FIXME this only works for files with one grouping condition
                                            alts_fd.write('{}\n'.format('\t'.join([chromosome, position, position, alt])))
                                            samp_fd.write('{}\n'.format('\t'.join([chromosome, position, position, sample])))
                                        else:
                                            del mut_sample_alt_dict[cohort][mutation_3][sample]
                                            del mut_sample_alt_dict[cohort][mutation][sample]
                                            logger.debug('sample {} has two SNVs in adjacent positions {} and {}'
                                                         'None is analyzed'.format(sample, mutation, mutation_3))
                                    else:
                                        del mut_sample_alt_dict[cohort][mutation_5][sample]
                                        del mut_sample_alt_dict[cohort][mutation][sample]
                                        logger.debug('sample {} has two SNVs in adjacent positions {} and {}.'
                                                     'None is analyzed'.format(sample, mutation_5, mutation))
                                else:
                                    if alt == mut_sample_alt_dict[cohort][mutation][sample]:
                                        logger.debug('sample {} has duplicated annotation of SNV at chr{}:{}>{}'
                                                       'This is analyzed a SNV'.format(
                                            sample, chromosome, position, alt))
                                    else:
                                        logger.debug('sample {0} has different nucleotide substitutions at '
                                                     'chr{1}:{2}>{3}|{4}'.format(
                                            sample, chromosome, position, alt, mut_sample_alt_dict[mutation][sample])
                                        )
                                        mut_sample_alt_dict[mutation][sample] = mut_sample_alt_dict[mutation][sample] + '|' + alt
                                        # FIXME this only works for files with one grouping condition
                                        alts_fd.write('{}\n'.format('\t'.join([chromosome, position, position, alt])))

                        else:     # This is an indel of any length
                            # if len(ref) == 1 and len(alt) == 1:
                            indel_type = 'insertion' if ref == '-' else 'deletion'
                            indels_fd.write('{}\n'.format(
                                '\t'.join([chromosome, position, position, '{}>{}'.format(ref, alt), indel_type])))
                            if mut_sample_alt_dict[cohort][mutation][sample]:
                                logger.debug('sample {} has an SNVs and an indel in the same position {}'
                                                .format(sample, position))
        # Save
        for k, v in mut_sample_alt_dict.items():
            if len(v) > 0:
                substitutions_dict[k] = len(v)
        self.json_convert_save(substitutions_dict, self.tmp_output_file_subs)

        for file in [self.tmp_output_file_samples, self.tmp_output_file_alts, self.tmp_output_file_indels]:
            self.tabix_convert_save(file)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input-file', default=None, required=True, type=click.Path(exists=True),
              help='File containing somatic mutations in TSV format')
@click.option('-o', '--output-path', default=None, required=True, help='Output directory')
@click.option('-mut', '--hotspot-mutations', type=click.IntRange(min=2, max=None, clamp=False), default=3,
              help='Cutoff of mutations to define a hotspot. Default is 3')
@click.option('-g', '--genome', default='hg19',
              type=click.Choice(['hg38', 'hg19']),
              help='Genome to use')
@click.option('-group', '--group_by', default=None, type=click.Choice(['GROUP', 'GROUP_BY', 'COHORT', 'CANCER_TYPE']),
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
    if not os.path.exists(output_path_tmp):
        os.makedirs(output_path_tmp, exist_ok=True)

    # Output file name
    output_file = os.path.join(output_path, 'hotspotfinder_results.txt.gz')

    # Log
    daiquiri.setup(level=LOGS[log_level], outputs=(
        daiquiri.output.STDERR,
        daiquiri.output.File(
            filename=output_file.split('.')[0] + '.log',
            directory=output_path
        )
    ))
    logger = daiquiri.getLogger()

    logger.info('HotspotFinder')
    logger.info('Initializing HotspotFinder...')
    logger.info('\n'.join([
        '\n'
        'Input data file: {}'.format(input_file),
        'Output results directory: {}'.format(output_path),
        'Genome: {}'.format(genome),
        'Cutoff hotspots mutations: {}'.format(hotspot_mutations)
    ]))



    # # Generate stage 1 hotspots (no annotations)
    # RawHotspots(
    #     input_path,
    #     output_path,
    #     hotspot_mutations,
    #     genome,
    #     cohort_level,
    #     cores
    # ).run()


if __name__ == '__main__':
    main()
