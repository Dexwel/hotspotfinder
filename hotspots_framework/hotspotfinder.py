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




@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input-file', default=None, required=True, type=click.Path(exists=True),
              help='File containing somatic mutations in TSV format')
@click.option('-o', '--output-path', default=None, required=True, help='Output directory')
@click.option('-mut', '--hotspot-mutations', type=click.IntRange(min=2, max=None, clamp=False), default=3,
              help='Cutoff of mutations to define a hotspot. Default is 3')
@click.option('-g', '--genome', default='hg19',
              type=click.Choice(['hg38', 'hg19']),
              help='Genome to use')
@click.option('-group', '--group_by', default=None, type=click.Choice(['GROUP', 'CANCER_TYPE']),
              help='Header of the column to group hotspots identification')
@click.option('-c', '--cores', type=click.IntRange(min=1, max=os.cpu_count(), clamp=False), default=os.cpu_count(),
              help='Number of cores to use in the computation. By default it uses all available cores.')
@click.option('--log-level', default='info', help='Verbosity of the logger',
              type=click.Choice(['debug', 'info', 'warning', 'error', 'critical']))
def main(input_file, output_path, hotspot_mutations, genome, cores, group_by, log_level):
    """
    Look for mutational hotspots across the genome
    Args:
        input_file (str): path to input data
        output_path (str): path to output directory
        hotspot_mutations (int): cutoff of mutations to define a hotspot
        genome (str): reference genome
        cores (int): number of cpu
        group_by (str): name of the column to group hotspots identification
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
