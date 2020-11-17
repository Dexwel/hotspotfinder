import logging
import os
import sys
from collections import defaultdict

import click
import daiquiri

from hotspots_framework import __logger_name__
from hotspots_framework.hotspotfinder import HotspotFinder
from hotspots_framework import configuration

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

LOG_LEVELS = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}


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

    # File name
    file_name = input_mutations.split('/')[-1].split('.')[0]

    # Output directories
    if not os.path.exists(output_directory):
        os.makedirs(output_directory, exist_ok=True)

    # Log
    daiquiri.setup(level=LOG_LEVELS[log_level], outputs=(
        daiquiri.output.STDERR,
        daiquiri.output.File(
            filename=f'{file_name}.log',
            directory=output_directory
        )
    ))
    logger = daiquiri.getLogger(__logger_name__)
    logger.info('HotspotFinder')
    logger.info('Initializing HotspotFinder...')

    # Override configuration file parameters by the CLI arguments
    dd = lambda: defaultdict(dd)
    config_override = dd()
    if genome:
        print(genome)
        config_override['genome']['build'] = genome
    if mutations_cutoff:
        config_override['hotspot_mutations']['cutoff'] = mutations_cutoff
    if group_by:
        config_override['group']['groupby'] = group_by
    if cores:
        if cores > os.cpu_count():
            cores = os.cpu_count()
            logger.warning('Input cores greater than maximum CPU cores. Using maximum CPU cores instead')
        config_override['settings']['cores'] = cores

    # Read configuration file
    config = configuration.load(config_file=configuration_file, override=None)

    genome = config['genome']['build']
    mutations_cutoff = config['hotspot_mutations']['cutoff']
    group_by = config['group']['groupby']
    cores = config['settings']['cores']


    # Mappability
    if os.path.isfile(config['mappability']['mappable_regions']):
        mappable_regions = config['mappability']['mappable_regions']
    else:
        if config['mappability']['mappable_regions'] == 'bgdata':
            mappable_regions = 'bgdata'
        else:
            logger.error(f"Mappable regions file does not exist: {config['mappability']['mappable_regions']}")
            sys.exit(-1)

    if os.path.isfile(config['mappability']['blacklisted_regions']):
        blacklisted_regions = config['mappability']['blacklisted_regions']
    else:
        if config['mappability']['blacklisted_regions'] == 'bgdata':
            blacklisted_regions = 'bgdata'
        else:
            logger.error(f"Blacklisted regions file does not exist: {config['mappability']['blacklisted_regions']}")
            sys.exit(-1)

    # Population variants
    if os.path.isdir(config['polymorphisms']['population_variants']):
        population_variants = config['polymorphisms']['population_variants']
    else:
        if config['polymorphisms']['population_variants'] == 'bgdata':
            population_variants = 'bgdata'
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