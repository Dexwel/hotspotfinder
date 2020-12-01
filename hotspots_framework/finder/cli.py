import logging
import os
from collections import defaultdict

import click
import daiquiri

from hotspots_framework import __logger_name__
from hotspots_framework.finder import HotspotFinderError
from hotspots_framework.finder.hotspotfinder import HotspotFinder
from hotspots_framework import configuration

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

LOG_LEVELS = {
    'debug': logging.DEBUG,
    'info': logging.INFO,
    'warning': logging.WARNING,
    'error': logging.ERROR,
    'critical': logging.CRITICAL
}

# TODO remove genome and cores
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-i', '--input-mutations', default=None, required=True, type=click.Path(exists=True),
              help='User input file containing somatic mutations in TSV format')
@click.option('-o', '--output-directory', default=None, required=True, help='Output directory')
@click.option('-g', '--genome', default=None, type=click.Choice(['hg38']), help='Genome to use')
@click.option('-mcutoff', '--mutations-cutoff', type=click.IntRange(min=3, max=None, clamp=False), default=None,
              help='Cutoff of mutations to define a hotspot. Default is 3')
@click.option('-scutoff', '--samples-cutoff', type=click.IntRange(min=3, max=None, clamp=False), default=None,
              help='Cutoff of number of mutated samples to define a hotspot. Default is 3')
@click.option('-group', '--group-by', default=None,
              type=click.Choice(['GROUP', 'GROUP_BY', 'COHORT', 'CANCER_TYPE', 'PLATFORM']),
              help='Header of the column to group hotspots identification')
@click.option('-c', '--cores', type=click.IntRange(min=1, max=os.cpu_count(), clamp=False), default=None,
              help='Number of cores to use in the computation. By default it uses all available cores.')
@click.option('-conf', '--configuration-file', default=None,
              required=False, type=click.Path(exists=True), help='User input configuration file')
@click.option('--log-level', default='info', help='Verbosity of the logger',
              type=click.Choice(['debug', 'info', 'warning', 'error', 'critical']))
def main(
        input_mutations,
        output_directory,
        mutations_cutoff,
        samples_cutoff,
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
        samples_cutoff (int): cutoff of number of mutated samples to define a hotspot
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

    # Suppress log messages from some libraries
    daiquiri.getLogger('bgdata').setLevel(logging.WARNING)

    logger = daiquiri.getLogger(__logger_name__)
    logger.info('HotspotFinder')
    logger.info('Initializing HotspotFinder...')

    # Override configuration file parameters by the CLI arguments
    dd = lambda: defaultdict(dd)
    config_override = dd()
    if genome:
        config_override['genome'] = genome
    if mutations_cutoff:
        config_override['finder']['mutations_cutoff'] = mutations_cutoff
    if samples_cutoff:
        config_override['finder']['samples_cutoff'] = samples_cutoff
    if group_by:
        config_override['finder']['groupby'] = group_by
    if cores:
        if cores > os.cpu_count():
            cores = os.cpu_count()
            logger.warning('Input cores greater than maximum CPU cores. Using maximum CPU cores instead')
        config_override['cores'] = cores

    # Read configuration file
    config = configuration.load(config_file=configuration_file, override=config_override)

    # TODO if not config['finder']['annotate'], skip next lines?
    #config['finder']['annotate']

    # Mappability
    if not os.path.isfile(config['mappable_regions']):
        raise HotspotFinderError(f"Mappable regions file does not exist: {config['mappable_regions']}")

    if not os.path.isfile(config['blacklisted_regions']):
        raise HotspotFinderError(f"Blacklisted regions file does not exist: {config['blacklisted_regions']}")

    # Population variants
    if not os.path.isfile(config['population_variants']):
        raise HotspotFinderError(f"Population variants file does not exist: {config['population_variants']}")

    # Repeats
    if not os.path.isfile(config['repeats']):
        raise HotspotFinderError(
            f"Repeats file does not exist: {config['repeats']}")

    # Ig and TR regions
    if not os.path.isfile(config['ig_tr_regions']):
        raise HotspotFinderError(
            f"IG-TR regions file does not exist: {config['ig_tr_regions']}")

    # Genomic elements
    if not os.path.isfile(config['genomic_elements']):
        raise HotspotFinderError(f"Genomic regions file does not exist: {config['genomic_elements']}")

    # Output file names
    compression = '.gz' if config['output_format'].endswith('.gz') else ''
    output_file_results = os.path.join(output_directory, f'{file_name}.results.{config["output_format"]}')
    output_file_warning = os.path.join(output_directory, f'{file_name}.warningpositions.txt{compression}')

    logger.info('Analysis parameters loaded')

    logger.info('\n'.join([
        '\n'
        f"* Input mutations file: {input_mutations}",
        f"* Annotate hotspots: {config['finder']['annotate']}",
        f"* Mapable regions file: {config['mappable_regions']}",
        f"* Blacklisted regions file: {config['blacklisted_regions']}",
        f"* Population variants file: {config['population_variants']}",
        f"* Repeats file: {config['repeats']}",
        f"* Ig-TR regions file: {config['ig_tr_regions']}",
        f"* Genomic elements file: {config['genomic_elements']}",
        f"* Output results directory: {output_directory}",
        f"* Output format: {config['output_format']}",
        f"* Cutoff hotspot mutations: {config['finder']['mutations_cutoff']}",
        f"* Cutoff hotspot mutated samples: {config['finder']['samples_cutoff']}",
        f"* Hotspots split alternates: {config['finder']['split_alternates']}",
        f"* Remove unknown nucleotides: {config['finder']['remove_unknown_reference_nucleotides']}",
        f"* Remove nonannotated hotspots: {config['finder']['remove_nonannotated_hotspots']}",
        f"* Genome: {config['genome']}",
        f"* Group analysis by: {config['finder']['groupby']}",
        f"* Cores: {config['cores']}",
    ]))

    logger.info('Loading annotations files...')
    # Run HotspotFinder
    experiment = HotspotFinder(
        input_mutations,
        output_file_results,
        output_file_warning,
        config
    )

    experiment.run()
    logger.info('HotspotFinder analysis finished!')


if __name__ == '__main__':
    main()
