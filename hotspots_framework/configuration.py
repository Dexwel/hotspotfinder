import logging
from os import path

from bgconfig import BGConfig

from hotspots_framework import __logger_name__


LOGGER = logging.getLogger(__logger_name__)


def load(config_file, override=None):
    """
    Load the configuration file and checks the format.

    Args:
        config_file: configuration file path
        override: override values

    Returns:
        :class:`bgconfig.BGConfig`: configuration as a :obj:`dict`

    """
    if not config_file:
        config_file = path.join(path.dirname(__file__), "hotspot.conf")
    config_template = path.join(path.dirname(__file__), "hotspot.conf.template")
    config_spec = path.join(path.dirname(__file__), "hotspot.conf.spec")

    config = BGConfig(
        config_template, config_file=config_file, config_spec=config_spec, use_env_vars=True,
        override_values=override, unrepr=False, use_bgdata=True)

    return config
