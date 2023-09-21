from importlib.resources import files
from pathlib import Path
from typing import Any, Dict

import yaml

from sciensnp.app.utils.loggingutils import logger


def load_config_data() -> Dict[str, Any]:
    """
    Loads the configuration data.
    :return: None
    """
    path_config = files('sciensnp').joinpath('config/config.yml')
    logger.debug(f'Loading configuration data from: {path_config}')
    with Path(str(path_config)).open() as handle:
        return yaml.safe_load(handle)


config = load_config_data()
