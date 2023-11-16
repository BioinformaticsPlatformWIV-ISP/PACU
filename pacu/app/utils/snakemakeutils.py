import re
from pathlib import Path
from typing import Dict, Any, List, Union

import yaml

from pacu.app.command import Command
from pacu.app.utils.loggingutils import logger


def generate_config_file(config_data: Dict[str, Any], output_dir: Path, output_basename: str = 'config.yml') -> str:
    """
    Generates a configuration file for Snakemake in YAML file format.
    :param config_data: Configuration data
    :param output_dir: Output directory
    :param output_basename: Output basename
    :return: Path to config file
    """
    config_path = output_dir / output_basename
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    with config_path.open('w') as handle:
        yaml.dump(config_data, handle)
    logger.info(f"Configuration file created: {config_path}")
    return str(config_path)


def run_snakemake(snakefile: Path, config_path: Path, targets: List[Path], dir_: Path, threads: int = 8) -> Command:
    """
    Helper function to run snakemake workflows.
    :param snakefile: Workflow snakefile
    :param config_path: Path to configuration file
    :param targets: Target output files
    :param dir_: Working directory
    :param threads: Number of threads to use
    :return: Snakemake command
    """
    # Create working directory
    if not dir_.exists():
        logger.debug(f'Creating working directory: {dir_}')
        dir_.mkdir(parents=True)

    # Create and run command
    command = Command(' '.join([
        'snakemake', *[str(x) for x in targets], '--snakefile', str(snakefile), '--configfile', str(config_path),
        '--cores', str(threads)]))
    command.run(dir_)

    # Check if command executed successfully
    if command.exit_code != 0:
        rule_failed = __get_failed_rule(command.stderr)
        logger.error(f"Failed at rule: {rule_failed if rule_failed is not None else 'n/a'}")
        raise RuntimeError(f'Error executing Snakemake')
    return command


def __get_failed_rule(stderr: str) -> Union[str, None]:
    """
    Returns the name of the rule that failed during Snakemake execution.
    :return: Name of the failed rule
    """
    for line in reversed(stderr.splitlines()):
        m = re.match(r'Error in rule (\w+):', line.strip())
        if not m:
            continue
        return m.group(1)
