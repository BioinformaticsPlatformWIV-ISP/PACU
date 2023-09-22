from importlib.resources import files
from pathlib import Path
from typing import Optional, Dict

from sciensnp.app.command import Command
from sciensnp.app.utils.loggingutils import logger

OPTIONS = {
    'branch_swap_filter': ['None', 'Very Weak', 'Weak', 'Moderate', 'Strong', 'Very Strong'],
    'missing_data_treatment': ['Complete deletion', 'Use all sites', 'Partial deletion'],
}

SUBSTITUTION_MODELS = {
    'JC': 'Jukes-Cantor model',
    'K2': 'Kimura 2-parameter model',
    'T92': 'Tamura 3-parameter model',
    'HKY': 'Hasegawa-Kishino-Yano model',
    'TN93': 'Tamura-Nei model',
    'GTR': 'General Time Reversible model'
}

HEURISTIC_METHODS = {
    'NNI': 'Nearest-Neighbor-Interchange (NNI)',
    'SPR3': 'Subtree-Pruning-Regrafting - Fast (SPR level 3)',
    'SPR5': 'Subtree-Pruning-Regrafting - Extensive (SPR level 5)'
}

INITIAL_TREE = {
    'NJ_BioNJ': 'Make initial tree automatically (Default - NJ/BioNJ)',
    'Max_Parsimony': 'Make initial tree automatically (Maximum Parsimony)',
    'NJ': 'Make initial tree automatically (Neighbor Joining)',
    'BioNJ': 'Make initial tree automatically (BioNJ)'
}

RATES_AMONG_SITES = {
    'U': 'Uniform Rates',
    'G': 'Gamma Distributed (G)',
    'I': 'Has Invariant Sites (I)',
    'G+I': 'Gamma Distributed With Invariant Sites (G+I)'
}


def __create_config_file_model_selection(
        dir_: Path, branch_swap_filter: str, missing_data_treatment: str, site_cov_cutoff: Optional[int],
        threads: int = 4) -> Path:
    """
    Generates the config file for the model selection.
    :param dir_: Working directory
    :param missing_data_treatment: Missing data treatment option
    :param site_cov_cutoff: Site coverage cut-off threshold
    :param threads: Number of threads to use
    :return: Path to config file
    """
    # Parse the template
    path_template = files('sciensnp').joinpath('resources/mega/model_sel_ml_nucleotide_template.mao')
    with path_template.open() as handle:
        template = handle.read()
    logger.info(f'Template parsed')

    # Create filled-in template
    path_out = dir_ / 'config-model_selection.mao'
    with path_out.open('w') as handle:
        handle.write(template.format(
            branch_swap_filter=branch_swap_filter,
            missing_data_treatment=missing_data_treatment,
            site_coverage_cutoff=site_cov_cutoff if site_cov_cutoff is not None else 'Not Applicable',
            threads=threads
        ))
    return path_out


def run_model_selection(path_fasta: Path, path_out: Path, dir_: Path, branch_swap_filter: str,
                        missing_data_treatment: str, site_cov_cutoff: Optional[int], threads: int = 4) -> None:
    """
    Runs model selection on the given SNP matrix.
    :param path_fasta: Input FASTA file
    :param path_out: Output path
    :param dir_: Working directory
    :param branch_swap_filter: Branch swap filter option
    :param missing_data_treatment: Missing data treatment option
    :param site_cov_cutoff: Site coverage cut-off
    :param threads: Number of threads to use
    :return: None
    """
    # Create configuration file
    path_config = __create_config_file_model_selection(
        dir_, branch_swap_filter, missing_data_treatment, site_cov_cutoff, threads)

    # Run MEGA model selection
    command = Command(' '.join([
        'megacc', f'-d {path_fasta}', f'-a {path_config}', f'-o {path_out}'
    ]))
    command.run(path_out.parent)
    if not command.exit_code == 0:
        raise RuntimeError(f'Error running model selection: {command.stderr}')


def parse_model_selection_csv(path_csv: Path) -> Dict:
    """
    Parses the model selection output CSV file.
    :return: Information for the best selected model
    """
    model_info = {}
    with path_csv.open() as handle:
        first_line = handle.readlines()[1]
        # The line has the following structure: K2+G+I, ...
        # The model is the first part (K2), +G means Gamma categories per site, +I means invariant sites
        complete_model = first_line.split(',')[0]
        model_info['model'] = first_line.split(',')[0].split('+')[0]
        model_info['model_full'] = SUBSTITUTION_MODELS[model_info['model']]

        # Rates
        complete_rates = '+'.join(complete_model.split('+')[1:])
        model_info['rates_among_sites'] = 'U' if complete_rates == '' else complete_rates
        model_info['rates_among_sites_full'] = RATES_AMONG_SITES[model_info['rates_among_sites']]
    return model_info


def __create_config_file_tree_construction(
        dir_: Path, model: str, branch_swap_filter: str, missing_data_treatment: str, site_cov_cutoff: Optional[int],
        bootstrap_replicates: int, heuristic_method: str, initial_tree: str, rates_among_sites: str,
        gamma_categories: int = 5, threads: int = 4) -> Path:
    """
    Generates the config file for the model selection.
    :param dir_: Working directory
    :param model: Nucleotide substitution model
    :param missing_data_treatment: Missing data treatment option
    :param site_cov_cutoff: Site coverage cut-off threshold
    :param bootstrap_replicates: Number of bootstrap replications
    :param heuristic_method: Heuristic method
    :param initial_tree: Initial tree to start from
    :param rates_among_sites: Rates among sites
    :param gamma_categories: Number of Gamma categories
    :param threads: Number of threads to use
    :return: Path to config file
    """
    # Parse the template
    path_template = files('sciensnp').joinpath('resources/mega/infer_ml_nucleotide_template.mao')
    with path_template.open() as handle:
        template = handle.read()
    logger.info(f'Template parsed')

    # Create filled-in template
    path_out = dir_ / 'config-model_selection.mao'
    with path_out.open('w') as handle:
        handle.write(template.format(
            bootstrap_replications=bootstrap_replicates,
            model=model,
            heuristic_method=HEURISTIC_METHODS[heuristic_method],
            initial_tree=INITIAL_TREE[initial_tree],
            rates_among_sites=RATES_AMONG_SITES[rates_among_sites],
            branch_swap_filter=branch_swap_filter,
            gamma_categories=gamma_categories,
            missing_data_treatment=missing_data_treatment,
            site_coverage_cutoff=site_cov_cutoff,
            threads=threads
        ))
    return path_out


def run_tree_building(path_fasta: Path, path_csv: Path, path_out: Path, dir_: Path, branch_swap_filter: str,
                      missing_data_treatment: str, site_cov_cutoff: Optional[int], bootstrap_replicates: int,
                      heuristic_method: str, initial_tree: str, gamma_categories: int = 5,
                      threads: int = 4) -> None:
    """
    Runs the ML tree reconstruction using MEGA.
    :return: None
    """
    # Parse model
    data_model_selection = parse_model_selection_csv(path_csv)

    # Create configuration file
    path_config = __create_config_file_tree_construction(
        dir_, data_model_selection['model_full'], branch_swap_filter, missing_data_treatment, site_cov_cutoff,
        bootstrap_replicates, heuristic_method, initial_tree, data_model_selection['rates_among_sites'],
        gamma_categories, threads)

    # Run MEGA model selection
    command = Command(' '.join([
        'megacc', f'-d {path_fasta}', f'-a {path_config}', f'-o {path_out}'
    ]))
    command.run(path_out.parent)
    if not command.exit_code == 0:
        raise RuntimeError(f'Error running tree construction: {command.stderr}')
