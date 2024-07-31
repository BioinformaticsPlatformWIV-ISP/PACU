import datetime
from pathlib import Path
from typing import Any, Dict

import pandas as pd
from matplotlib import pyplot
from upsetplot import plot

from pacu import version
from pacu.app.report.htmlcitation import HtmlCitation
from pacu.app.report.htmlelement import HtmlElement
from pacu.app.report.htmlreportsection import HtmlReportSection
from pacu.app.report.htmltablecell import HtmlTableCell


def create_upsetplot_overlap(data_overlap: pd.DataFrame, png_out: Path) -> None:
    """
    Creates an UpSet plot for the region filtering overlaps.
    :param data_overlap: Overlap data
    :param png_out: Output PNG file
    :return: None
    """
    plot(data_overlap.groupby(by=['phages', 'gubbins', 'depth']).size(), show_counts=True)
    pyplot.savefig(png_out, dpi=300)


def create_analysis_info_section(config_: Dict[str, Any]) -> HtmlReportSection:
    """
    Creates the analysis info section.
    :param config_: Configuration data
    :return: Analysis info section
    """
    section = HtmlReportSection('Analysis info')

    # Analysis info
    section.add_table([
        ['Workflow version:', version.__version__],
        ['Analysis date:', datetime.datetime.now().strftime('%d/%m/%Y - %X')],
        ['Nb. of samples:', str(len(config_['input'].keys()))],
        ['Phylogenetic method:', config_.get('phylogeny_method', 'iqtree')],
    ], table_attributes=[('class', 'information')])

    # Citation disclaimer
    section.add_header('Disclaimer', 2)
    section.add_paragraph('If you use this pipeline for your scientific work, please cite:')
    section.add_html_object(HtmlCitation.parse_from_json('Bogaerts_2023-ont_outbreak'))

    return section


def create_parameter_section(config_: Dict[str, Any]) -> HtmlReportSection:
    """
    Creates the parameter section.
    :param config_: Configuration data
    :return: Section
    """
    section = HtmlReportSection('Parameters')
    section.add_table([
        ['Min. global depth:', f"{config_['depth']['min_depth']:.2f}"],
        ['Min. SNP allele frequency:', f"{config_['filters']['min_af']:.2f}"],
        ['Min. SNP depth:', str(config_['filters']['min_depth'])],
        ['Min. SNP quality:', str(config_['filters']['min_qual'])],
        ['Min. SNP distance:', str(config_['filters']['min_dist'])],
    ], table_attributes=[('class', 'information')])
    return section


def __get_colored_cell_depth(value: int) -> HtmlTableCell:
    """
    Returns a colored cell for the mapping output.
    :param value: Depth value
    :return: Table cell
    """
    value_str = f'{value:,}'
    if value >= 20:
        color = 'green'
    elif value < 10:
        color = 'red'
    else:
        color = 'orange'
    return HtmlTableCell(value_str, color=color)


def __get_colored_cell_coverage(value: int) -> HtmlTableCell:
    """
    Returns a colored cell for the mapping output.
    :param value: Coverage value
    :return: Table cell
    """
    value_str = f'{value:.2f}%'
    if value >= 95:
        color = 'green'
    elif value < 90:
        color = 'red'
    else:
        color = 'orange'
    return HtmlTableCell(value_str, color=color)


def create_mapping_section(path_stats: Path, path_ref: Path) -> HtmlReportSection:
    """
    Creates the read mapping section.
    :param path_stats: Path to the TSV stats file
    :param path_ref: Path to reference genome
    :return: Section
    """
    section = HtmlReportSection('Read mapping')
    data_vc = pd.read_table(path_stats)
    header = ['Isolate', 'Median  depth', '% covered']
    section.add_paragraph(f'Reference genome: {path_ref.name}')
    section.add_table([[
        row['key'],
        __get_colored_cell_depth(row['median_depth']),
        __get_colored_cell_coverage(row['perc_covered']),
    ] for row in data_vc.to_dict('records')], header, [('class', 'data')])
    section.add_paragraph(
        """
        Datasets with a median depth between 20x-10x are displayed in orange (warning). Datasets with a median depth 
        <10x are displayed in red and should be re-sequenced. The genome coverage cut-offs are 95% (warning) and 90% 
        (error).
        """
    )
    return section


def create_variant_calling_section(path_stats: Path) -> HtmlReportSection:
    """
    Creates the variant calling section.
    :param path_stats: Path to the TSV stats file
    :return: Section
    """
    section = HtmlReportSection('Variant calling and filtering')
    data_vc = pd.read_table(path_stats)
    columns = [
        {'key': 'key', 'header': 'Isolate'},
        {'key': 'nb_variants', 'header': 'Number of variants', 'fmt': lambda x: f'{x:,}'},
        {'key': 'nb_variants_filtered', 'header': 'Number of variants (Filtered)', 'fmt': lambda x: f'{x:,}'},
        {'key': 'perc_passed', 'header': '% passing filtering', 'fmt': lambda x: f'{x:.2f}'},
    ]
    # noinspection PyCallingNonCallable
    section.add_table([
        [col['fmt'](row[col['key']]) if col.get('fmt') is not None else row[col['key']] for col in columns] for
        row in data_vc.to_dict('records')], [c['header'] for c in columns], [('class', 'data')])
    return section


def create_region_filtering_section(path_stats: Path, path_png: Path, skip_gubbins: bool = False) -> HtmlReportSection:
    """
    Creates the section with the region filtering results.
    :param path_stats: Path to TSV stats file
    :param path_png: Visualization of region overlap
    :param skip_gubbins: Boolean to indicate whether Gubbins was skipped.
    :return: Section
    """
    section = HtmlReportSection('Region filtering')

    # Add stats
    data_regions = pd.read_table(path_stats)
    columns = [
        {'key': 'key', 'header': 'BED file'},
        {'key': 'nb_covered_positions', 'header': 'Positions covered', 'fmt': lambda x: f'{x:,}'},
        {'key': 'perc_ref_covered', 'header': '% covered', 'fmt': lambda x: f'{x:.2f}'},
    ]
    # noinspection PyCallingNonCallable
    section.add_table([
        [col['fmt'](row[col['key']]) if col.get('fmt') is not None else row[col['key']] for col in columns] for
        row in data_regions.to_dict('records')], [c['header'] for c in columns], [('class', 'data')])

    if skip_gubbins:
        section.add_warning_message('Gubbins was disabled')

    # Add the path to the image
    relative_path_png = Path(path_png.name)
    section.add_file(path_png, relative_path_png)
    section.add_line_break()
    return section


def create_tree_section(path_nwk: Path, path_png: Path, path_snp_matrix: Path, model_name: str) -> HtmlReportSection:
    """
    Creates the section with the visualization of the tree.
    :param path_nwk: Newick tree
    :param path_png: Tree visualization
    :param path_snp_matrix: SNP matrix
    :param model_name: Name of the model used
    :return: Section
    """
    section = HtmlReportSection('Phylogeny')

    # Add SNP matrix
    relative_path_snp_matrix = Path(path_snp_matrix.name)
    section.add_header('SNP matrix', 4)
    section.add_link_to_file('Download (FASTA)', relative_path_snp_matrix)
    section.add_file(path_snp_matrix, relative_path_snp_matrix)

    # Add visualization
    section.add_header('Tree', 4)
    section.add_paragraph(f'Selected model: {model_name}')
    relative_path_png = Path(path_png.name)
    section.add_html_object(HtmlElement('img', attributes=[('src', str(relative_path_png)), ('border', '1')]))
    section.add_file(path_png, relative_path_png)
    section.add_line_break()

    # Add download link (Newick file)
    relative_path = Path(path_nwk.name)
    section.add_file(path_nwk, relative_path)
    section.add_link_to_file('Download (NWK)', relative_path)
    return section


def create_snp_distances_section(path_distances: Path) -> HtmlReportSection:
    """
    Creates the section with pair-wise SNP distances.
    :param path_distances: Path with pair-wise SNP distances
    :return: Section
    """
    section = HtmlReportSection(f'Pair-wise SNP distances')
    data_snp_dists = pd.read_table(path_distances)
    header = [HtmlElement('th', '')] + [
        HtmlElement('th', col, [('style', 'writing-mode: vertical-lr;')]) for col in data_snp_dists.columns[1:]]
    section.add_table(
        [header] +
        [[record[col] for col in data_snp_dists.columns] for record in data_snp_dists.to_dict('records')],
        None, [('class', 'data')])
    return section


def create_citations_section() -> HtmlReportSection:
    """
    Creates the report section with the citations.
    :return: Section
    """
    section = HtmlReportSection('Citations')
    keys = [
        'Croucher_2015-gubbins',
        'Danecek_2021-bcftools',
        'Li_2018-minimap2',
        'Quinlan_2010-bedtools',
        'Nguyen_2015-iq_tree',
        'Kalyaanamoorthy_2017-iq_tree_model',
        'Kumar_2018-megax'
    ]
    for citation_key in keys:
        section.add_html_object(HtmlCitation.parse_from_json(citation_key))
    return section
