import datetime
from pathlib import Path
from typing import Any, Dict

import pandas as pd
from matplotlib import pyplot
from upsetplot import plot

from sciensnp.app.report.htmlcitation import HtmlCitation
from sciensnp.app.report.htmlelement import HtmlElement
from sciensnp.app.report.htmlreportsection import HtmlReportSection


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
        ['Analysis date:', datetime.datetime.now().strftime('%d/%m/%Y - %X')],
        ['Nb. of samples:', str(len(config_['input'].keys()))],
        ['Reference genome:', Path(config_['reference']['fasta']).name],
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
        ['Min. allele frequency:', f"{config_['filters']['min_af']:.2f}"],
        ['Min. depth:', str(config_['filters']['min_depth'])],
        ['Min. SNP quality:', str(config_['filters']['min_qual'])],
    ], table_attributes=[('class', 'information')])
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


def create_region_filtering_section(path_stats: Path, path_png: Path) -> HtmlReportSection:
    """
    Creates the section with the region filtering results.
    :param path_stats: Path to TSV stats file
    :param path_png: Visualization of region overlap
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

    # Add path to image
    relative_path_png = Path(path_png.name)
    section.add_file(path_png, relative_path_png)
    section.add_line_break()
    return section


def create_tree_section(path_nwk: Path, path_png: Path, path_snp_matrix: Path) -> HtmlReportSection:
    """
    Creates the section with the visualization of the tree.
    :param path_nwk: Newick tree
    :param path_png: Tree visualization
    :param path_snp_matrix: SNP matrix
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
    keys = ['Croucher_2015-gubbins', 'Danecek_2021-bcftools', 'Li_2018-minimap2', 'Quinlan_2010-bedtools']
    for citation_key in keys:
        section.add_html_object(HtmlCitation.parse_from_json(citation_key))
    return section
