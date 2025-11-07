import binascii
import re
import tempfile
from importlib.resources import files
from pathlib import Path
from typing import Union

import numpy
import pandas as pd
import vcf
from Bio import Phylo

from pacu.app.command import Command
from pacu.app.utils.loggingutils import logger


def is_new_region(record: pd.Series) -> bool:
    """
    Function to check if a region overlaps the previous one.
    :param record: BED file record
    :return: True if it is a new region (not overlapping the previous one), False otherwise
    """
    if numpy.isnan(record['shift_pos']):
        return False
    elif int(record['pos']) - 1 != int(record['shift_pos']):
        return True
    elif record['chr'] != record['chr']:
        return True
    return False


def calculate_overlaps(size: int, bed_phages: Path, bed_gubbins: Path, bed_depth: Path) -> pd.DataFrame:
    """
    Calculates the size of the overlaps between the BED files.
    :param bed_phages: BED file with phage regions
    :param bed_gubbins: BED file with recombination detected by Gubbins
    :param bed_depth: BED file with low-depth positions
    :param size: Reference genome size
    """
    data_overlap = pd.DataFrame(data={'pos': range(size)})

    bed_dict = {'phages': bed_phages, 'gubbins': bed_gubbins, 'depth': bed_depth}
    for key, path_bed in bed_dict.items():
        data_overlap[key] = False
        data_bed = pd.read_table(path_bed, names=['seq_id', 'start', 'end', 'name'])
        logger.info(f'{len(data_bed)} regions parsed from {path_bed.name}')
        for region in data_bed.itertuples():
            # noinspection PyUnresolvedReferences
            pos = (data_overlap['pos'][(data_overlap['pos'] > region.start) & (data_overlap['pos'] <= region.end)])
            data_overlap.loc[pos, key] = True

    data_overlap = data_overlap[data_overlap[bed_dict.keys()].any(axis=1)]
    logger.info(f'{len(data_overlap):,} positions removed in total')
    return data_overlap


def count_covered_positions(bed_file: Path) -> int:
    """
    Counts the number of genomic positions from a BED file.
    - *Important!*: Only works for BED files with non-overlapping features
    :param bed_file: Input BED file
    :return: Number of covered positions
    """
    data_in = pd.read_table(bed_file, usecols=[0, 1, 2], names=['chr', 'start', 'end'])
    data_in['interval_size'] = data_in['end'] - data_in['start']
    return sum(data_in['interval_size'])


def count_overlap(bed_file_a: Path, bed_file_b: Path) -> int:
    """
    Counts the overlap between two bed files.
    :param bed_file_a: Input bed file A
    :param bed_file_b: Input bed file B
    :return: Number of overlaps
    """
    with tempfile.TemporaryDirectory(prefix='pacu_') as dir_temp:
        command = Command(' '.join([
            f'bedtools multiinter -i {bed_file_a} {bed_file_b}',
            '|', "awk -F'\t' 'BEGIN{SUM=0}{SUM+=$3-$2 }END{print SUM}'"
        ]))
        command.run(Path(dir_temp))
    return int(command.stdout)


def calculate_distance(row: pd.Series) -> Union[int, None]:
    """
    Calculates the distance to the closest SNP.
    :param row: Table row
    :return: Distance to SNP
    """
    if pd.isna(row['chr_next']) or pd.isna(row['chr_prev']):
        return None

    # Distance to previous SNP
    dist_to_prev = None
    if row['chr_prev'] == row['chr']:
        dist_to_prev = row['pos'] - row['pos_prev']

    # Distance to next SNP
    dist_to_next = None
    if row['chr_next'] == row['chr']:
        dist_to_next = row['pos_next'] - row['pos']

    # Single variant on chromosome
    if dist_to_next is None and dist_to_prev is None:
        return None

    return min(x for x in (dist_to_prev, dist_to_next) if x is not None)


def filter_snp_distance(vcf_in: Path, vcf_out: Path, min_dist: int = 10) -> None:
    """
    Soft-filters SNPs that are located within less than other SNPs.
    Note that the filter considers all SNPs (including those at filtered positions)
    :param min_dist: Minimum distance between SNPs
    :param vcf_in: Input VCF file
    :param vcf_out: Output VCF file
    :return: None
    """
    # Parse input SNPs
    records_snps = []
    with vcf_in.open() as handle:
        for row in vcf.Reader(handle):
            records_snps.append({'chr': row.CHROM, 'pos': row.POS})
    data_snps = pd.DataFrame(records_snps)
    logger.info(f'{len(data_snps):,} SNPs parsed')

    # Calculate distance to closest SNP
    if len(data_snps) > 0:
        data_snps['pos_next'] = data_snps['pos'].shift(-1)
        data_snps['chr_next'] = data_snps['chr'].shift(-1)
        data_snps['pos_prev'] = data_snps['pos'].shift(1)
        data_snps['chr_prev'] = data_snps['chr'].shift(1)
        data_snps['dist'] = data_snps.apply(lambda x: calculate_distance(x), axis=1)

        # Extract the keys of the variants that need to be filtered
        keys_removed = [(x['chr'], x['pos']) for x in data_snps[data_snps['dist'].apply(
            lambda x: x is not None and x <= min_dist)].to_dict('records')]
    else:
        keys_removed = []
    logger.info(f'Filtering {len(keys_removed):,} variants based on distance')

    # Save output file
    with vcf_in.open() as handle_in, vcf_out.open('w') as handle_out:
        vcf_reader = vcf.Reader(handle_in)
        vcf_writer = vcf.Writer(handle_out, template=vcf_reader)

        # Update filtering records
        for record in vcf_reader:
            if (record.CHROM, record.POS) in keys_removed:
                record.add_filter('distance')

            # Write to output file
            vcf_writer.write_record(record)
        vcf_writer.close()
    logger.info(f'Output VCF file created: {vcf_out}')


def sort_snp_dist_matrix(path_nwk: Path, path_tsv: Path, path_out: Path) -> None:
    """
    Sorts the input SNP distance matrix based on the order of the nodes in the input phylogeny.
    :param path_nwk: Path to the input Newick file
    :param path_tsv: Path to the input TSV file
    :param path_out: Output path
    :return: None
    """
    with path_nwk.open() as handle:
        tree = Phylo.read(handle, 'newick')
        tree.root_at_midpoint()
        node_order = [n.name for n in tree.find_clades() if n.is_terminal()]

    # Parse distance information
    data_in_dist = pd.read_table(path_tsv, index_col=0)

    # Sort rows & columns
    data_in_dist = data_in_dist.loc[node_order]
    data_in_dist = data_in_dist[node_order]
    data_in_dist.to_csv(path_out, sep='\t', index=True, index_label='isolate')


def plot_newick_phylogeny(path_nwk: Path, path_out: Path, width: int = 600, height: int = 600) -> None:
    """
    Plots the input Newick phylogeny in PNG format.
    :param path_nwk: Input Newick file
    :param path_out: Output PNG path
    :param width: Image width
    :param height: Image height
    :return: None
    """
    with tempfile.NamedTemporaryFile(prefix='pacu_') as file_:
        # Convert the tree to Nexus format
        Phylo.convert(str(path_nwk), 'newick', file_.name, 'nexus')

        # Add the figtree code
        path_template = Path(str(files('pacu').joinpath('resources/figtree_template.txt')))
        with open(file_.name, 'a') as handle_out, open(path_template) as handle_in:
            handle_out.write('\n')
            for line in handle_in.readlines():
                handle_out.write(line)

        # Create the visualization
        command = Command(' '.join([
            'figtree',
            '-graphic PNG',
            f'-width {width}',
            f'-height {height}',
            str(file_.name),
            str(path_out.absolute())
        ]))
        command.run(path_out.parent)
        if not command.exit_code == 0:
            raise RuntimeError(f"Error visualizing tree: {command.stderr}")
        logger.info(f'Visualization exported to: {path_out}')


def sanitize_input_name(name: str, extension: str) -> str:
    """
    Sanitizes the input file name.
    :param name: Name
    :param extension: Expected file extension (e.g., 'bam' or 'fasta')
    :return: None
    """
    invalid_chars = '/!@#$\\'

    # Replace spaces by dashes
    name = name.replace(' ', '_')

    # Avoid double dot before the extension
    if name.endswith('.'):
        name = name[:-1]

    # Add extension
    if not name.endswith(f'.{extension}'):
        return f'{name}.{extension}'

    # Return sample name without invalid characters
    return ''.join(c for c in name if c not in invalid_chars)


PATTERNS_FQ_PE = [
    r'(.*)_1P?\.(fastq|fq)(\.gz)?',
    r'(.*)_L\d+_R1_\d+\.(fastq|fq)(\.gz)?',
    r'(.*)_R1.(fastq|fq)(\.gz)?',
]

def determine_name_from_fq(fq_ont: Path = None, fq_illumina_1p: Path = None) -> str:
    """
    Determines the sample name from the FASTQ input.
    :param fq_ont: Input ONT FASTQ file
    :param fq_illumina_1p: Input illumina forward FASTQ file
    :return: Sample name
    """
    if fq_ont is not None:
        return re.sub(r'\.(fastq|fq)(\.gz)?', '', fq_ont.name)
    elif fq_illumina_1p is not None:
        for pattern in PATTERNS_FQ_PE:
            m = re.match(pattern, fq_illumina_1p.name)
            if not m:
                continue
            return m.group(1)
        logger.warning(
            'The input FASTQ filename does not match any of the supported formats for sample name determination')
        return re.sub(r'.(fastq|fq)(\.gz)?', '', fq_illumina_1p.name)
    else:
        raise ValueError('No FASTQ file provided')


def is_gzipped(path: Path) -> bool:
    """
    Checks if the given file is compressed with gzip.
    :param path: Path
    :return: True if gzipped, False otherwise
    """
    with path.open('rb') as handle:
        magic_number = binascii.hexlify(handle.read(2))
    return magic_number == b'1f8b'
