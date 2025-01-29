from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Optional

import pandas as pd
import vcf
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# noinspection PyProtectedMember
from vcf.model import _Record as VCFRecord

from pacu.app.utils.loggingutils import logger


@dataclass(unsafe_hash=True, frozen=True, order=True)
class SNPPosition:
    contig: str
    position: int
    reference_base: chr


def parse_vcf_file(path_vcf: Path, include_filtered: bool) -> List[VCFRecord]:
    """
    Parses a single VCF file.
    :param path_vcf: VCF path
    :param include_filtered: If True, filtered variants are included
    :return: List of SNP positions
    """
    vcf_records = []
    logger.debug(f'Parsing VCF file: {path_vcf}')
    for record in list(vcf.Reader(filename=str(path_vcf))):
        # Remove non-SNP positions
        if not record.is_snp:
            continue
        # Remove filtered records (unless specified otherwise)
        if (not include_filtered) and len(record.FILTER) > 0:
            continue
        vcf_records.append(record)
    logger.info(f"{path_vcf.name} parsed: {len(vcf_records):,} variant positions")
    return vcf_records


def __get_nucleotides_per_position(paths_vcf: List[Path], names: List[str], include_filtered: bool = False) -> \
        Dict[SNPPosition, Dict[str, str]]:
    """
    Returns a dictionary with the nucleotide for each sample at each variant position.
    :param paths_vcf: List of input VCF files
    :param include_filtered: If True, filtered positions are retained
    :return: Sample_names, nucleotide per position per sample
    """
    nucl_by_position = defaultdict(dict)
    for path_vcf, name in zip(paths_vcf, names):
        for record in parse_vcf_file(path_vcf, include_filtered):
            position = SNPPosition(record.CHROM, record.POS, record.REF)
            if len(record.FILTER) > 0:
                nucl_by_position[position][name] = 'N'
            else:
                nucl_by_position[position][name] = str(record.ALT[0])

    # Remove positions with only N
    nucl_by_position = {pos: nucl for pos, nucl in sorted(nucl_by_position.items()) if not all(
        [x == 'N' for x in nucl.values()])}

    logger.info(f"{len(nucl_by_position):,} SNP positions found across all samples")
    return nucl_by_position


def create_snp_matrix(paths_vcf: List[Path], names: List[str], path_out: Path, path_tsv: Optional[Path] = None,
                      include_ref: bool = False) -> None:
    """
    Creates an SNP matrix from the input VCF file.
    :param paths_vcf: List of input VCF files
    :param names: Isolate names
    :param path_out: Output path
    :param path_tsv: TSV file with SNP positions
    :param include_ref: If True, the reference genome is included
    :return: None
    """
    nucl_by_pos = __get_nucleotides_per_position(paths_vcf, names, True)

    # Export TSV file with SNP positions
    if path_tsv is not None:
        data_positions = pd.DataFrame([{'seq_id': pos.contig, 'position': pos.position} for pos in nucl_by_pos.keys()])
        data_positions.to_csv(path_tsv, sep='\t', index=False)
        logger.info(f'SNP positions exported to: {path_tsv}')

    # Create FASTA records
    seq_by_sample_name = {name: [] for name in names}
    if include_ref is True:
        seq_by_sample_name['reference'] = []
    for snp_pos, nucl_by_sample in sorted(nucl_by_pos.items(), key=lambda x: x[0].position):
        for name in seq_by_sample_name.keys():
            seq_by_sample_name[name].append(nucl_by_sample.get(name, snp_pos.reference_base))

    # Write output file
    seqs = [SeqRecord(Seq(''.join(seq)), name, description='') for name, seq in seq_by_sample_name.items()]
    with path_out.open('w') as handle:
        SeqIO.write(seqs, handle, 'fasta')
