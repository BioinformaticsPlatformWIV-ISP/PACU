from pathlib import Path


def is_bt2_indexed(path_fasta: Path) -> bool:
    """
    Checks if the FASTA file is indexed for Bowtie2.
    :param path_fasta: Input FASTA file
    :return: True if indexed, False otherwise
    """
    return (path_fasta.parent / f'{path_fasta.name}.1.bt2').exists()

def is_mm2_indexed(path_fasta: Path) -> bool:
    """
    Checks if the FASTA file is indexed for Bowtie2.
    :param path_fasta: Input FASTA file
    :return: True if indexed, False otherwise
    """
    return (path_fasta.parent / f'{path_fasta.name}.mni').exists()
