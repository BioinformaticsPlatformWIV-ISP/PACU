from pathlib import Path


def is_indexed(path_fasta: Path) -> bool:
    """
    Checks if the FASTA file is indexed for Bowtie2.
    :param path_fasta: Input FASTA file
    :return: True if indexed, False otherwise
    """
    return (path_fasta.parent / f'{path_fasta.name}.1.bt2').exists()
