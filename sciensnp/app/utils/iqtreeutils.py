import unittest
from pathlib import Path

from sciensnp import Command


def run_ml_tree_construction(path_fasta: Path, path_out: Path, threads: int = 4) -> None:
    """
    Runs the ML tree reconstruction using IQ-TREE.
    :param path_fasta: Input SNP matrix FASTA file
    :param path_out: Output path
    :param threads: Number of threads to use
    :return: None
    """
    if not path_out.name.endswith('.nwk'):
        raise ValueError("Output path should end with the '.nwk' extension")
    basename = path_out.name.replace('.nwk', '')
    prefix = path_out.parent / basename
    command = Command(f"iqtree2 -s {path_fasta} --boot 100 -m MFP --prefix {prefix} -T {threads}")
    command.run(path_out.parent)
    if not command.exit_code == 0:
        raise RuntimeError(f'Error running IQ-TREE: {command.stderr}')
    Path(path_out.parent, f'{basename}.treefile').rename(path_out)


if __name__ == '__main__':
    unittest.main()
