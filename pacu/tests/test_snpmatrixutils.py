import tempfile
import unittest
from importlib.resources import files
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from pacu.app.utils import snpmatrixutils


class TestSNPMatrixUtils(unittest.TestCase):
    """
    Tests for the SNP matrix utils.
    """

    paths_vcf = [p for p in Path(str(files('pacu').joinpath('resources/testdata/vcfs'))).glob('*.vcf')]

    def test_create_snp_matrix(self) -> None:
        """
        Tests the create_snp_matrix function.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            # Create SNP matrix
            path_out = Path(dir_, 'snp_matrix.fasta')
            snpmatrixutils.create_snp_matrix(
                TestSNPMatrixUtils.paths_vcf,
                [p.name for p in TestSNPMatrixUtils.paths_vcf],
                path_out)

            # Verify output file
            with path_out.open() as handle:
                seqs = list(SeqIO.parse(handle, 'fasta'))
                self.assertTrue(all(len(s) > 0 for s in seqs))
                self.assertTrue(all(len(s) == len(seqs[0]) for s in seqs))

    def test_create_snp_matrix_with_tsv(self) -> None:
        """
        Tests the create_snp_matrix function.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            path_tsv_out = Path(dir_, 'snp_positions.tsv')

            # Create SNP matrix
            path_out = Path(dir_, 'snp_matrix.fasta')
            snpmatrixutils.create_snp_matrix(
                TestSNPMatrixUtils.paths_vcf,
                [p.name for p in TestSNPMatrixUtils.paths_vcf],
                path_out,
                path_tsv=path_tsv_out
            )

            # Verify output SNP matrix
            with path_out.open() as handle:
                seqs = list(SeqIO.parse(handle, 'fasta'))
                self.assertTrue(all(len(s) > 0 for s in seqs), 'SNP matrix contains empty sequence(s)')
                self.assertTrue(
                    all(len(s) == len(seqs[0]) for s in seqs),
                    'SNP matrix contains sequences of varying lengths')
                sequence_length = len(seqs[0])

            # Verify output SNP positions TSV file
            self.assertTrue(path_tsv_out.exists(), 'SNP positions TSV file not generated')
            data_snp_positions = pd.read_table(path_tsv_out)
            self.assertEqual(len(data_snp_positions), sequence_length)


if __name__ == '__main__':
    unittest.main()
