import tempfile
import unittest
from importlib.resources import files
from pathlib import Path

from Bio import SeqIO

from sciensnp.app.utils import snpmatrixutils


class TestSNPMatrixUtils(unittest.TestCase):
    """
    Tests for the SNP matrix utils.
    """

    paths_vcf = [p for p in Path(str(files('sciensnp').joinpath('resources/testdata/vcfs'))).glob('*.vcf')]

    def test_create_snp_matrix(self) -> None:
        """
        Tests the create_snp_matrix function.
        """
        with tempfile.TemporaryDirectory(prefix='sciensnp') as dir_:
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


if __name__ == '__main__':
    unittest.main()
