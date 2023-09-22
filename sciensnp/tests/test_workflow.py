import tempfile
import unittest
from importlib.resources import files
from pathlib import Path

from sciensnp import ScienSNP


class TestWorkflow(unittest.TestCase):
    """
    Tests for the SNP workflow.
    """

    def test_snp_workflow(self) -> None:
        """
        Tests the SNP workflow with mixed input (Illumina + R10)
        """
        with tempfile.TemporaryDirectory(prefix='sciensnp') as dir_:
            dir_out = Path(dir_, 'output')
            workflow = ScienSNP([
                '--ref-fasta', str(files('sciensnp').joinpath('resources/testdata/NC_002695.2-subset.fasta')),
                '--ilmn-in', str(files('sciensnp').joinpath('resources/testdata/bam/ilmn')),
                '--ont-in', str(files('sciensnp').joinpath('resources/testdata/bam/ont')),
                '--dir-working', str(dir_),
                '--output', str(dir_out),
                '--threads', '8'
            ])
            workflow.run()


if __name__ == '__main__':
    unittest.main()
