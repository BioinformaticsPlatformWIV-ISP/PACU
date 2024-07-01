import tempfile
import unittest
from importlib.resources import files
from pathlib import Path

from pacu import PACU


class TestWorkflow(unittest.TestCase):
    """
    Tests for the SNP workflow.
    """

    def test_snp_workflow_iqtree(self) -> None:
        """
        Tests the SNP workflow with mixed input (Illumina + R10) and IQ-TREE tree construction.
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            dir_out = Path(dir_, 'output')
            workflow = PACU([
                '--ref-fasta', str(files('pacu').joinpath('resources/testdata/NC_002695.2-subset.fasta')),
                '--ilmn-in', str(files('pacu').joinpath('resources/testdata/bam/ilmn')),
                '--ont-in', str(files('pacu').joinpath('resources/testdata/bam/ont')),
                '--dir-working', str(dir_),
                '--output', str(dir_out),
                '--threads', '8'
            ])
            workflow.run()

    def test_snp_workflow_mega(self) -> None:
        """
        Tests the SNP workflow with mixed input (Illumina + R10) and MEGA tree construction.
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            dir_out = Path(dir_, 'output')
            workflow = PACU([
                '--ref-fasta', str(files('pacu').joinpath('resources/testdata/NC_002695.2-subset.fasta')),
                '--ilmn-in', str(files('pacu').joinpath('resources/testdata/bam/ilmn')),
                '--ont-in', str(files('pacu').joinpath('resources/testdata/bam/ont')),
                '--use-mega',
                '--dir-working', str(dir_),
                '--output', str(dir_out),
                '--threads', '8'
            ])
            workflow.run()

    def test_snp_workflow_no_gubbins(self) -> None:
        """
        Tests the SNP workflow with mixed input (Illumina + R10) without the 'gubbins' step.
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            dir_ = Path('/scratch/temp/tmp-wf')
            dir_out = Path(dir_, 'output')
            workflow = PACU([
                '--ref-fasta', str(files('pacu').joinpath('resources/testdata/NC_002695.2-subset.fasta')),
                '--ilmn-in', str(files('pacu').joinpath('resources/testdata/bam/ilmn')),
                '--ont-in', str(files('pacu').joinpath('resources/testdata/bam/ont')),
                '--skip-gubbins',
                '--dir-working', str(dir_),
                '--output', str(dir_out),
                '--threads', '8'
            ])
            workflow.run()


if __name__ == '__main__':
    unittest.main()
