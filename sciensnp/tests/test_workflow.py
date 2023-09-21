import unittest
from pathlib import Path

from sciensnp.run_sciensnp import ScienSNP


class TestWorkflow(unittest.TestCase):
    """
    Tests for the SNP workflow.
    """

    DIR_TEST = Path('/scratch/temp/snp_wf_in')

    def test_snp_workflow(self) -> None:
        """
        Tests the SNP workflow with mixed input (Illumina + R10)
        """
        workflow = ScienSNP([
            '--ref-fasta', '/db/refgenomes/Escherichia_coli/NC_002695.2.fasta',
            '--ilmn-in', str(TestWorkflow.DIR_TEST / 'ilmn'),
            '--ont-in', str(TestWorkflow.DIR_TEST / 'r10'),
            '--dir-working', '/scratch/temp/scien_snp/work',
            '--output', '/scratch/temp/'
        ])
        workflow.run()


if __name__ == '__main__':
    unittest.main()
