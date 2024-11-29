import tempfile
import unittest
from importlib.resources import files
from pathlib import Path

from pacu.app.utils import workflowutils


class TestWorkflowUtils(unittest.TestCase):
    """
    Tests for the workflow utils.
    """

    def test_plot_newick_phylogeny(self) -> None:
        """
        Tests the plot_newick_phylogeny function.
        """
        path_nwk = Path(str(files('pacu').joinpath('resources/testdata/phylogeny.nwk')))
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            # Create visualization
            path_out = Path(dir_, 'tree.png')
            workflowutils.plot_newick_phylogeny(path_nwk, path_out)

            # Verify output file
            self.assertTrue(path_out.exists())
            self.assertGreater(path_out.stat().st_size, 0)

    def test_determine_name_from_fq_ont(self) -> None:
        """
        Tests the 'determine_name_from_fq' function with ONT input.
        :return: None
        """
        self.assertEqual(workflowutils.determine_name_from_fq(fq_ont=Path('/path/to/reads.fastq')), 'reads')
        self.assertEqual(workflowutils.determine_name_from_fq(fq_ont=Path('/path/to/reads.fastq.gz')), 'reads')
        self.assertEqual(workflowutils.determine_name_from_fq(fq_ont=Path('/path/to/reads.fq')), 'reads')
        self.assertEqual(workflowutils.determine_name_from_fq(fq_ont=Path('/path/to/reads.fq.gz')), 'reads')

    def test_determine_name_from_fq_illumina(self) -> None:
        """
        Tests the 'determine_name_from_fq' function with Illumina input.
        :return: None
        """
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/my-reads_1P.fastq')), 'my-reads')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/input_file_1.fastq.gz')), 'input_file')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/reads_1.fq')), 'reads')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/reads_1P.fq.gz')), 'reads')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/reads_S22_L001_R1_001.fastq.gz')), 'reads_S22')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/reads_S22_L001_R1_001.fq.gz')), 'reads_S22')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/unknown_pattern.fastq')), 'unknown_pattern')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/unknown_pattern.fastq.gz')), 'unknown_pattern')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/reads_R1.fq')), 'reads')
        self.assertEqual(workflowutils.determine_name_from_fq(
            fq_illumina_1p=Path('/path/to/reads_R1.fastq.gz')), 'reads')


if __name__ == '__main__':
    unittest.main()
