import tempfile
import unittest
from importlib.resources import files
from pathlib import Path

from pacu.app.utils import iqtreeutils


class TestIqTreeUtils(unittest.TestCase):
    """
    Tests for the IQ-TREE utils.
    """

    path_fasta = Path(str(files('pacu').joinpath('resources/testdata/snp_matrix.fasta')))

    def test_tree_construction(self) -> None:
        """
        Tests the model selection functionality with regular input.
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            path_out_nwk = Path(dir_, 'iqtree.nwk')
            command = iqtreeutils.run_ml_tree_construction(
                path_fasta=TestIqTreeUtils.path_fasta,
                path_out=path_out_nwk,
                threads=4
            )

            # Verify that the model selection worked
            model = iqtreeutils.extract_selected_model(command.stdout)
            self.assertIsNotNone(model)

            # Verify that the output file is OK
            self.assertTrue(path_out_nwk.exists())
            self.assertGreater(path_out_nwk.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
