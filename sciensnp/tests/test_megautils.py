import tempfile
import unittest
from importlib.resources import files
from pathlib import Path

from sciensnp.app.utils import megautils


class TestMEGAUtils(unittest.TestCase):
    """
    Tests for the MEGA utils.
    """

    path_fasta = Path(str(files('sciensnp').joinpath('resources/testdata/snp_matrix.fasta')))

    def test_model_selection(self) -> None:
        """
        Tests the model selection functionality with regular input.
        """
        with tempfile.TemporaryDirectory(prefix='sciensnp') as dir_:
            path_out = Path(dir_, 'model_selection.csv')
            megautils.run_model_selection(
                path_fasta=TestMEGAUtils.path_fasta,
                path_out=path_out,
                dir_=Path(dir_),
                branch_swap_filter='Very weak',
                missing_data_treatment='Partial deletion',
                site_cov_cutoff=50,
                threads=4
            )

            # Parse the output file
            model_info = megautils.parse_model_selection_csv(path_out)
            self.assertIn('model', model_info)
            self.assertIn('rates_among_sites', model_info)

    def test_tree_building(self) -> None:
        """
        Tests the tree building functionality with regular input.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='sciensnp') as dir_:
            path_model_selection_csv = Path(str(files('sciensnp').joinpath(
                'resources/testdata/mega/model_selection_output.csv')))

            # Run tree building
            path_out_nwk = Path(dir_, 'tree.nwk')
            megautils.run_tree_building(
                path_fasta=TestMEGAUtils.path_fasta,
                path_csv=path_model_selection_csv,
                path_out=path_out_nwk,
                dir_=Path(dir_),
                branch_swap_filter='Weak',
                missing_data_treatment='Use all sites',
                site_cov_cutoff=None,
                bootstrap_replicates=100,
                heuristic_method='SPR3',
                initial_tree='NJ',
                gamma_categories=5,
                threads=4
            )

            # Verify that the output file is OK
            self.assertTrue(path_out_nwk.exists())
            self.assertGreater(path_out_nwk.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
