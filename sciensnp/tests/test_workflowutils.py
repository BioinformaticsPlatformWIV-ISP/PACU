import tempfile
import unittest
from importlib.resources import files
from pathlib import Path

from sciensnp.app.utils import workflowutils
from sciensnp.config import config


class TestWorkflowUtils(unittest.TestCase):
    """
    Tests for the workflow utils.
    """

    def test_plot_newick_phylogeny(self) -> None:
        """
        Tests the plot_newick_phylogeny function.
        """
        path_nwk = Path(str(files('sciensnp').joinpath('resources/testdata/phylogeny.nwk')))
        with tempfile.TemporaryDirectory(prefix='mega_', dir=config['dir_temp']) as dir_:
            # Create visualization
            path_out = Path(dir_, 'tree.png')
            workflowutils.plot_newick_phylogeny(path_nwk, path_out)

            # Verify output file
            self.assertTrue(path_out.exists())
            self.assertGreater(path_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
