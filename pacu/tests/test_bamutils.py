import tempfile
import unittest
from importlib.resources import files
from pathlib import Path

from pacu.app.utils import bamutils
from pacu.app.utils.loggingutils import initialize_logging


class TestBAMUtils(unittest.TestCase):
    """
    Tests the BAM utils module.
    """

    def test_add_custom_tag(self) -> None:
        """
        Tests add_custom_tag function.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            path_bam_in = Path(str(files('pacu').joinpath('resources/testdata/bam/ilmn/TIAC1151.bam')))
            path_bam_out = Path(dir_, 'updated.bam')
            bamutils.add_custom_tag('PACU_name', 'custom_sample', path_bam_in, path_bam_out)
            self.assertTrue(path_bam_out.exists())
            self.assertEqual(bamutils.read_custom_tag(path_bam_out, 'PACU_name'), 'custom_sample')


if __name__ == '__main__':
    initialize_logging()
    unittest.main()
