import tempfile
import unittest
from importlib.resources import files
from pathlib import Path

from sciensnp.map_to_ref import MapToRef


class TestMapToRef(unittest.TestCase):
    """
    Tests the map to ref script.
    """

    def test_map_to_ref_illumina(self) -> None:
        """
        Tests the map to ref script with Illumina data.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='sciensnp') as dir_:
            path_out = Path(dir_, 'mapped_reads.bam')
            map_to_ref = MapToRef([
                '--ref-fasta', str(files('sciensnp').joinpath('resources/testdata/NC_002695.2-subset.fasta')),
                '--data-type', 'illumina',
                '--fastq-illumina',
                str(files('sciensnp').joinpath('resources/testdata/fastq/TIAC1151_1P.fastq.gz')),
                str(files('sciensnp').joinpath('resources/testdata/fastq/TIAC1151_2P.fastq.gz')),
                '--output', str(path_out),
                '--dir-working', dir_,
                '--threads', '4'
            ])
            map_to_ref.run()

            # Verify the output file
            self.assertTrue(path_out.exists())
            self.assertGreater(path_out.stat().st_size, 0)

    def test_map_to_ref_ont(self) -> None:
        """
        Tests the map to ref script with ONT data.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='sciensnp') as dir_:
            path_out = Path(dir_, 'mapped_reads.bam')
            map_to_ref = MapToRef([
                '--ref-fasta', str(files('sciensnp').joinpath('resources/testdata/NC_002695.2-subset.fasta')),
                '--data-type', 'ont',
                '--fastq-ont',
                str(files('sciensnp').joinpath('resources/testdata/fastq/TIAC1151-ont.gz')),
                '--output', str(path_out),
                '--dir-working', dir_,
                '--threads', '4'
            ])
            map_to_ref.run()

            # Verify the output file
            self.assertTrue(path_out.exists())
            self.assertGreater(path_out.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
