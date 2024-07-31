import tempfile
import unittest
from importlib.resources import files
from pathlib import Path

from pacu import initialize_logging
from pacu.map_to_ref import MapToRef


class TestMapToRef(unittest.TestCase):
    """
    Tests the map to ref script.
    """

    def test_map_to_ref_illumina(self) -> None:
        """
        Tests the map to ref script with Illumina data.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            path_out = Path(dir_, 'mapped_reads.bam')
            map_to_ref = MapToRef([
                '--ref-fasta', str(files('pacu').joinpath('resources/testdata/NC_002695.2-subset.fasta')),
                '--read-type', 'illumina',
                '--fastq-illumina',
                str(files('pacu').joinpath('resources/testdata/fastq/TIAC1151_1P.fastq.gz')),
                str(files('pacu').joinpath('resources/testdata/fastq/TIAC1151_2P.fastq.gz')),
                '--output', str(path_out),
                '--dir-working', dir_,
                '--threads', '4'
            ])
            map_to_ref.run()

            # Verify the output file
            self.assertTrue(path_out.exists())
            self.assertGreater(path_out.stat().st_size, 0)

    def test_map_to_ref_illumina_with_trim(self) -> None:
        """
        Tests the map to ref script with Illumina data and trimming enabled.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            path_out = Path(dir_, 'mapped_reads.bam')
            map_to_ref = MapToRef([
                '--ref-fasta', str(files('pacu').joinpath('resources/testdata/NC_002695.2-subset.fasta')),
                '--read-type', 'illumina',
                '--fastq-illumina',
                str(files('pacu').joinpath('resources/testdata/fastq/TIAC1151_1P.fastq.gz')),
                str(files('pacu').joinpath('resources/testdata/fastq/TIAC1151_2P.fastq.gz')),
                '--output', str(path_out),
                '--trim',
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
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            path_out = Path(dir_, 'mapped_reads.bam')
            map_to_ref = MapToRef([
                '--ref-fasta', str(files('pacu').joinpath('resources/testdata/NC_002695.2-subset.fasta')),
                '--read-type', 'ont',
                '--fastq-ont',
                str(files('pacu').joinpath('resources/testdata/fastq/TIAC1151-ont.fastq.gz')),
                '--output', str(path_out),
                '--dir-working', dir_,
                '--threads', '4'
            ])
            map_to_ref.run()

            # Verify the output file
            self.assertTrue(path_out.exists())
            self.assertGreater(path_out.stat().st_size, 0)

    def test_map_to_ref_ont_with_trim(self) -> None:
        """
        Tests the map to ref script with ONT data and trimming enabled.
        :return: None
        """
        with tempfile.TemporaryDirectory(prefix='pacu') as dir_:
            path_out = Path(dir_, 'mapped_reads.bam')
            map_to_ref = MapToRef([
                '--ref-fasta', str(files('pacu').joinpath('resources/testdata/NC_002695.2-subset.fasta')),
                '--read-type', 'ont',
                '--fastq-ont',
                str(files('pacu').joinpath('resources/testdata/fastq/TIAC1151-ont.fastq.gz')),
                '--trim',
                '--output', str(path_out),
                '--dir-working', dir_,
                '--threads', '4'
            ])
            map_to_ref.run()

            # Verify the output file
            self.assertTrue(path_out.exists())
            self.assertGreater(path_out.stat().st_size, 0)


if __name__ == '__main__':
    initialize_logging()
    unittest.main()
