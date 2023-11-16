import argparse
from pathlib import Path
from typing import Sequence, Optional

from pacu import initialize_logging, Command, logger


class MapToRef(object):
    """
    Main script to map reads to the reference genome.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main script.
        :param args: Arguments (optional)
        :return: None
        """
        self._args = MapToRef._parse_arguments(args)

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        logger.info(f'Starting mapping helper ({self._args.read_type})')
        self._check_dependencies()
        if self._args.read_type == 'illumina':
            path_ref = self._illumina_idx_ref()
            self._map_illumina(path_ref, self._args.output)
        else:
            path_ref = self._ont_idx_ref()
            self._map_ont(path_ref, self._args.output)

    @property
    def ref_name(self) -> str:
        """
        Name of the input reference genome.
        :return: Reference name
        """
        if self._args.ref_fasta_name is None:
            return self._args.ref_fasta.name
        return self._args.ref_fasta_name

    def __symlink_ref_fasta(self, dir_: Path) -> Path:
        """
        Creates a symlink for the reference FASTA file.
        :param dir_: Directory
        :return: Path to symlink
        """
        dir_.mkdir(exist_ok=True, parents=True)
        path_ref_link = dir_ / self.ref_name
        path_ref_link.symlink_to(self._args.ref_fasta)
        logger.debug(f'Creating symlink for reference genome: {path_ref_link} -> {self._args.ref_fasta}')
        return path_ref_link

    def _illumina_idx_ref(self) -> Path:
        """
        Creates an index for the reference genome for Illumina mapping.
        :return: Path to indexed FASTA file
        """
        dir_idx = self._args.dir_working / 'bt2_idx'
        path_ref_link = self.__symlink_ref_fasta(dir_idx)
        command = Command(' '.join(['bowtie2-build', str(path_ref_link.name), str(path_ref_link.name)]))
        command.run(dir_idx)
        if not command.exit_code == 0:
            raise RuntimeError(f'Error creating Bowtie2 index: {command.stderr}')
        logger.info(f'Bowtie2 index created in: {dir_idx}')
        return path_ref_link

    def _map_illumina(self, path_ref: Path, path_out: Path) -> None:
        """
        Maps the Illumina reads
        :param path_ref: Path to reference genome
        :param path_out: Path to output BAM file
        :return: None
        """
        command = Command(' '.join([
            'bowtie2',
            '--end-to-end',
            '--sensitive',
            f'-x {path_ref}',
            f'-1 {self._args.fastq_illumina[0]}',
            f'-2 {self._args.fastq_illumina[1]}',
            f'-p {self._args.threads}',
            f'| samtools view -b | samtools sort - > {path_out}'
        ]))
        command.run(self._args.dir_working)
        if not command.exit_code == 0:
            raise RuntimeError(f'Error mapping reads: {command.stderr}')

    def _ont_idx_ref(self) -> Path:
        """
        Creates an index for the reference genome for ONT mapping.
        :return: Path to indexed FASTA file
        """
        dir_idx = self._args.dir_working / 'mm2_idx'
        path_ref_link = self.__symlink_ref_fasta(dir_idx)
        path_mni = path_ref_link.parent / f'{path_ref_link.name}.mni'
        command = Command(f'minimap2 -x map-ont -d {path_mni} {path_ref_link}')
        command.run(dir_idx)
        if not command.exit_code == 0:
            raise RuntimeError(f'Error creating minimap2 index: {command.stderr}')
        logger.info(f'Minimap2 index created in: {dir_idx}')
        return path_mni

    def _map_ont(self, path_ref: Path, path_out: Path) -> None:
        """
        Maps the ONT reads
        :param path_ref: Path to reference genome
        :param path_out: Path to output BAM file
        :return: None
        """
        command = Command(' '.join([
            'minimap2',
            f'-a {path_ref}',
            str(self._args.fastq_ont),
            f'-t {self._args.threads}',
            f'| samtools view -b | samtools sort - > {path_out}'
        ]))
        command.run(self._args.dir_working)
        if not command.exit_code == 0:
            raise RuntimeError(f'Error mapping reads: {command.stderr}')

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments
        :return: Parsed arguments.
        """
        parser = argparse.ArgumentParser()
        # FASTQ input
        parser.add_argument('--read-type', choices=['ont', 'illumina'], required=True)
        parser.add_argument('--fastq-ont', type=Path)
        parser.add_argument('--fastq-ont-name', type=str, help='Original FASTQ name (used for Galaxy)')
        parser.add_argument('--fastq-illumina', type=Path, nargs=2)
        parser.add_argument(
            '--fastq-illumina-names', type=str, nargs=2, help='Original FASTQ names (used for Galaxy)')

        # FASTA input
        parser.add_argument('--ref-fasta', required=True, help='Reference FASTA file', type=Path)
        parser.add_argument(
            '--ref-fasta-name',  type=Path, help='Original FASTA file name (used for Galaxy)')

        # Other options
        parser.add_argument(
            '--dir-working', type=Path, help='Working directory', default=Path.cwd())
        parser.add_argument('--output', required=True, type=Path, help='Output BAM file')
        parser.add_argument('--threads', type=int, default=4)
        return parser.parse_args(args)

    def _check_dependencies(self) -> None:
        """
        Checks if the required dependencies are available.
        :return: None
        """
        commands = {
            'bowtie2': 'bowtie2 --version',
            'minimap2': 'minimap2 --version',
            'samtools': 'samtools --version',
        }
        logger.info(f'Checking dependencies')
        for tool, command in commands.items():
            command = Command(command)
            command.run(self._args.dir_working, disable_logging=True)
            if not command.exit_code == 0:
                raise RuntimeError(f"Dependency '{tool}' not available")
            logger.info(f"{tool}: OK")


if __name__ == '__main__':
    initialize_logging()
    main = MapToRef()
    main.run()
