#! /usr/bin/env python
import argparse
from pathlib import Path
from typing import Sequence, Optional, Dict

from pacu import initialize_logging, Command, logger
from pacu.app.utils import workflowutils, trimmingutils, bamutils


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
        self._name = self._determine_name()
        logger.info(f'Input name: {self._name}')

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        logger.info(f'Starting mapping helper ({self._args.read_type})')
        self._check_dependencies()
        self._rename_galaxy_input()
        path_bam_temp = Path(self._args.dir_working, f'{self._name}.bam')

        # Illumina reads
        if self._args.read_type == 'illumina':
            path_ref = self._illumina_idx_ref()
            if self._args.trim:
                fq_dict = self._illumina_trim()
            else:
                fq_dict = {'1P': self._args.fastq_illumina[0], '2P': self._args.fastq_illumina[1]}
            self._illumina_map(fq_dict, path_ref, path_bam_temp)

        # ONT reads
        else:
            path_ref = self._ont_idx_ref()
            path_fq = self._ont_trim() if self._args.trim else self._args.fastq_ont
            self._ont_map(path_fq, path_ref, path_bam_temp)

        # Add a custom tag with the original dataset name (used for Galaxy)
        bamutils.add_custom_tag('PACU_name', self._name, path_bam_temp, self._args.output)

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

    def _illumina_trim(self) -> Dict[str, Path]:
        """
        Trims the Illumina reads.
        :return: None
        """
        logger.info(f'Trimming Illumina reads')

        # Create directory
        dir_trim = Path(self._args.dir_working, 'trim')
        dir_trim.mkdir(exist_ok=True, parents=True)

        # Run trimming
        command = Command(' '.join([
            'trimmomatic', 'PE',
            '-baseout', f'{self._name}.fastq.gz',
            f'-threads {self._args.threads}',
            str(self._args.fastq_illumina[0]),
            str(self._args.fastq_illumina[1]),
            f"ILLUMINACLIP:{Path(trimmingutils.trimmomatic_dir_adapters(), 'NexteraPE-PE.fa')}:2:30:10",
            'LEADING:10',
            'TRAILING:10',
            'SLIDINGWINDOW:4:20',
            'MINLEN:40',
            '-phred33'
        ]))
        command.run(dir_trim)
        if not command.exit_code == 0:
            raise RuntimeError(f'Error trimming reads: {command.stderr}')
        return trimmingutils.trimmomatic_collect_output(dir_trim)

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

    def _illumina_map(self, fq_dict: Dict[str, Path], path_ref: Path, path_out: Path) -> None:
        """
        Maps the Illumina reads to the reference genome.
        :param fq_dict: Input FASTQ dictionary
        :param path_ref: Path to reference genome
        :param path_out: Path to output BAM file
        :return: None
        """
        # Construct bowtie2 command
        parts_bt2 = [
            'bowtie2', '--end-to-end', '--sensitive', f'-x {path_ref}', f"-1 {fq_dict['1P']}", f"-2 {fq_dict['2P']}"]

        # Add orphaned reads (if available)
        orphaned_reads = [fq_dict.get(k) for k in ('1U', '2U') if fq_dict.get(k) is not None]
        if len(orphaned_reads) > 0:
            parts_bt2.append(f"-U {','.join(str(p) for p in orphaned_reads)}")
        parts_bt2.append(f'-p {self._args.threads}')

        # Execute the command
        command = Command(' '.join([
            *parts_bt2,
            f'| samtools view -b | samtools sort - > {path_out}'
        ]))
        command.run(self._args.dir_working)
        if not command.exit_code == 0:
            raise RuntimeError(f'Error mapping reads: {command.stderr}')

    def _ont_trim(self, min_qual: int = 7, min_len: int = 500) -> Path:
        """
        Trims the input ONT reads.
        :param min_qual: Minimum read quality
        :param min_len: Minimum read length
        :return: Path to trimmed reads
        """
        logger.info(f'Trimming ONT reads')

        # Create directory
        dir_trim = Path(self._args.dir_working, 'trim')
        dir_trim.mkdir(exist_ok=True, parents=True)

        # Run the command
        path_out = dir_trim / f'{self._name}.fastq.gz'
        command = Command(' '.join([
            'seqkit seq',
            '--min-qual', str(min_qual),
            '--min-len', str(min_len),
            str(self._args.fastq_ont),
            f'> {path_out}'
        ]))
        command.run(path_out.parent)
        if not command.exit_code == 0:
            raise RuntimeError(f'Error trimming reads: {command.stderr}')
        return path_out

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

    def _ont_map(self, path_fq: Path, path_ref: Path, path_out: Path) -> None:
        """
        Maps the ONT reads
        :param path_fq: Path to the input FASTQ file
        :param path_ref: Path to reference genome
        :param path_out: Path to output BAM file
        :return: None
        """
        command = Command(' '.join([
            'minimap2',
            f'-a {path_ref}',
            str(path_fq),
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
        parser.add_argument('--trim', action='store_true', help='Trim reads prior to mapping')
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
        if self._args.trim:
            commands['trimmomatic'] = 'trimmomatic -version'
            commands['seqkit'] = 'seqkit version'
        logger.info(f'Checking dependencies')
        for tool, command_str in commands.items():
            command = Command(command_str)
            command.run(self._args.dir_working, disable_logging=True)
            if not command.exit_code == 0:
                raise RuntimeError(f"Dependency '{tool}' not available")
            logger.info(f"{tool}: OK")

    def _rename_galaxy_input(self) -> None:
        """
        Renames the input files when the script is executed through Galaxy.
        :return: None
        """
        if self._args.fastq_illumina_names is not None:
            # Symlink the input files
            Path(self._args.dir_working, 'input').mkdir(parents=True, exist_ok=True)
            for idx in (0, 1):
                is_gzipped = workflowutils.is_gzipped(self._args.fastq_illumina[idx])
                path_link = Path(
                    self._args.dir_working, 'input', f'{self._name}_{idx+1}P.fastq' + ('.gz' if is_gzipped else ''))
                logger.debug(f'Creating link: {path_link} -> {self._args.fastq_illumina[idx]}')
                path_link.absolute().symlink_to(self._args.fastq_illumina[idx].absolute())
                self._args.fastq_illumina[idx] = path_link
        elif self._args.fastq_ont_name is not None:
            Path(self._args.dir_working, 'input').mkdir(parents=True, exist_ok=True)
            is_gzipped = workflowutils.is_gzipped(self._args.fastq_ont)
            path_link = Path(
                self._args.dir_working, 'input', f'{self._name}.fastq' + ('.gz' if is_gzipped else '')).absolute()
            path_link.symlink_to(self._args.fastq_ont.absolute())
        else:
            logger.debug(f'Not renaming inputs')

    def _determine_name(self) -> str:
        """
        Determines the input name.
        :return: Input name
        """
        if self._args.read_type == 'ont':
            if self._args.fastq_ont_name is not None:
                return workflowutils.determine_name_from_fq(fq_ont=Path(self._args.fastq_ont_name))
            return workflowutils.determine_name_from_fq(fq_ont=self._args.fastq_ont)
        else:
            if self._args.fastq_illumina_names is not None:
                return workflowutils.determine_name_from_fq(fq_illumina_1p=Path(self._args.fastq_illumina_names[0]))
            return workflowutils.determine_name_from_fq(fq_illumina_1p=self._args.fastq_illumina[0])

def main() -> None:
    """
    Wrapper around the main script.
    :return: None
    """
    initialize_logging()
    map_to_ref = MapToRef()
    map_to_ref.run()


if __name__ == '__main__':
    main()
