import argparse
import shutil
from importlib.resources import files
from pathlib import Path
from typing import Optional, Sequence

from Bio import SeqIO

from pacu.app.command import Command
from pacu.app.utils import snakemakeutils
from pacu.app.utils.loggingutils import logger
from .app.utils.cliutils import path_to_absolute
from .version import __version__


class PACU(object):
    """
    Main class to run the PACU pipeline.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        :return: None
        """
        self._args = PACU._parse_arguments(args)
        self._path_html_out = self._args.output.absolute() / (
            self._args.output_html if self._args.output_html else 'report.html')
        if self._path_html_out.exists():
            self._path_html_out.unlink()

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        logger.info('Starting PACU workflow')
        self._check_dependencies()
        self._validate_input_files()
        path_snakefile = Path(str(files('pacu').joinpath('resources/snp_workflow.smk')))
        if self._args.ref_bed is None:
            path_bed = self._args.dir_working / 'empty_phages.bed'
            path_bed.touch()
            logger.info(f'Creating empty phage BED file: {path_bed}')
            self._args.ref_bed = path_bed
        # Run Snakemake workflow
        path_config = self.__create_config_file()
        snakemakeutils.run_snakemake(
            path_snakefile, path_config, [self._path_html_out], self._args.dir_working, self._args.threads)
        logger.info(f'Workflow finished successfully, output available in: {self._args.output}')

        # Copy the log file if it exists
        if (self._args.dir_working / 'pacu.log').exists():
            shutil.copyfile(self._args.dir_working / 'pacu.log', self._args.output / 'pacu.log')

    @staticmethod
    def _parse_arguments(args: Optional[Sequence[str]]) -> argparse.Namespace:
        """
        Parses the command line arguments.
        :param args: Arguments
        :return: Parsed arguments.
        """
        parser = argparse.ArgumentParser()

        # Input and output
        parser.add_argument('--ilmn-in', help='Directory with Illumina input BAM files')
        parser.add_argument('--ont-in', help='Directory with ONT input BAM files')
        parser.add_argument('--ref-fasta', required=True, help='Reference FASTA file', type=path_to_absolute)
        parser.add_argument('--ref-bed', type=path_to_absolute, help='BED file with phage regions')
        parser.add_argument(
            '--dir-working', type=path_to_absolute, help='Working directory', default=Path.cwd())
        parser.add_argument('--output', required=True, type=path_to_absolute, help='Output directory')
        parser.add_argument('--output-html', type=path_to_absolute, help='Output report name')

        # Parameters
        parser.add_argument(
            '--use-mega', action='store_true',
            help='If set, MEGA is used for the construction of the phylogeny (instead of IQ-TREE)')
        parser.add_argument(
            '--include-ref', action='store_true', help='If set, the reference genome is included in the phylogeny')
        parser.add_argument('--min-snp-af', type=float, default=0.90, help='Minimum allele frequency for variants')
        parser.add_argument('--min-snp-qual', type=int, default=25, help='Minimum SNP quality')
        parser.add_argument('--min-snp-depth', type=int, default=5, help='Minimum SNP depth')
        parser.add_argument('--min-snp-dist', type=int, default=10, help='Minimum distance between SNPs')
        parser.add_argument('--skip-gubbins', action='store_true', help='If set, gubbins is skipped')
        # parser.add_argument('--min-sb-pv', type=float, default=0, help='Minimum P-value for strand bias')
        parser.add_argument(
            '--min-global-depth', type=int, default=5,
            help='Minimum depth for all samples to include positions in SNP analysis')
        parser.add_argument('--min-mq-depth', type=int, default=5, help='MQ cutoff for samtools depth')
        parser.add_argument(
            '--bcftools-filt-af1', action='store_true',
            help='If enabled, allele frequency filtering also considers the VAF value')

        # Image
        parser.add_argument('--image-width', type=int, default=480, help='Image width')
        parser.add_argument('--image-height', type=int, default=480, help='Image height')

        # Threads
        parser.add_argument('--threads', type=int, default=4)
        parser.add_argument(
            '--version', help='Print version and exit', action='version', version=f'PACU {__version__}')

        return parser.parse_args(args)

    def __create_config_file(self) -> Path:
        """
        Creates the config file for the Snakemake workflow.
        :return: Path to config file
        """
        input_dict = {}
        if self._args.ilmn_in is not None:
            for bam_file in sorted(Path(self._args.ilmn_in).glob('*.bam')):
                input_dict[bam_file.stem] = {'bam': str(bam_file.absolute()), 'tech': 'ilmn'}
        if self._args.ont_in is not None:
            for bam_file in sorted(Path(self._args.ont_in).glob('*.bam')):
                input_dict[bam_file.stem] = {'bam': str(bam_file.absolute()), 'tech': 'ont'}
        if len(input_dict) == 0:
            raise FileNotFoundError(
                f'No input datasets found (searched in Illumina={self._args.ilmn_in}, ONT={self._args.ont_in})')
        if len(input_dict) < 4:
            raise ValueError(
                f'At least four input BAM files are required (for bootstrapping), {len(input_dict)} found.')

        config_data = {
            'depth': {
                'min_mq': self._args.min_mq_depth,
                'min_depth': 5
            },
            'filters': {
                'min_af': self._args.min_snp_af,
                'min_depth': self._args.min_snp_depth,
                'min_qual': self._args.min_snp_qual,
                'min_dist': self._args.min_snp_dist
            },
            'image': {
                'width': self._args.image_width,
                'height': self._args.image_height
            },
            'skip_gubbins': self._args.skip_gubbins,
            'include_ref': self._args.include_ref,
            'input': input_dict,
            'phylogeny_method': 'mega' if self._args.use_mega else 'iqtree',
            'output': {
                'html': str(self._path_html_out),
                'dir': str(self._args.output.absolute())
            },
            'reference': {
                'fasta': str(self._args.ref_fasta.absolute()),
                'bed_phages': str(self._args.ref_bed.absolute()),
            }
        }
        return Path(snakemakeutils.generate_config_file(config_data, self._args.dir_working.absolute()))

    def _check_dependencies(self) -> None:
        """
        Checks if the required dependencies are available.
        :return: None
        """
        commands = {
            'bcftools': 'bcftools --version',
            'samtools': 'samtools --version',
            'bedtools': 'bedtools --version',
            'gubbins': 'gubbins -h',
            'snp-dists': 'snp-dists -v',
            'figtree': 'figtree -help',
        }
        if self._args.skip_gubbins:
            commands.pop('gubbins')
        if self._args.use_mega:
            commands['MEGA'] = 'megacc --version'
        else:
            commands['IQ-TREE'] = 'iqtree2 --version'
        logger.info('Checking dependencies')
        for tool, command in commands.items():
            command = Command(command)
            command.run(self._args.dir_working, disable_logging=True)
            if not command.exit_code == 0:
                raise RuntimeError(f"Dependency '{tool}' not available")
            logger.info(f"{tool}: OK")

    def _validate_input_files(self) -> None:
        """
        Checks if the provided input files are valid.
        :return: None
        """
        if self._args.skip_gubbins:
            return
        with open(self._args.ref_fasta) as handle:
            seqs = list(SeqIO.parse(handle, "fasta"))
            if len(seqs) > 1:
                raise ValueError(
                    f"Reference genome must contain exactly 1 sequence "
                    f'(found {len(seqs):,}). If your reference is multi-contig, run with "--skip-gubbins."'
                )
