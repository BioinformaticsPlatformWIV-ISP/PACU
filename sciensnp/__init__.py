import argparse
from importlib.resources import files
from pathlib import Path
from typing import Optional, Sequence

from sciensnp.app.command import Command
from sciensnp.app.utils import snakemakeutils
from sciensnp.app.utils.loggingutils import initialize_logging, logger
from .version import __version__


class ScienSNP(object):
    """
    Main class to run the ScienSNP pipeline.
    """

    def __init__(self, args: Optional[Sequence[str]] = None) -> None:
        """
        Initializes the main class.
        :param args: Arguments (optional)
        :return: None
        """
        self._args = ScienSNP._parse_arguments(args)

    def run(self) -> None:
        """
        Runs the main script.
        :return: None
        """
        logger.info('Starting ScienSNP workflow')
        self._check_dependencies()
        path_snakefile = Path(str(files('sciensnp').joinpath('resources/snp_workflow.smk')))
        targets = [self._args.output.absolute() / 'report.html']
        if self._args.ref_bed is None:
            path_bed = self._args.dir_working / 'empty_phages.bed'
            path_bed.touch()
            logger.info(f'Creating empty phage BED file: {path_bed}')
            self._args.ref_bed = path_bed
        path_config = self.__create_config_file()
        snakemakeutils.run_snakemake(path_snakefile, path_config, targets, self._args.dir_working, self._args.threads)
        logger.info(f'Workflow finished successfully, output available in: {self._args.output}')

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
        parser.add_argument('--ref-fasta', required=True, help='Reference FASTA file', type=Path)
        parser.add_argument('--ref-bed', type=Path, help='BED file with phage regions')
        parser.add_argument(
            '--dir-working', type=Path, help='Working directory', default=Path.cwd())
        parser.add_argument('--output', required=True, type=Path, help='Output directory')

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
            '--version', help='Print version and exit', action='version', version=f'ScienSNP {__version__}')

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
            'include_ref': self._args.include_ref,
            'input': input_dict,
            'phylogeny_method': 'mega' if self._args.use_mega else 'iqtree',
            'output': {
                'html': str(self._args.output.absolute() / 'report.html'),
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
            'gubbins': 'ml gubbins; gubbins -h',
            'snp-dists': 'snp-dists -v',
            'figtree': 'figtree -help',

        }
        if self._args.use_mega:
            commands['MEGA'] = 'megacc --version'
        else:
            commands['IQ-TREE'] = 'iqtree2 --version'
        logger.info(f'Checking dependencies')
        for tool, command in commands.items():
            command = Command(command)
            command.run(self._args.dir_working, disable_logging=True)
            if not command.exit_code == 0:
                raise RuntimeError(f"Dependency '{tool}' not available")
            logger.info(f"{tool}: OK")
