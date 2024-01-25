#!/usr/bin/env python
import argparse
from pathlib import Path
from typing import Tuple, List

from pacu import initialize_logging, logger, PACU
from pacu.app.utils import workflowutils, bamutils


def parse_galaxy_args() -> Tuple[argparse.Namespace, List[str]]:
    """
    Parses the Galaxy arguments.
    :return: Parsed arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref-fasta', required=True, type=Path)
    parser.add_argument('--ref-fasta-name', type=str)
    parser.add_argument('--bam-ilmn', nargs=2, action='append')
    parser.add_argument('--bam-ont', nargs=2, action='append')
    parser.add_argument('--dir-working', default=Path.cwd(), type=Path, help='Working directory')
    return parser.parse_known_args()


def _extract_dataset_name(bam_in: Path, name_orig: str) -> str:
    """
    Extracts the dataset name from the input BAM file.
    :param bam_in: Input BAM file
    :param name_orig: Galaxy name
    :return: Extracted name
    """
    try:
        return bamutils.read_custom_tag(bam_in, 'PACU_name')
    except ValueError:
        logger.debug(f"'PACU_name' not found in BAM header: {name_orig}")
    return workflowutils.sanitize_bam_input(name_orig)


if __name__ == '__main__':
    initialize_logging()
    logger.info(f'Running PACU through Galaxy')
    args, unparsed_args = parse_galaxy_args()

    # Create input directory
    dir_in = args.dir_working / 'input'
    dir_in.mkdir(exist_ok=True)

    # Symlink the reference genome
    path_ref_link = dir_in / args.ref_fasta_name
    if not path_ref_link.exists():
        path_ref_link.symlink_to(Path(args.ref_fasta))
    unparsed_args.extend(['--ref-fasta', str(path_ref_link.absolute())])

    # Symlink the BAM files (Illumina)
    if args.bam_ilmn is not None:
        logger.info(f'Creating symlinks for Illumina input ({len(args.bam_ilmn)} datasets)')
        dir_bam = dir_in / 'bam_ilmn'
        dir_bam.mkdir(exist_ok=True)
        for path_bam, galaxy_name in args.bam_ilmn:
            try:
                basename = f'{_extract_dataset_name(Path(path_bam), galaxy_name)}.bam'
                (dir_bam / basename).symlink_to(path_bam)
            except FileExistsError:
                logger.debug(f"Symlink for '{galaxy_name}' already exists")
        unparsed_args.extend(['--ilmn-in', str(dir_bam.absolute())])

    # Symlink the BAM files (ONT)
    if args.bam_ont is not None:
        logger.info(f'Creating symlinks for ONT input ({len(args.bam_ont)} datasets)')
        dir_bam = dir_in / 'bam_ont'
        dir_bam.mkdir(exist_ok=True)
        for path_bam, galaxy_name in args.bam_ont:
            try:
                basename = f'{_extract_dataset_name(Path(path_bam), galaxy_name)}.bam'
                (dir_bam / basename).symlink_to(path_bam)
            except FileExistsError:
                logger.debug(f"Symlink for '{galaxy_name}' already exists")
        unparsed_args.extend(['--ont-in', str(dir_bam.absolute())])

    # Run the main script
    workflow = PACU(unparsed_args)
    workflow.run()
