#!/usr/bin/env python
import argparse
from pathlib import Path
from typing import Tuple, List

from sciensnp import initialize_logging, logger, ScienSNP
from sciensnp.app.utils import workflowutils


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


if __name__ == '__main__':
    initialize_logging()
    logger.info(f'Running ScienSNP through Galaxy')
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
        for path, name in args.bam_ilmn:
            try:
                (dir_bam / workflowutils.sanitize_bam_input(name)).symlink_to(path)
            except FileExistsError:
                logger.debug(f"Symlink for '{name}' already exists")
        unparsed_args.extend(['--ilmn-in', str(dir_bam.absolute())])

    # Symlink the BAM files (ONT)
    if args.bam_ont is not None:
        logger.info(f'Creating symlinks for ont input ({len(args.bam_ont)} datasets)')
        dir_bam = dir_in / 'bam_ont'
        dir_bam.mkdir(exist_ok=True)
        for path, name in args.bam_ont:
            try:
                (dir_bam / workflowutils.sanitize_bam_input(name)).symlink_to(path)
            except FileExistsError:
                logger.debug(f"Symlink for '{name}' already exists")
        unparsed_args.extend(['--ont-in', str(dir_bam.absolute())])

    # Run the main script
    workflow = ScienSNP(unparsed_args)
    workflow.run()
