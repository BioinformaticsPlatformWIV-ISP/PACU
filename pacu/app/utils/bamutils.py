import re
import tempfile
from pathlib import Path

from pacu import Command


def add_custom_tag(name: str, value: str, bam_in: Path, bam_out: Path) -> None:
    """
    Adds the custom tag to the input BAM file.
    :param name: Tag name
    :param value: Tag value
    :param bam_in: Input BAM file
    :param bam_out: Output BAM file
    :return: None
    """
    with tempfile.TemporaryDirectory(prefix='pacu') as dir_temp:
        # Extract the original header
        command = Command(f'samtools view -H {bam_in}')
        command.run(Path(dir_temp))
        if not command.exit_code == 0:
            raise RuntimeError(f'Error extracting header from BAM file: {bam_in}')
        header = command.stdout

        # Add custom tag and save header
        header += f'@CO\t{name}:{value}\n'
        path_header_updated = Path(dir_temp, f'header_updated.txt')
        with path_header_updated.open('w') as handle:
            handle.write(header)

        # Create output BAM file with custom tag
        bam_out.parent.mkdir(parents=True, exist_ok=True)
        command = Command(f'samtools reheader {path_header_updated} {bam_in.absolute()} > {bam_out.absolute()}')
        command.run(Path(dir_temp))
        if not command.exit_code == 0:
            raise RuntimeError(f'Error replacing BAM header: {command.stderr}')


def read_custom_tag(bam_in: Path, name: str) -> str:
    """
    Reads a custom tag from the input BAM file.
    :param bam_in: Input BAM file
    :param name: Tag name
    :return: Tag value
    """
    command = Command(f'samtools view -H {bam_in}')
    command.run(bam_in.parent)
    if not command.exit_code == 0:
        raise RuntimeError(f'Cannot extract header from: {bam_in}')
    for line in command.stdout.splitlines():
        if not line.startswith('@CO'):
            continue
        m = re.match(r'@CO\t(.*):(.*)', line)
        if (m is None) or (m.group(1) != name):
            continue
        return m.group(2)
    raise ValueError(f"Cannot extract tag '{name}'")
