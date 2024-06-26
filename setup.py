from pathlib import Path

from setuptools import setup, find_packages

from pacu.version import __version__

with open(Path(__file__).parent / 'README.md', encoding='utf-8') as handle:
    long_description = handle.read()


setup(
    name='pacu_snp',
    version=__version__,
    description='Workflow for whole genome sequencing based phylogeny of Illumina and ONT data.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/BioinformaticsPlatformWIV-ISP/PACU/',
    author='Bert Bogaerts',
    author_email='bioit@sciensano.be',
    license='GPL',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3",
    ],
    keywords='nanopore illumina sequencing phylogeny',
    packages=find_packages() + ['pacu.resources'],
    python_requires='>=3',
    install_requires=[
        'PuLP==2.7.0',
        'PyYAML>=6.0',
        'UpSetPlot>=0.8.0',
        'beautifulsoup4>=4.11.1',
        'biopython',
        'matplotlib>=3.7.1',
        'numpy>=1.26.4',
        'pandas>=2.1.0',
        'pytest>=8.2.2',
        'pyvcf3>=1.0.3',
        'snakemake==7.18.2',
        'yattag>=1.14.0'
    ],
    package_data={
        'pacu': [
            'resources/citations/*.json',
            'resources/figtree_template.txt',
            'resources/mega/*.mao',
            'resources/snp_workflow.smk',
            'resources/style.css',
            'resources/testdata/bam/ilmn/*.bai',
            'resources/testdata/bam/ilmn/*.bam',
            'resources/testdata/bam/ont/*.bai',
            'resources/testdata/bam/ont/*.bam',
            'resources/testdata/fastq/*.fastq.gz',
            'resources/testdata/mega/model_selection_output.csv',
            'resources/testdata/vcfs/*.vcf',
            'resources/testdata/NC_002695.2-subset.fasta',
            'resources/testdata/phylogeny.nwk',
            'resources/testdata/snp_matrix.fasta',
        ]},
    entry_points={
        'console_scripts': [
            'PACU=pacu.run_pacu:main',
        ],
    }
)
