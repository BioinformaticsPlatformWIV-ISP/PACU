from setuptools import setup, find_packages

from pacu.version import __version__


setup(
    name='PACU',
    version=__version__,
    description='Workflow for whole genome sequencing based phylogeny of Illumina and ONT data.',
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
        # 'PyVCF==0.6.8',
        # 'PyYAML==6.0',
        # 'UpSetPlot==0.8.0',
        # 'beautifulsoup4==4.11.1',
        # 'biopython==1.79',
        # 'matplotlib==3.7.1',
        # 'pandas==2.1.0',
        # 'pyvcf3==1.0.3',
        # 'snakemake==7.18.2',
        # 'yattag==1.14.0'
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
