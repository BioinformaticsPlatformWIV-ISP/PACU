from setuptools import setup, find_packages

from sciensnp.version import __version__

import pprint
pprint.pprint(find_packages())

setup(
    name='ScienSNP',
    version=__version__,
    description='Workflow for whole genome sequencing based phylogeny of Illumina and ONT data.',
    url="https://github.com/TODO",
    author="Bert Bogaerts",
    author_email="bert.bogaerts@sciensano.be",
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    keywords="nanopore sequencing plotting quality control",
    packages=find_packages() + ['sciensnp.resources'],
    python_requires=">=3",
    install_requires=[
        # 'PyVCF==0.6.8',
        # 'PyYAML==6.0',
        # 'beautifulsoup4==4.11.1',
        # 'biopython==1.79',
        # 'matplotlib==3.7.1',
        # 'pandas==2.1.0',
        # 'snakemake==7.18.2',
        # 'yattag==1.14.0',
    ],
    package_data={
        'sciensnp': [
            'resources/citations/*.json',
            'resources/figtree_template.txt',
            'resources/mega/*.mao',
            'resources/snp_workflow.smk'
            'resources/style.css',
            'resources/testdata/bam/ilmn/*.bai',
            'resources/testdata/bam/ilmn/*.bam',
            'resources/testdata/bam/ont/*.bai',
            'resources/testdata/bam/ont/*.bam',
            'resources/testdata/figtree_template.txt',
            'resources/testdata/mega/model_selection_output.tsv',
            'resources/testdata/phylogeny.nwk',
            'resources/testdata/snp_matrix.fasta',
            'resources/testdata/vcfs/*.vcf',
        ]},
    entry_points={
        'console_scripts': [
            'ScienSNP=sciensnp.run_sciensnp:main',
        ],
    }
)
