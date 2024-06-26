# PACU
PACU is a workflow for whole genome sequencing based phylogeny of Illumina and ONT R9/R10 data.

PACU stands for the Prokaryotic Awesome variant Calling Utility and is named after an omnivorous fish (that eats both 
Illumina and ONT reads).

### PACU is also available on our public [Galaxy instance](https://galaxy.sciensano.be/) (registration required).

----

## INSTALLATION

```
pip install pacu_snp
```

The PACU dependencies can be installed via Conda or manually.

### CONDA installation

**Tip**: It is recommended to use `mamba` to install the dependencies, but `conda` can also be used.

```
cd PACU/
mamba env create -f environment.yml;
conda activate pacu_env;
pip install pacu_snp;
```

**Note:** `MEGA` is currently not available through Conda, it can be installed manually from the link below, or 
`IQ-Tree` can be used instead.

### Manual installation

The PACU workflow has the following dependencies:
- [BEDTools 2.27.1](https://github.com/arq5x/bedtools2/releases/tag/v2.27.1)
- [bcftools 1.17](https://github.com/samtools/bcftools/releases/tag/1.17)
- [FigTree 1.4.4](http://tree.bio.ed.ac.uk/software/figtree/)
- [samtools 1.17](https://github.com/samtools/samtools/releases/tag/1.17)
- [snpdists 0.8.2](https://github.com/tseemann/snp-dists)
- [MEGA 10.0.4](https://www.megasoftware.net/)
- [IQ-Tree 2.2.5](https://github.com/iqtree/iqtree2)

The mapping script has the following additional dependencies:
- [Trimmomatic 0.39](https://github.com/usadellab/Trimmomatic)
- [SeqKit 2.3.1](https://github.com/shenwei356/seqkit)
- [Bowtie2 2.5.1](https://github.com/BenLangmead/bowtie2)
- [Minimap2 2.26](https://github.com/lh3/minimap2)

The corresponding binaries should be in your PATH to run the workflow. 
Other versions of these tools may work, but have not been tested.

The required Python packages are listed in the `requirements.txt` file. 
Python 3.9 or 3.10 is recommended for a manual installation.

```
virtualenv pacu_env --python=python3.10;
. pacu_env/bin/activate;
pip install pacu_snp;
```

Note: Make sure you are in the directory containing the `setup.py` script when running the `pip install` command.

## USAGE

```
usage: PACU [-h] [--ilmn-in ILMN_IN] [--ont-in ONT_IN] --ref-fasta REF_FASTA [--ref-bed REF_BED] [--dir-working DIR_WORKING] --output OUTPUT
            [--output-html OUTPUT_HTML] [--use-mega] [--include-ref] [--min-snp-af MIN_SNP_AF] [--min-snp-qual MIN_SNP_QUAL]
            [--min-snp-depth MIN_SNP_DEPTH] [--min-snp-dist MIN_SNP_DIST] [--min-global-depth MIN_GLOBAL_DEPTH] [--min-mq-depth MIN_MQ_DEPTH]
            [--bcftools-filt-af1] [--image-width IMAGE_WIDTH] [--image-height IMAGE_HEIGHT] [--threads THREADS] [--version]

options:
  -h, --help            show this help message and exit
  --ilmn-in ILMN_IN     Directory with Illumina input BAM files
  --ont-in ONT_IN       Directory with ONT input BAM files
  --ref-fasta REF_FASTA
                        Reference FASTA file
  --ref-bed REF_BED     BED file with phage regions
  --dir-working DIR_WORKING
                        Working directory
  --output OUTPUT       Output directory
  --output-html OUTPUT_HTML
                        Output report name
  --use-mega            If set, MEGA is used for the construction of the phylogeny (instead of IQ-TREE)
  --include-ref         If set, the reference genome is included in the phylogeny
  --min-snp-af MIN_SNP_AF
                        Minimum allele frequency for variants
  --min-snp-qual MIN_SNP_QUAL
                        Minimum SNP quality
  --min-snp-depth MIN_SNP_DEPTH
                        Minimum SNP depth
  --min-snp-dist MIN_SNP_DIST
                        Minimum distance between SNPs
  --min-global-depth MIN_GLOBAL_DEPTH
                        Minimum depth for all samples to include positions in SNP analysis
  --min-mq-depth MIN_MQ_DEPTH
                        MQ cutoff for samtools depth
  --bcftools-filt-af1   If enabled, allele frequency filtering also considers the VAF value
  --image-width IMAGE_WIDTH
                        Image width
  --image-height IMAGE_HEIGHT
                        Image height
  --threads THREADS
  --version             Print version and exit
```

**Note:** The location of the temporary directory can be changed by setting the `TMPDIR` environment variable.

### Basic usage example

The PACU workflow requires BAM files as input with reads mapped to a reference genome. 
Illumina data can be provided using the `--ilmn-in` option, ONT data can be provided using the `--ont-in` option.

```
PACU \
    --ilmn-in in/ilmn/ \
    --ont-in in/ont/ \
    --ref-fasta ref.fasta \
    --output output/ \
    --dir-working work/ \
    --threads 8
```

### Read mapping

A script is included to map reads to a reference genome in FASTA format for both ONT and Illumina data.
The resulting BAM files can be used as input for the SNP workflow. The `--trim` option can be used to perform read
trimming before mapping.

*Illumina data*
```
PACU_map \
    --ref-fasta genome.fasta \
    --data-type illumina \
    --fastq-illumina reads_1.fastq.gz reads_2.fastq.gz \
    --output mapped.bam \
    --threads 4
```
*ONT data*
```
PACU_map \
    --ref-fasta genome.fasta \
    --data-type ont \
    --fastq-ont reads_ont.fastq.gz \
    --output mapped.bam \
    --threads 4
```

## TESTING

A test dataset is available under `resources/testdata/bam`, these files contain *Escherichia coli* reads mapped to a 
small part of the *E. coli* NC_002695.2 genome. This is a not a real dataset, and should only be used for testing.

The complete workflow can be tested using the following command:
```
pytest --log-cli-level=DEBUG pacu/tests/test_workflow.py
```

## CONTACT
[Create an issue](https://github.com/BioinformaticsPlatformWIV-ISP/PACU/issues) to report bugs, propose new functions or ask for help.

## CITATION
If you use this tool, please consider citing our [publication](https://pubmed.ncbi.nlm.nih.gov/38441926/).

-----

Copyright - 2024 Bert Bogaerts <bert.bogaerts@sciensano.be>
