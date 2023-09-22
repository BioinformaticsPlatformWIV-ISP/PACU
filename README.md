# ScienSNP
Workflow for whole genome sequencing based phylogeny of Illumina and ONT data.

### ScienSNP is also available on our public [Galaxy instance](https://galaxy.sciensano.be/).

----

## INSTALLATION

```
git clone https://bioit-git.sciensano.be/bioit-shared/bb_sciensnp.git
```

The ScienSNP dependencies can be installed via Conda or manually.

### CONDA installation

**Tip**: It is recommended to use `mamba` to install the dependencies, but `conda` can also be used.

```
mamba env create -f environment.yml;
conda activate sciensnp_env;
pip install .;
```

### Manual installation

The ScienSNP workflow has the following dependencies:
- [BEDTools 2.27.1](https://github.com/arq5x/bedtools2/releases/tag/v2.27.1)
- [bcftools 1.17](https://github.com/samtools/bcftools/releases/tag/1.17)
- [FigTree 1.4.4](http://tree.bio.ed.ac.uk/software/figtree/)
- [MEGA 10.0.4](https://www.megasoftware.net/)
- [samtools 1.17](https://github.com/samtools/samtools/releases/tag/1.17)
- [snpdists 0.8.2](https://github.com/tseemann/snp-dists)

The corresponding binaries should be in your PATH to run the workflow. 
Other versions of these tools may work, but have not been tested.

The required Python packages are listed in the `requirements.txt` file. 
Python 3.9 is recommended for a manual installation.

```
virtualenv sciensnp_env --python=python3.9
. sciensnp_env/bin/activate;
pip install -r requirements.txt 
pip install .;
```

## USAGE

```
usage: run_sciensnp.py [-h] [--ilmn-in ILMN_IN] [--ont-in ONT_IN] --ref-fasta REF_FASTA [--ref-bed REF_BED]
                       [--dir-working DIR_WORKING] [--output OUTPUT] [--calling-method {clair3,samtools}] [--include-ref]
                       [--min-snp-af MIN_SNP_AF] [--min-snp-qual MIN_SNP_QUAL] [--min-snp-depth MIN_SNP_DEPTH]
                       [--min-snp-dist MIN_SNP_DIST] [--min-global-depth MIN_GLOBAL_DEPTH] [--min-mq-depth MIN_MQ_DEPTH]
                       [--bcftools-filt-af1] [--image-width IMAGE_WIDTH] [--image-height IMAGE_HEIGHT] [--threads THREADS]
                       [--version]

optional arguments:
  -h, --help            show this help message and exit
  --ilmn-in ILMN_IN     Directory with Illumina input BAM files
  --ont-in ONT_IN       Directory with ONT input BAM files
  --ref-fasta REF_FASTA
                        Reference FASTA file
  --ref-bed REF_BED     BED file with phage regions
  --dir-working DIR_WORKING
                        Working directory
  --output OUTPUT       Output directory
  --calling-method {clair3,samtools}
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

The ScienSNP workflow requires BAM files as input with reads mapped to a reference genome. 
Illumina data can be provided using the `--ilmn-in` option, ONT data can be provided using the `--ont-in` option.

```
run_sciensnp.py \
    --ilmn-in in/ilmn/ \
    --ont-in in/ont/ \
    --ref-fasta ref.fasta \
    --output output/ \
    --dir-working work/ \
    --threads 8
```

## TESTING

A test dataset is available under `resources/testdata/bam`, these files contain *Escherichia coli* reads mapped to a 
small part of the *E. coli* NC_002695.2 genome. This is a not a real dataset, and should only be used for testing.

The complete workflow can be tested using the following command:
```
pytest --log-cli-level=DEBUG sciensnp/tests/test_workflow.py
```

## CITATION
If you use this tool, please consider citing our [TODO](https://example.com).

-----

Copyright - 2023 Bert Bogaerts <bert.bogaerts@sciensano.be>
