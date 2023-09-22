from pathlib import Path

##############
# Input data #
##############
tech_by_key = {}
# Illumina data
for path_fq_fwd in Path('ilmn').glob('*_1P.fastq.gz'):
    tech_by_key[path_fq_fwd.name.replace('_1P.fastq.gz', '')] = 'ilmn'

# ONT data
for path_fq in Path('ont').glob('*.fastq.gz'):
    tech_by_key[path_fq.name.replace('.fastq.gz', '')] = 'ont'

ref_name = Path(config['ref']).name.replace('.fasta', '')

#########
# Rules #
#########
rule all:
    input:
        BAM = [f'mapping/{key}.bam' for key in tech_by_key.keys()]

rule prep_reference:
    input:
        FASTA = config['ref']
    output:
        FASTA = f'ref/{ref_name}.fasta',
        MNI = f'ref/{ref_name}.mni'
    shell:
        """
        cp {input.FASTA} {output.FASTA};
        bowtie2-build {output.FASTA} {output.FASTA};
        minimap2 -x map-ont -d {output.MNI} {output.FASTA};
        """

rule map_illumina:
    """
    Maps the Illumina reads.
    """
    input:
        FQ_fwd = 'ilmn/{key}_1P.fastq.gz',
        FQ_rev = 'ilmn/{key}_2P.fastq.gz',
        FASTA = rules.prep_reference.output.FASTA
    output:
        SAM = temporary('mapping/ilmn/{key}.sam')
    threads: 8
    shell:
        """
        bowtie2 -x {input.FASTA} --sensitive -1 {input.FQ_fwd} -2 {input.FQ_rev} -S {output.SAM} -p {threads};
        """

rule map_minimap2:
    """
    Maps the ONT reads using Minimap2.
    """
    input:
        FQ = 'ont/{key}.fastq.gz',
        MNI = rules.prep_reference.output.MNI
    output:
        SAM = temporary('mapping/ilmn/{key}.sam')
    threads: 8
    shell:
        """
        minimap2 -a {input.MNI} {input.FQ} -t {threads} > {output.SAM};
        """

rule sam_to_bam:
    """
    Converts the SAM file to sorted BAM format.
    """
    input:
        SAM = lambda wildcards: rules.map_illumina.output.SAM if tech_by_key[wildcards.key] == 'illumina' else rules.map_minimap2.output.SAM
    output:
        BAM = 'mapping/{key}.bam'
    shell:
        """
        samtools view -b -F 4 {input.SAM} | samtools sort - > {output.BAM};
        samtools index {output.BAM};
        """
