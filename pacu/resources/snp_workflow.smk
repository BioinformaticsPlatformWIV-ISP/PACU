import logging
import shutil
from importlib.resources import files
from pathlib import Path
import json
import pandas as pd
import vcf

from pacu.app.command import Command
from pacu.app.report.htmlreport import HtmlReport
from pacu.app.utils import workflowutils
from pacu.app.utils.loggingutils import initialize_logging

initialize_logging()

rule samtools_depth:
    """
    Determines the median depth of an input BAM file.
    """
    input:
        BAM = lambda wildcards: str(config['input'][wildcards.key]['bam'])
    output:
        TSV = 'depth/{key}/depth_{key}.tsv'
    params:
        min_mq = config['depth']['min_mq']
    threads: 1
    shell:
        """
        samtools depth {input.BAM} -a --min-MQ {params.min_mq} > {output.TSV}
        """

rule samtools_depth_parse:
    """
    Parses the samtools depth output.
    """
    input:
        TSV = rules.samtools_depth.output.TSV
    output:
        JSON = 'depth/{key}/depth_{key}.json'
    params:
        key = lambda wildcards: wildcards.key
    run:
        data_depth = pd.read_table(input.TSV,  names=['chr', 'pos', 'depth'], dtype={'chr': str})
        with open(output.JSON, 'w') as handle:
            json.dump({
                'key': params.key,
                'positions': len(data_depth),
                'median_depth': data_depth['depth'].median()
            }, handle, indent=2)

rule samtools_depth_combine:
    """
    Combines the depth across all datasets.
    """
    input:
        TSV = lambda wildcards: expand(rules.samtools_depth.output.TSV, key=config['input'].keys())
    output:
        TSV = 'combined/depth_all.tsv'
    params:
        names = config['input'].keys()
    run:
        records_out = []
        # noinspection PyTypeChecker
        for name, tsv_file in zip(params.names, [Path(x) for x in input.TSV]):
            data_depth_isolate = pd.read_table(tsv_file, names=['chr', 'pos', 'depth'])
            records_out.append({
                'key': name,
                'median_depth': data_depth_isolate['depth'].median(),
                'perc_covered': 100 * (sum(data_depth_isolate['depth'] > 0) / len(data_depth_isolate))
            })
        data_out = pd.DataFrame(records_out)
        data_out.to_csv(output.TSV, sep='\t', index=False)

rule variant_calling_faidx_ref:
    """
    Creates a FASTA index for the reference genome.
    """
    input:
        FASTA = str(config['reference']['fasta'])
    output:
        FASTA = f"ref/{Path(config['reference']['fasta']).name}"
    shell:
        """
        cp {input.FASTA} {output.FASTA}
        samtools faidx {output.FASTA}
        """

rule variant_calling_bcftools_mpileup:
    """
    Creates a pileup using bcftools based on a BAM input file.
    """
    input:
        BAM = lambda wildcards: str(config['input'][wildcards.key]['bam']),
        FASTA = rules.variant_calling_faidx_ref.output.FASTA
    output:
        VCF_GZ = 'variant_calling/{key}/pileup_{key}.vcf.gz'
    params:
        profile = lambda wildcards: 'illumina' if config['input'][wildcards.key]['tech'] == 'ilmn' else 'ont'
    shell:
        """
        bcftools mpileup {input.BAM} --fasta-ref {input.FASTA} --output {output.VCF_GZ} --output-type v \
            --config {params.profile}
        """

rule variant_calling_bcftools_call:
    """
    Creates a pileup using bcftools based on a BAM input file.
    """
    input:
        VCF_GZ = rules.variant_calling_bcftools_mpileup.output.VCF_GZ
    output:
        VCF_GZ = 'variant_calling/{key}/unfiltered_{key}.vcf.gz'
    params:
        mutation_rate = lambda wildcards: '1.1e-3' if config['input'][wildcards.key]['tech'] == 'ilmn' else '0.01'
    shell:
        """
        bcftools call {input.VCF_GZ} --output {output.VCF_GZ} --output-type v -P {params.mutation_rate} \
            --variants-only --skip-variants indels --ploidy 1 --consensus-caller
        """

rule variant_filtering_allele_freq:
    """
    Soft filters variants based on allele frequency.
    """
    input:
        VCF_GZ = rules.variant_calling_bcftools_call.output.VCF_GZ
    output:
        VCF_GZ = 'variant_filtering/{key}/filtered_af-{key}.vcf.gz'
    params:
        expr = f"((DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])) < {config['filters']['min_af']}"
    shell:
        """
        bcftools filter --output-type z --soft-filter af --exclude "{params.expr}" {input.VCF_GZ} \
            --output {output.VCF_GZ}
        """

rule variant_filtering_depth:
    """
    Soft filters variants based on depth.
    """
    input:
        VCF_GZ = rules.variant_filtering_allele_freq.output.VCF_GZ
    output:
        VCF_GZ = 'variant_filtering/{key}/filtered_depth-{key}.vcf.gz'
    params:
        min_dp = config['filters']['min_depth']
    shell:
        """
        bcftools filter --output-type z --soft-filter dp --exclude "DP<{params.min_dp}" {input.VCF_GZ} \
            --output {output.VCF_GZ}
        """

rule variant_filtering_qual:
    """
    Soft filters variants based on SNP quality.
    """
    input:
        VCF_GZ = rules.variant_filtering_depth.output.VCF_GZ
    output:
        VCF_GZ = 'variant_filtering/{key}/filtered_qual-{key}.vcf.gz'
    params:
        min_qual = config['filters']['min_qual']
    shell:
        """
        bcftools filter --output-type z --soft-filter qual --exclude "QUAL<{params.min_qual}" {input.VCF_GZ} \
            --output {output.VCF_GZ}
        """

rule gubbins_consensus_create:
    """
    Creates the consensus sequence.
    """
    input:
        VCF_GZ = rules.variant_filtering_qual.output.VCF_GZ,
        FASTA = str(config['reference']['fasta'])
    output:
        FASTA = temporary('gubbins/{key}/consensus_{key}.fasta.tmp')
    params:
        VCF_GZ = lambda wildcards: f'gubbins/{wildcards.key}/tmp_variants_{wildcards.key}.vcf.gz'
    shell:
        """
        # Compress and index input VCF file
        cp {input.VCF_GZ} {params.VCF_GZ}
        bcftools index -f {params.VCF_GZ}

        # Create consensus sequence
        bcftools consensus --exclude 'FILTER!="PASS"' --fasta-ref {input.FASTA} {params.VCF_GZ} > {output.FASTA}

        # Remove temporary files
        rm {params.VCF_GZ}*
        """

rule gubbins_consensus_update_header:
    """
    Renames the header of the consensus FASTA file. 
    """
    input:
        FASTA = rules.gubbins_consensus_create.output.FASTA
    output:
        FASTA = 'gubbins/{key}/consensus_{key}.fasta'
    params:
        key = lambda wildcards: wildcards.key
    run:
        from Bio import SeqIO
        with open(input.FASTA) as handle:
            seqs = list(SeqIO.parse(handle, 'fasta'))
        if len(seqs) > 1:
            raise ValueError("More than one sequence found")
        seqs[0].description = ''
        seqs[0].id = params.key
        with open(output.FASTA, 'w') as handle:
            SeqIO.write(seqs, handle, 'fasta')

rule gubbins_merge_alignment:
    """
    Merges the alignment.
    """
    input:
        FASTA = expand(rules.gubbins_consensus_update_header.output.FASTA, key=config['input'].keys())
    output:
        FASTA = 'gubbins/wga_all.fasta'
    run:
        from Bio import SeqIO
        seqs_out = []
        for fasta_file in [Path(x) for x in input.FASTA]:
            with fasta_file.open() as handle:
                seq_in = next(SeqIO.parse(handle, 'fasta'))
                # Check for duplicates
                if any(str(s.seq) == seq_in.seq for s in seqs_out):
                    continue
                seqs_out.append(seq_in)

        # Save output file
        with open(output.FASTA, 'w') as handle:
            SeqIO.write(seqs_out, handle, 'fasta')

rule gubbins_run:
    """
    Runs gubbins on the SNP matrix.
    """
    input:
        FASTA = rules.gubbins_merge_alignment.output.FASTA
    output:
        GFF = 'gubbins/wga_all.recombination_predictions.gff'
    threads: 8
    run:
        from Bio import SeqIO

        # Check if all sequences are the same
        with open(input.FASTA) as handle:
            all_seqs = list(SeqIO.parse(handle, 'fasta'))
            logging.info(f'Running Gubbins on {len(all_seqs)} sequences')

        # Run Gubbins
        command = Command(' '.join(['run_gubbins.py', str(Path(input.FASTA).absolute()), f'--threads {threads}']))
        command.run(Path(output.GFF).parent)
        if not command.exit_code == 0:
            if 'Failed while building the tree' in command.stderr:
                logging.warning("Gubbins tree building failed, creating empty output file")
                Path(output.GFF).touch()
            elif 'Three or more sequences are required' in command.stderr:
                logging.warning("Gubbins tree building failed, creating empty output file")
                Path(output.GFF).touch()
            else:
                raise RuntimeError(f"Unknown Error executing gubbins: {command.stderr}")

rule gubbins_to_bed:
    """
    Converts the gubbins output to BED format.
    """
    input:
        FASTA = str(config['reference']['fasta']),
        GFF = rules.gubbins_run.output.GFF
    output:
        BED = 'gubbins/gubbins.bed'
    run:
        from Bio import SeqIO

        # Get sequence name
        with open(input.FASTA) as handle:
            seq_name = SeqIO.read(handle, 'fasta').id

        # Construct the command
        command = Command(' '.join([
            f'bedtools sort -i {Path(input.GFF).absolute()} |',
            'bedtools merge -i - |',
            f"sed 's/SEQUENCE/{seq_name}/'",
            f'> {Path(output.BED).absolute()}'
        ]))
        command.run(Path(output.BED).parent)
        if not command.exit_code == 0:
            raise RuntimeError(f"Error creating BED file: {command.stderr}")

rule gubbins_empty_bed:
    """
    Creates an empty BED file when gubbins is disabled. 
    """
    output:
        BED = 'gubbins/empty.bed'
    shell:
        """
        touch {output.BED}
        """

rule gubbins_select_bed:
    """
    Rule that selects the BED file for gubbins dependent on the pipeline setting.
    """
    input:
        BED_gubbins = rules.gubbins_to_bed.output.BED if config.get('skip_gubbins', False) is False else [],
        BED_empty = rules.gubbins_empty_bed.output.BED if config.get('skip_gubbins', False) is True else []
    output:
        BED = 'gubbins/selected.bed'
    params:
        skip_gubbins = config.get('skip_gubbins', False)
    run:
        if len(input.BED_gubbins) > 0:
            shutil.copyfile(input.BED_gubbins, output.BED)
        else:
            shutil.copyfile(input.BED_empty, output.BED)

rule region_filtering_collect_low_depth_regions:
    """
    Extracts the low coverage regions
    """
    input:
        FASTA = str(config['reference']['fasta']),
        TSV = lambda wildcards: expand(rules.samtools_depth.output.TSV, key=config['input'].keys())
    output:
        BED = 'region_filtering/low_depth.bed'
    params:
        min_depth = config['depth']['min_depth']
    run:
        # Parse input data
        data_depth_all = pd.DataFrame()
        # noinspection PyTypeChecker
        for tsv_file in [Path(x) for x in input.TSV]:
            data_depth_isolate = pd.read_table(tsv_file, names=['chr', 'pos', 'depth'], dtype={'chr': str})
            logging.info(f"Parsing: {tsv_file}")
            data_depth_isolate['isolate'] = tsv_file.parent.name
            data_depth_all = pd.concat([data_depth_all, data_depth_isolate])

        # Determine the minimum depth at each position across all samples
        data_min_depth = data_depth_all.groupby(by=['chr', 'pos'])['depth'].min()
        df = pd.DataFrame([{
            'chr': chr,
            'pos': pos,
            'min_depth': min_depth} for (chr, pos), min_depth in data_min_depth.items()
        ])
        # noinspection PyTypeChecker
        logging.info(f"Positions passed depth filtering: {sum(df['min_depth'] > params.min_depth):,}")

        # Create data frame containing only the low depth positions
        df_low_depth = df[df['min_depth'] < params.min_depth].copy()

        # Determine the start of each region
        df_low_depth['shift_pos'] = df_low_depth['pos'].shift()
        df_low_depth['shift_chr'] = df_low_depth['chr'].shift()
        df_low_depth['is_new'] = df_low_depth.apply(lambda x: workflowutils.is_new_region(x), axis=1)

        # Create column with region number using cumsum
        df_low_depth['cumsum'] = df_low_depth['is_new'].cumsum()

        # Extract the start and stop position of each region
        df_regions = df_low_depth.groupby([df_low_depth['cumsum'], df_low_depth['chr']]).agg({'pos': ['min', 'max']})

        # Create output file
        with Path(output.BED).absolute().open('w') as handle:
            for (_, chr), row in df_regions.iterrows():
                # Minus one on the start coordinate because samtools is 1 based and BED 0 based
                handle.write('\t'.join([chr, str(row['pos']['min'] - 1), str(row['pos']['max'])]))
                handle.write('\n')

rule region_filtering_merge_bed_files:
    """
    Merges the BED files with phage regions and recombinant regions detected by Gubbins.
    """
    input:
        BED_a = rules.gubbins_select_bed.output.BED,
        BED_b = config['reference']['bed_phages'],
        BED_c = rules.region_filtering_collect_low_depth_regions.output.BED
    output:
        BED = 'region_filtering/merged.bed'
    run:
        paths_bed = input.BED_a, input.BED_b, input.BED_c
        if all(workflowutils.count_regions(p) == 0 for p in paths_bed):
            logging.warning("No regions to merge, creating empty output file")
            Path(output.BED).touch()
        else:
            command = Command(' '.join([
                'bedtools', 'multiinter', '-i', *[str(Path(i).absolute()) for i in paths_bed],
                f'> {Path(output.BED).absolute()}'
            ]))
            command.run(Path(output.BED).parent.absolute())
            if not command.exit_code == 0:
                raise RuntimeError(f"Error merging BED files: {command.stderr}")

rule region_filtering_plot:
    """
    Plots the overlap between the BED files used for filtering variants.
    """
    input:
        FASTA = config['reference']['fasta'],
        BED_phages = config['reference']['bed_phages'],
        BED_gubbins = rules.gubbins_select_bed.output.BED,
        BED_depth = rules.region_filtering_collect_low_depth_regions.output.BED
    output:
        PNG = 'stats/region_filtering.png',
        TSV = 'stats/filtered_positions.tsv'
    run:
        from Bio import SeqIO
        from pacu.app.utils import reportutils

        # Get the size of the reference genome
        with open(input.FASTA) as handle:
            ref_size = sum([len(s) for s in SeqIO.parse(handle, 'fasta')])

        # Calculate data
        data_overlaps = workflowutils.calculate_overlaps(
            ref_size, Path(input.BED_phages), Path(input.BED_gubbins), Path(input.BED_depth))
        data_overlaps.to_csv(output.TSV, sep='\t', index=False)

        # Create the plot
        reportutils.create_upsetplot_overlap(data_overlaps, Path(output.PNG))

rule region_filtering_combine_stats:
    """
    Collect statistics on the input BED files.
    """
    input:
        FASTA = config['reference']['fasta'],
        BED_phages = config['reference']['bed_phages'],
        BED_gubbins = rules.gubbins_select_bed.output.BED,
        BED_depth = rules.region_filtering_collect_low_depth_regions.output.BED,
        BED_merged = rules.region_filtering_merge_bed_files.output.BED
    output:
        TSV = 'stats/stats_regions.tsv',
        TSV_overlap = 'stats/stats_region_overlap.tsv'
    run:
        from Bio import SeqIO
        import itertools

        # Determine the reference genome size
        with open(input.FASTA) as handle:
            ref_genome_size = sum([len(s) for s in SeqIO.parse(handle, 'fasta')])

        # Determine the reference genome size
        with open(input.FASTA) as handle:
            reference_genome_size = sum([len(s) for s in SeqIO.parse(handle, 'fasta')])

        # Parse BED statistics
        records_out = []
        nb_pos_by_bed_file = {}
        # noinspection PyUnresolvedReferences
        for key, bed_file in input.items():
            if not key.startswith('BED'):
                continue
            positions_covered = workflowutils.count_covered_positions(Path(bed_file).absolute())
            records_out.append({
                'key': key,
                'nb_covered_positions': positions_covered,
                'perc_ref_covered': 100 * positions_covered / reference_genome_size
            })
            nb_pos_by_bed_file[str(bed_file)] = positions_covered

        # Determine overlaps
        records_overlap = []
        for bed_file_a, bed_file_b in itertools.combinations([Path(x) for x in (
                input.BED_phages, input.BED_depth, input.BED_merged, input.BED_gubbins)], r=2):
            records_overlap.append({
                'bed_file_a': bed_file_a.stem.split('-')[0],
                'size_a': nb_pos_by_bed_file[str(bed_file_a)],
                'bed_file_b': bed_file_b.stem.split('-')[0],
                'size_b': nb_pos_by_bed_file[str(bed_file_b)],
                'overlap': workflowutils.count_overlap(bed_file_a.absolute(), bed_file_b.absolute())
            })
        pd.DataFrame(records_overlap).to_csv(output.TSV_overlap, sep='\t', index=False)

        # Create output file
        pd.DataFrame(records_out).to_csv(output.TSV, sep='\t', index=False)

rule region_filtering_remove_snps:
    """
    Removes SNPs located in problematic regions.
    """
    input:
        VCF_GZ = rules.variant_filtering_qual.output.VCF_GZ,
        BED = rules.region_filtering_merge_bed_files.output.BED
    output:
        VCF = 'region_filtering/{key}/snps_{key}.filtered_reg.vcf'
    params:
        tempfile = lambda wildcards: Path('region_filtering') / wildcards.key / f'tmp_{wildcards.key}.filtered_reg.vcf.gz'
    shell:
        """
        # Compress and index input file
        bcftools view --output-type z {input.VCF_GZ} > {params.tempfile}
        bcftools index -f {params.tempfile}

        # Filter VCF file
        if [ -s {input.BED} ]; then
            bcftools filter {params.tempfile} --targets-file ^{input.BED} > {output.VCF}
        else
            # If empty BED file, copy the original VCF
            bcftools view --output-type v {input.VCF_GZ} > {output.VCF}
        fi

        # Remove temporary file
        rm {params.tempfile}*
        """

rule variant_filtering_distance:
    """
    Filters SNPs based on distance to other SNPs.
    """
    input:
        VCF_GZ = rules.region_filtering_remove_snps.output.VCF
    output:
        VCF = 'region_filtering/{key}/snps_{key}.filtered_reg_dist.vcf'
    params:
        min_dist = config['filters']['min_dist']
    run:
        workflowutils.filter_snp_distance(Path(input.VCF_GZ), Path(output.VCF), params.min_dist)

rule collect_region_filtering_stats:
    """
    Collects the statistics for the region filtering.
    """
    input:
        VCF_before = lambda wildcards: expand(rules.variant_calling_bcftools_call.output.VCF_GZ, key=config['input'].keys()),
        VCF_after = lambda wildcards: expand(rules.variant_filtering_distance.output.VCF, key=config['input'].keys())
    output:
        TSV = 'stats/stats_region_filtering.tsv'
    run:
        import gzip

        records_out = []
        # noinspection PyTypeChecker
        for vcf_before, vcf_after in zip([Path(x) for x in input.VCF_before], [Path(x) for x in input.VCF_after]):
            with gzip.open(vcf_before) as handle:
                snps_before = list(vcf.Reader(handle))

            with vcf_after.open() as handle:
                snps_after = list(vcf.Reader(handle))

            records_out.append({
                'sample': vcf_before.parent.name,
                'nb_snps_in': len(snps_before),
                'nb_snps_out': len(snps_after),
                'perc_passed': 100 * len(snps_after) / max(len(snps_before), 1),
                'nb_snps_filt_in': sum(len(s.FILTER) == 0 for s in snps_before),
                'nb_snps_filt_out': sum(len(s.FILTER) == 0 for s in snps_after),
            })
        data_out = pd.DataFrame(records_out)
        data_out['perc_passed_filt'] = data_out.apply(
            lambda x: 100 * x['nb_snps_filt_out'] / max(x['nb_snps_filt_in'], 1), axis=1)
        data_out.to_csv(output.TSV, sep='\t', index=False)

rule collect_variant_calling_stats:
    """
    Collects the stats for the variant calling.
    """
    input:
        VCF = rules.variant_filtering_distance.output.VCF
    output:
        JSON = 'variant_filtering/{key}/stats_{key}.filtered.json'
    params:
        key = lambda wildcards: wildcards.key
    run:
        # noinspection PyTypeChecker
        with open(input.VCF) as handle:
            vcf_records = list(vcf.Reader(handle))
        with open(output.JSON, 'w') as handle:
            json.dump({
                'key': wildcards.key,
                'nb_variants': len(vcf_records),
                'nb_variants_filtered': sum(len(rec.FILTER) == 0 for rec in vcf_records)
            }, handle, indent=2)

rule combine_variant_calling_stats:
    """
    Combines the variant calling stats for all samples.
    """
    input:
        JSON = lambda wildcards: expand(rules.collect_variant_calling_stats.output.JSON, key=config['input'].keys()),
        JSON_depth = lambda wildcards: expand(rules.samtools_depth_parse.output.JSON, key=config['input'].keys())
    output:
        TSV = 'stats/stats_filtering.tsv'
    run:
        records_out = []
        # noinspection PyTypeChecker
        for json_stats, json_depth in zip([Path(x) for x in input.JSON], [Path(x) for x in input.JSON_depth]):
            with json_depth.open() as handle:
                median_depth = int(json.load(handle)['median_depth'])
            with json_stats.open() as handle:
                records_out.append({'depth': median_depth, **json.load(handle)})
        data_vc = pd.DataFrame(records_out)
        data_vc['perc_passed'] = 100 * data_vc['nb_variants_filtered'] / data_vc['nb_variants']
        data_vc.insert(0, 'key', data_vc.pop('key'))
        data_vc.to_csv(output.TSV, sep='\t', index=False)

rule create_snp_matrix:
    """
    Creates a SNP matrix from the filtered VCF files.
    """
    input:
        VCF = expand(rules.variant_filtering_distance.output.VCF, key=config['input'].keys())
    output:
        FASTA = 'tree/snp_matrix.fasta',
        TSV = 'combined/snp_positions.tsv'
    params:
        include_ref = False,
        names = list(config['input'].keys())
    run:
        from pacu.app.utils import snpmatrixutils
        snpmatrixutils.create_snp_matrix(
            paths_vcf=[Path(x) for x in input.VCF],
            names=params.names,
            path_out=Path(output.FASTA).absolute(), path_tsv=Path(output.TSV).absolute(),
            include_ref=params.include_ref)

rule mega_model_selection:
    """
    Performs model selection using MEGA.
    """
    input:
        FASTA = rules.create_snp_matrix.output.FASTA
    output:
        CSV = 'tree/mega/model_selection.csv',
        JSON = 'tree/mega/model_selection.json'
    params:
        branch_swap_filter = 'Very weak',
        missing_data_treatment= 'partial_deletion',
        site_cov_cutoff=50
    threads: 8
    run:
        from pacu.app.utils import megautils

        # Run model selection
        megautils.run_model_selection(
            path_fasta=Path(input.FASTA).absolute(),
            path_out=Path(output.CSV).absolute(),
            dir_=Path(output.CSV).absolute().parent,
            branch_swap_filter=params.branch_swap_filter,
            missing_data_treatment=params.missing_data_treatment,
            site_cov_cutoff=params.site_cov_cutoff,
            threads=threads
        )

        # Extract model name
        model_info = megautils.parse_model_selection_csv(Path(output.CSV).absolute())
        with open(output.JSON, 'w') as handle:
            json.dump({'model': model_info['model']}, handle, indent=2)

rule mega_construct_tree:
    """
    Constructs the phylogeny using MEGA.
    """
    input:
        FASTA = rules.create_snp_matrix.output.FASTA,
        CSV = rules.mega_model_selection.output.CSV
    output:
        NWK = 'tree/mega/phylogeny.nwk'
    params:
        branch_swap_filter = 'Very weak',
        missing_data_treatment = 'partial_deletion',
        site_cov_cutoff = 50,
        bootstrap_replicates = 100,
        heuristic_method = 'SPR3',
        initial_tree = 'NJ',
        gamma_categories = 5
    threads: 8
    run:
        from pacu.app.utils import megautils
        megautils.run_tree_building(
            path_fasta=Path(input.FASTA).absolute(),
            path_csv=Path(input.CSV).absolute(),
            path_out=Path(output.NWK).absolute(),
            dir_=Path(output.NWK).absolute().parent,
            branch_swap_filter=params.branch_swap_filter,
            missing_data_treatment=params.missing_data_treatment,
            site_cov_cutoff=params.site_cov_cutoff,
            bootstrap_replicates=params.bootstrap_replicates,
            heuristic_method=params.heuristic_method,
            initial_tree=params.initial_tree,
            gamma_categories=params.gamma_categories,
            threads=threads
        )

rule iqtree_construct_tree:
    """
    Constructs the phylogenetic tree using IQ-TREE.
    """
    input:
        FASTA = rules.create_snp_matrix.output.FASTA
    output:
        NWK = 'tree/iqtree/phylogeny.nwk',
        JSON = 'tree/iqtree/model.json'
    params:
        bootstrap_replicates = 100
    threads: 8
    run:
        from pacu.app.utils import iqtreeutils

        # Construct tree
        command = iqtreeutils.run_ml_tree_construction(
            Path(input.FASTA).absolute(), Path(output.NWK).absolute(), threads)

        # Extract selected model
        model = iqtreeutils.extract_selected_model(command.stdout)
        with open(output.JSON, 'w') as handle:
            json.dump({'model': model}, handle, indent=2)

rule select_tree:
    """
    Selects the tree based on the configuration.
    """
    input:
        NWK = rules.mega_construct_tree.output.NWK if config.get('phylogeny_method') == 'mega' else
            rules.iqtree_construct_tree.output.NWK,
        JSON = rules.mega_model_selection.output.JSON if config.get('phylogeny_method') == 'mega' else
            rules.iqtree_construct_tree.output.JSON
    output:
        NWK = 'tree/phylogeny.nwk',
        JSON = 'tree/model.json'
    shell:
        """
        cp {input.NWK} {output.NWK}
        cp {input.JSON} {output.JSON}
        """

rule snp_dists_extract:
    """
    Calculates the SNP distances based on the SNP matrix.
    """
    input:
        FASTA = rules.create_snp_matrix.output.FASTA
    output:
        TSV = temporary('tree/distances-unsorted.tsv')
    shell:
        """
        snp-dists {input.FASTA} > {output.TSV}
        """

rule snp_dists_sort:
    """
    Sorts the SNP distances based on the order in the phylogeny.
    """
    input:
        TSV = rules.snp_dists_extract.output.TSV,
        NWK = rules.select_tree.output.NWK
    output:
        TSV = 'tree/distances.tsv'
    run:
        try:
            workflowutils.sort_snp_dist_matrix(Path(input.NWK), Path(input.TSV), Path(output.TSV))
        except UnboundLocalError:
            shutil.copyfile(Path(input.TSV), Path(output.TSV))

rule visualize_tree:
    """
    Visualizes the tree using FigTree.
    """
    input:
        NWK = rules.select_tree.output.NWK
    output:
        PNG = 'tree/tree.png'
    params:
        width = config['image']['width'],
        height = config['image']['height']
    run:
        workflowutils.plot_newick_phylogeny(
            Path(input.NWK).absolute(), Path(output.PNG).absolute(), params.width, params.height)

rule create_report:
    """
    Creates the output report of the workflow.
    """
    input:
        # Statistics
        TSV_stats = rules.combine_variant_calling_stats.output.TSV,
        TSV_depth = rules.samtools_depth_combine.output.TSV,
        TSV_regions = rules.region_filtering_combine_stats.output.TSV,
        TSV_snp_pos = rules.create_snp_matrix.output.TSV,
        TSV_dist = rules.snp_dists_sort.output.TSV,
        # Tree
        FASTA = rules.create_snp_matrix.output.FASTA,
        NWK = rules.select_tree.output.NWK,
        JSON = rules.select_tree.output.JSON,
        PNG = rules.visualize_tree.output.PNG,
        # Variant calling and filtering
        VCF = expand(rules.variant_filtering_distance.output.VCF, key=config['input'].keys()),
        PNG_regions = rules.region_filtering_plot.output.PNG,
        BED = [
            rules.region_filtering_collect_low_depth_regions.output.BED,
            rules.region_filtering_merge_bed_files.output.BED,
            rules.gubbins_select_bed.output.BED]
    output:
        HTML = Path(config['output']['html'])
    params:
        dir_out = Path(config['output']['dir']),
        config = config
    run:
        from pacu.app.utils import reportutils

        # Initialize report
        report = HtmlReport(Path(output.HTML), Path(params.dir_out))
        if not Path(params.dir_out).exists():
            Path(params.dir_out).mkdir(parents=True)
        report.initialize('PACU report', Path(str(files('pacu').joinpath('resources/style.css'))))
        report.add_pipeline_header('PACU')

        # Add sections
        report.add_html_object(reportutils.create_analysis_info_section(params.config))
        report.add_html_object(reportutils.create_parameter_section(params.config))
        report.add_html_object(reportutils.create_mapping_section(
            Path(input.TSV_depth), Path(params.config['reference']['fasta'])))
        report.add_html_object(reportutils.create_variant_calling_section(Path(input.TSV_stats)))
        section_region_filt = reportutils.create_region_filtering_section(
            Path(input.TSV_regions), Path(input.PNG_regions), config.get('skip_gubbins', False))
        report.add_html_object(section_region_filt)
        section_region_filt.copy_files(Path(params.dir_out))
        with open(input.JSON) as handle:
            model = json.load(handle)['model']
        section_tree = reportutils.create_tree_section(Path(input.NWK), Path(input.PNG), Path(input.FASTA),
            Path(input.TSV_snp_pos), model)
        report.add_html_object(section_tree)
        section_tree.copy_files(Path(params.dir_out))
        report.add_html_object(reportutils.create_snp_distances_section(Path(input.TSV_dist)))
        report.add_html_object(reportutils.create_citations_section())

        # Add VCF files to the output report
        dir_vcf = report.output_dir / 'vcf'
        dir_vcf.mkdir(exist_ok=True)
        for name, path_vcf in zip(config['input'].keys(), [Path(x) for x in input.VCF]):
            shutil.copyfile(path_vcf, dir_vcf / f'{name}.vcf')

        # Add BED files to the output report
        dir_bed = report.output_dir / 'region_filtering'
        dir_bed.mkdir(exist_ok=True)
        for path_bed in [Path(x) for x in input.BED]:
            shutil.copyfile(path_bed, dir_bed / path_bed.name)

        # Save report
        report.save()
