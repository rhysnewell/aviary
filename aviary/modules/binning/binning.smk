BASE_SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), '..', '..', 'scripts')
BINNING_SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), 'scripts')

from aviary.modules.common import pixi_run

localrules: vamb_skip, amber_checkm_output, recover_mags, recover_mags_no_singlem

ruleorder: prepare_binning_files_gather > prepare_binning_files
ruleorder: dereplicate_and_get_abundances_paired > dereplicate_and_get_abundances_interleaved
ruleorder: checkm_rosella > amber_checkm_output
ruleorder: checkm_metabat2 > amber_checkm_output
ruleorder: checkm_semibin > amber_checkm_output
ruleorder: checkm_completebin > amber_checkm_output

onstart:
    from snakemake.utils import min_version

    # minimum required snakemake version
    min_version("6.0")
    long_reads = config["long_reads"]
    fasta = config["fasta"]
    short_reads_1 = config["short_reads_1"]
    short_reads_2 = config["short_reads_2"]
    min_contig_size = config["min_contig_size"]
    min_bin_size = config["min_bin_size"]
    gtdbtk_folder = config["gtdbtk_folder"]
    busco_folder = config["busco_folder"]
    threads = config["max_threads"]
    ## pplacer deadlocks on too many threads
    pplacer_threads = min(48, int(config["pplacer_threads"]))

    if long_reads == "none" and short_reads_1 == "none":
        sys.exit("Need at least one of long_reads or short_reads_1")
    if long_reads != "none" and not os.path.exists(long_reads[0]):
        sys.exit("long_reads does not point to a file")
    if short_reads_1 != "none" and not os.path.exists(short_reads_1[0]):
        sys.exit("short_reads_1 does not point to a file")
    if short_reads_2 != "none" and not os.path.exists(short_reads_2[0]):
        sys.exit("short_reads_2 does not point to a file")
    if gtdbtk_folder != "none" and not os.path.exists(gtdbtk_folder):
        sys.stderr.write("gtdbtk_folder does not point to a folder\n")
    if busco_folder != "none" and not os.path.exists(busco_folder):
        sys.stderr.write("busco_folder does not point to a folder\n")
    if singlem_metapackage != "none" and not os.path.exists(singlem_metapackage):
        sys.stderr.write("singlem_metapackage does not point to a folder\n")

if config['fasta'] == 'none':
    config['fasta'] = 'assembly/final_contigs.fasta'

import os
import sys
import glob

def get_num_samples():
    """
    Get the number of samples in the config
    """
    num_samples = 0
    if config["long_reads"] != "none":
        num_samples += len(config["long_reads"])
    if config["short_reads_1"] != "none":
        num_samples += len(config["short_reads_1"])
    return num_samples

rule prepare_binning_files:
    input:
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        input_fasta = "data/large_contigs.fasta",
    output:
        maxbin_coverage = "data/maxbin.cov.list",
        metabat_coverage = "data/coverm.cov"
    params:
        tmpdir = f"--tmpdir {config['tmpdir']}" if 'tmpdir' in config and config['tmpdir'] else "",
        long_reads = config["long_reads"],
        long_read_type = config["long_read_type"][0],
        short_reads_1 = config["short_reads_1"],
        short_reads_2 = config["short_reads_2"],
        bam_cache = "data/binning_bams/",
        working_dir = "data/binning_cov/",
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 512*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60 + 24*60*attempt,
    log:
        "logs/coverm_prepare.log"
    shell:
        f'{pixi_run} -e coverm {BINNING_SCRIPTS_DIR}/'+\
        """get_coverage.py \
        --long-reads {config[long_reads]} \
        --short-reads-1 {config[short_reads_1]} \
        --short-reads-2 {config[short_reads_2]} \
        --long-read-type {config[long_read_type]} \
        --input-fasta {input.input_fasta} \
        --bam-cache {params.bam_cache} \
        --working-dir {params.working_dir} \
        --coverm-output {output.metabat_coverage} \
        --maxbin-output {output.maxbin_coverage} \
        {params.tmpdir} \
        --threads {threads} \
        --log {log[0]}
        """

def select_split_samples(wildcards, read_type):
    short = config["short_reads_1"]
    if read_type == "long":
        read_data = config["long_reads"]
    elif read_type == "short_1":
        read_data = config["short_reads_1"]
    elif read_type == "short_2":
        read_data = config["short_reads_2"]
        short = config["short_reads_2"]
    else:
        raise ValueError("Unknown read type")

    reads = []
    if config["short_reads_1"] != "none":
        reads.extend(short)
    if config["long_reads"] != "none":
        reads.extend(config["long_reads"])

    split = int(wildcards.split)
    split_size = config["coverage_samples_per_split"]
    start = split * split_size
    end = (split + 1) * split_size
    split_reads = reads[start:end]

    output = list(set(split_reads) & set(read_data))

    return output if output else "none"

rule prepare_binning_files_split:
    input:
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        input_fasta = "data/large_contigs.fasta",
    output:
        maxbin_coverage = "data/{split}/maxbin.cov.list",
        metabat_coverage = "data/{split}/coverm.cov"
    params:
        tmpdir = f"--tmpdir {config['tmpdir']}" if 'tmpdir' in config and config['tmpdir'] else "",
        long_reads = lambda wildcards: select_split_samples(wildcards, "long"),
        long_read_type = config["long_read_type"][0],
        short_reads_1 = lambda wildcards: select_split_samples(wildcards, "short_1"),
        short_reads_2 = lambda wildcards: select_split_samples(wildcards, "short_2"),
        bam_cache = "data/binning_bams/",
        working_dir = "data/{split}/binning_cov/",
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 512*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60 + 24*60*attempt,
    log:
        "logs/coverm_prepare_{split}.log"
    shell:
        f'{pixi_run} -e coverm {BINNING_SCRIPTS_DIR}/'+\
        """get_coverage.py \
        --long-reads {config[long_reads]} \
        --short-reads-1 {config[short_reads_1]} \
        --short-reads-2 {config[short_reads_2]} \
        --long-read-type {config[long_read_type]} \
        --input-fasta {input.input_fasta} \
        --bam-cache {params.bam_cache} \
        --working-dir {params.working_dir} \
        --coverm-output {output.metabat_coverage} \
        --maxbin-output {output.maxbin_coverage} \
        {params.tmpdir} \
        --threads {threads} \
        --log {log[0]}
        """

def get_number_of_splits():
    return (get_num_samples() + config["coverage_samples_per_split"] - 1) // config["coverage_samples_per_split"]

rule prepare_binning_files_gather:
    input:
        maxbin_coverages = expand("data/{split}/maxbin.cov.list", split=range(get_number_of_splits())),
        metabat_coverages = expand("data/{split}/coverm.cov", split=range(get_number_of_splits())),
    output:
        maxbin_coverage = "data/maxbin.cov.list" if config["coverage_split"] else [],
        metabat_coverage = "data/coverm.cov" if config["coverage_split"] else [],
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 16*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60*attempt,
    run:
        import pandas as pd

        maxbin_cov = []
        for maxbin in input.maxbin_coverages:
            with open(maxbin) as f:
                maxbin_cov.extend(f.readlines())
        with open(output.maxbin_coverage, "w") as f:
            f.writelines(maxbin_cov)

        # Load csvs from input.metabat_coverages and merge them on the contigName and contigLen columns, removing the totalAvgDepth column from each
        metabat_cov = []
        for metabat in input.metabat_coverages:
            cov = pd.read_csv(metabat, sep='\t')
            cov.drop(columns=["totalAvgDepth"], inplace=True)
            metabat_cov.append(cov)
        metabat_cov = pd.concat(metabat_cov, axis=1)
        metabat_cov = metabat_cov.loc[:, ~metabat_cov.columns.duplicated()]

        # Calculate totalAvgDepth as the average of all columns that don't end in -var or Variance
        depth_columns = [col for col in metabat_cov.columns if not col.endswith('-var') and not col.endswith(' Variance') and col not in ['contigName', 'contigLen']]
        metabat_cov['totalAvgDepth'] = metabat_cov[depth_columns].mean(axis=1)

        # Reorder columns to place totalAvgDepth as the third column
        cols = metabat_cov.columns.tolist()
        cols.insert(2, cols.pop(cols.index('totalAvgDepth')))
        metabat_cov = metabat_cov[cols]
        metabat_cov.to_csv(output.metabat_coverage, sep='\t', index=False)

rule get_bam_indices:
    input:
        coverage = "data/coverm.cov"
    output:
        bams = "data/binning_bams/done"
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60*attempt,
    shell:
        f"ls data/binning_bams/*.bam | {pixi_run}" + " -e coverm parallel -j 1 'samtools index -@ {threads} {{}} {{}}.bai' &&"
        "touch data/binning_bams/done"


rule maxbin2:
    input:
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta",
        maxbin_cov = ancient("data/maxbin.cov.list")
    params:
        min_contig_size = config["min_contig_size"],
        touch = "" if config["strict"] else "|| touch data/maxbin2_bins/done",
        really_done = "data/maxbin2_bins/really_done",
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 72*60 + 24*60*attempt,
    output:
        "data/maxbin2_bins/done"
    log:
        "logs/maxbin2.log"
    benchmark:
        "benchmarks/maxbin2.benchmark.txt"
    shell:
        "rm -rf data/maxbin2_bins/; "
        "mkdir -p data/maxbin2_bins && " + \
        pixi_run + " -e maxbin2 run_MaxBin.pl -contig {input.fasta} -thread {threads} -abund_list {input.maxbin_cov} "
        "-out data/maxbin2_bins/maxbin -min_contig_length {params.min_contig_size} > {log} 2>&1 "
        "&& touch {output[0]} {params.really_done} {params.touch}"


rule concoct:
    input:
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta",
        bam_done = ancient("data/binning_bams/done")
    params:
        min_contig_size = config["min_contig_size"],
        touch = "" if config["strict"] else "|| touch data/concoct_bins/done",
        really_done = "data/concoct_bins/really_done",
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 72*60 + 24*60*attempt,
    output:
        "data/concoct_bins/done"
    log:
        "logs/concoct.log"
    benchmark:
        "benchmarks/concoct.benchmark.txt"
    shell:
        f"{pixi_run} -e concoct bash -c '"
        "rm -rf data/concoct_*/; "
        "mkdir -p data/concoct_working && "
        "cut_up_fasta.py {input.fasta} -c 10000 -o 0 --merge_last -b data/concoct_working/contigs_10K.bed > data/concoct_working/contigs_10K.fa 2> {log} && "
        "concoct_coverage_table.py data/concoct_working/contigs_10K.bed data/binning_bams/*.bam > data/concoct_working/coverage_table.tsv 2>> {log} && "
        "concoct --threads {threads} -l {params.min_contig_size} --composition_file data/concoct_working/contigs_10K.fa --coverage_file data/concoct_working/coverage_table.tsv -b data/concoct_working/ 2>/dev/null && "
        "merge_cutup_clustering.py data/concoct_working/clustering_gt{params.min_contig_size}.csv > data/concoct_working/clustering_merged.csv 2>> {log} && "
        "mkdir -p data/concoct_bins && "
        "extract_fasta_bins.py {input.fasta} data/concoct_working/clustering_merged.csv --output_path data/concoct_bins/ >> {log} 2>&1 "
        "&& touch {output[0]} {params.really_done} {params.touch}"


rule vamb_jgi_filter:
    """
    vamb has to have to coverage file filtered prior to running otherwise it throws an error
    Outputs a coverage file containing no contigs smaller than minimum contig size
    """
    input:
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta",
        done = ancient("data/coverm.cov")
    output:
        vamb_bams_done = "data/coverm.filt.cov"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 16*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60*attempt,
    params:
        min_contig_size = config['min_contig_size']
    run:
        import pandas as pd
        coverm_out = pd.read_csv('data/coverm.cov', sep='\t')
        coverm_out = coverm_out[coverm_out['contigLen']>=int(params.min_contig_size)]
        coverm_out.to_csv("data/coverm.filt.cov", sep='\t', index=False)


rule filter_contigs_by_size:
    input:
        fasta = ancient(config["fasta"]),
    output:
        done = touch("data/done/filter_contigs_by_size.done"),
        fasta = "data/large_contigs.fasta",
    params:
        min_contig_size = config["min_contig_size"],
    log:
        "logs/filter_contigs_by_size.log"
    shell:
        # use seqtkit
        f"{pixi_run} -e seqkit "
        "seqkit seq -m {params.min_contig_size} {input.fasta} > {output.fasta} 2> {log}"


rule vamb:
    input:
        coverage = ancient("data/coverm.vamb.cov"),
        fasta = ancient("data/large_contigs.fasta"),
    params:
        min_bin_size = config["min_bin_size"],
        min_contig_size = config["min_contig_size"],
        touch = "" if config["strict"] else "|| touch data/vamb_bins/done",
        really_done = "data/vamb_bins/really_done",
    threads:
        min(config["max_threads"], 24)
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 32*1024*attempt),
        runtime = lambda wildcards, attempt: 48*60*attempt,
    output:
        "data/vamb_bins/done"
    log:
        "logs/vamb.log"
    benchmark:
        "benchmarks/vamb.benchmark.txt"
    shell:
        "rm -rf data/vamb_bins/; " + \
        pixi_run + " -e vamb bash -c 'OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} vamb bin default --outdir data/vamb_bins/ -p {threads} --abundance_tsv {input.coverage} --fasta {input.fasta} "
        "--minfasta {params.min_bin_size} -m {params.min_contig_size} > {log} 2>&1' "
        "&& touch {output[0]} {params.really_done} {params.touch} && mkdir -p data/vamb_bins/bins"


rule vamb_skip:
    output:
        "data/vamb_bins/skipped"
    shell:
        "mkdir -p data/vamb_bins/;"
        "touch data/vamb_bins/skipped"


rule taxvamb_abundance_tsv:
    """
    VAMB post v4 has to have a coverage file with the following format:
    A TSV file with the header being "contigname" followed by one samplename per sample,
    and the values in the TSV file being precomputed abundances.
    """
    input:
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta",
        done = ancient("data/coverm.filt.cov")
    output:
        vamb_bams_done = "data/coverm.vamb.cov"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 16*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60*attempt,
    run:
        # Short-only
        # contigName	contigLen	totalAvgDepth	assembly.fasta/wgsim.1.fq.gz.bam	assembly.fasta/wgsim.1.fq.gz.bam-var
        # Short + long
        # contigName	contigLen	totalAvgDepth	assembly.fasta/wgsim.1.fq.gz.bam	assembly.fasta/wgsim.1.fq.gz.bam-var	assembly.fasta/wgsim.css.fastq.gz Trimmed Mean	assembly.fasta/wgsim.css.fastq.gz Variance
        import pandas as pd
        coverm_out = pd.read_csv("data/coverm.filt.cov", sep='\t')
        coverm_out.drop(columns=["contigLen", "totalAvgDepth"], inplace=True)
        coverm_out.rename(columns={"contigName": "contigname"}, inplace=True)
        columns_to_drop = [col for col in coverm_out.columns if col.endswith('-var') or col.endswith(' Variance')]
        coverm_out.drop(columns=columns_to_drop, inplace=True)
        coverm_out.to_csv("data/coverm.vamb.cov", sep='\t', index=False)

rule metabuli_taxonomy:
    input:
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta",
    output:
        "data/metabuli_taxonomy/done"
    threads:
        min(config["max_threads"], 24)
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        mem_gb = lambda wildcards, attempt: min(int(config["max_memory"]), 128*attempt),
        runtime = lambda wildcards, attempt: 48*60*attempt,
    params:
        metabuli_db = config['metabuli_folder'],
    log:
        "logs/metabuli.log"
    benchmark:
        "benchmarks/metabuli.benchmark.txt"
    shell:
        "rm -rf data/metabuli_taxonomy/; "
        "mkdir -p data/metabuli_taxonomy && "
        f"{pixi_run} -e metabuli metabuli classify "
        "{input.fasta} {params.metabuli_db}/gtdb data/metabuli_taxonomy tax > {log} 2>&1 "
        "--seq-mode 1 --threads {threads} --max-ram {resources.mem_gb} "
        "&& touch {output[0]}"

rule convert_metabuli:
    input:
        "data/metabuli_taxonomy/done",
        filt_cov = ancient("data/coverm.filt.cov")
    output:
        "data/metabuli_taxonomy/taxonomy.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 16*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60*attempt,
    params:
        report = "data/metabuli_taxonomy/tax_report.tsv",
        classifications = "data/metabuli_taxonomy/tax_classifications.tsv",
    log:
        "logs/metabuli_convert.log"
    script:
        "scripts/convert_metabuli.py"

rule taxvamb:
    input:
        coverage = ancient("data/coverm.vamb.cov"),
        fasta = ancient("data/large_contigs.fasta"),
        taxonomy = "data/metabuli_taxonomy/taxonomy.tsv"
    params:
        min_bin_size = config["min_bin_size"],
        min_contig_size = config["min_contig_size"],
        gpu_flag = "--cuda" if config["request_gpu"] else "",
        pixi_env = "taxvamb-gpu" if config["request_gpu"] else "taxvamb",
        touch = "" if config["strict"] else "|| touch data/taxvamb_bins/done",
        really_done = "data/taxvamb_bins/really_done",
    threads:
        min(config["max_threads"], 24)
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 32*1024*attempt),
        runtime = lambda wildcards, attempt: 48*60*attempt,
        gpus = 1 if config["request_gpu"] else 0
    output:
        "data/taxvamb_bins/done"
    log:
        "logs/taxvamb.log"
    benchmark:
        "benchmarks/taxvamb.benchmark.txt"
    shell:
        # Specify -o since we are not doing binsplitting
        "rm -rf data/taxvamb_bins/; " + \
        pixi_run + " -e {params.pixi_env} bash -c 'OPENBLAS_NUM_THREADS={threads} OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} vamb bin taxvamb --outdir data/taxvamb_bins/ -p {threads} --fasta {input.fasta} "
        "--abundance_tsv {input.coverage} --taxonomy {input.taxonomy} "
        "--minfasta {params.min_bin_size} -m {params.min_contig_size} {params.gpu_flag} -o > {log} 2>&1' "
        "&& touch {output[0]} {params.really_done} {params.touch} && mkdir -p data/vamb_bins/bins"


rule metabat2:
    input:
        coverage = ancient("data/coverm.cov"),
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta"
    params:
        min_contig_size = max(int(config["min_contig_size"]), 1500),
        min_bin_size = config["min_bin_size"],
        touch = "" if config["strict"] else "|| touch data/metabat_bins_2/done",
        really_done = "data/metabat_bins_2/really_done",
    output:
        metabat_done = "data/metabat_bins_2/done"
    threads:
        min(config["max_threads"], 24)
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 32*1024*attempt),
        runtime = lambda wildcards, attempt: 12*60*attempt,
    log:
        "logs/metabat2.log"
    benchmark:
        "benchmarks/metabat_2.benchmark.txt"
    shell:
        "rm -rf data/metabat_bins_2/; " + \
        pixi_run + " -e metabat2 metabat -t {threads} -m {params.min_contig_size} -s {params.min_bin_size} --seed 89 -i {input.fasta} "
        "-a {input.coverage} -o data/metabat_bins_2/binned_contigs > {log} 2>&1 "
        "&& touch {output[0]} {params.really_done} {params.touch}"


rule metabat_spec:
    input:
        coverage = ancient("data/coverm.cov"),
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta"
    output:
        'data/metabat_bins_spec/done'
    params:
        min_contig_size = max(int(config["min_contig_size"]), 1500),
        min_bin_size = config["min_bin_size"],
        touch = "" if config["strict"] else "|| touch data/metabat_bins_spec/done",
        really_done = "data/metabat_bins_spec/really_done",
    log:
        "logs/metabat_spec.log"
    benchmark:
        "benchmarks/metabat_spec.benchmark.txt"
    threads:
        min(config["max_threads"], 24)
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60*attempt,
    shell:
        "rm -rf data/metabat_bins_spec; " + \
        pixi_run + " -e metabat2 metabat1 -t {threads} -m {params.min_contig_size} -s {params.min_bin_size} --seed 89 --specific -i {input.fasta} "
        "-a {input.coverage} -o data/metabat_bins_spec/binned_contigs > {log} 2>&1 "
        "&& touch {output[0]} {params.really_done} {params.touch}"

rule metabat_sspec:
    input:
        coverage = ancient("data/coverm.cov"),
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta"
    output:
        'data/metabat_bins_sspec/done'
    params:
        min_contig_size = max(int(config["min_contig_size"]), 1500),
        min_bin_size = config["min_bin_size"],
        touch = "" if config["strict"] else "|| touch data/metabat_bins_sspec/done",
        really_done = "data/metabat_bins_sspec/really_done",
    log:
        "logs/metabat_sspec.log"
    benchmark:
        "benchmarks/metabat_sspec.benchmark.txt"
    threads:
        min(config["max_threads"], 24)
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60*attempt,
    shell:
        "rm -rf data/metabat_bins_sspec; " + \
        pixi_run + " -e metabat2 metabat1 -t {threads} -m {params.min_contig_size} -s {params.min_bin_size} --seed 89 --superspecific "
        "-i {input.fasta} -a {input.coverage} -o data/metabat_bins_sspec/binned_contigs > {log} 2>&1 "
        "&& touch {output[0]} {params.really_done} {params.touch}"

rule metabat_sens:
    input:
        coverage = ancient("data/coverm.cov"),
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta"
    output:
        'data/metabat_bins_sens/done'
    params:
        min_contig_size = max(int(config["min_contig_size"]), 1500),
        min_bin_size = config["min_bin_size"],
        touch = "" if config["strict"] else "|| touch data/metabat_bins_sens/done",
        really_done = "data/metabat_bins_sens/really_done",
    log:
        "logs/metabat_sens.log"
    benchmark:
        "benchmarks/metabat_sens.benchmark.txt"
    threads:
        min(config["max_threads"], 24)
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60*attempt,
    shell:
        "rm -rf data/metabat_bins_sens; " + \
        pixi_run + " -e metabat2 metabat1 -t {threads} -m {params.min_contig_size} -s {params.min_bin_size} --seed 89 --sensitive "
        "-i {input.fasta} -a {input.coverage} -o data/metabat_bins_sens/binned_contigs > {log} 2>&1 "
        "&& touch {output[0]} {params.really_done} {params.touch}"

rule metabat_ssens:
    input:
        coverage = ancient("data/coverm.cov"),
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta"
    output:
        'data/metabat_bins_ssens/done'
    params:
        min_contig_size = max(int(config["min_contig_size"]), 1500),
        min_bin_size = config["min_bin_size"],
        touch = "" if config["strict"] else "|| touch data/metabat_bins_ssens/done",
        really_done = "data/metabat_bins_ssens/really_done",
    log:
        "logs/metabat_ssens.log"
    benchmark:
        "benchmarks/metabat_ssens.benchmark.txt"
    threads:
        min(config["max_threads"], 24)
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60*attempt,
    shell:
        "rm -rf data/metabat_bins_ssens; " + \
        pixi_run + " -e metabat2 metabat1 -t {threads} -m {params.min_contig_size} -s {params.min_bin_size} --seed 89 --supersensitive "
        "-i {input.fasta} -a {input.coverage} -o data/metabat_bins_ssens/binned_contigs > {log} 2>&1 "
        "&& touch {output[0]} {params.really_done} {params.touch}"

rule rosella:
    """
    Runs Rosella.
    """
    input:
        coverage = ancient("data/coverm.cov"),
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta"
    params:
        min_contig_size = config["min_contig_size"],
        min_bin_size = config["min_bin_size"],
        touch = "" if config["strict"] else "|| touch data/rosella_bins/done",
        really_done = "data/rosella_bins/really_done",
    output:
        # kmers = "data/rosella_bins/kmer_frequencies.tsv",
        done = "data/rosella_bins/done"
    threads:
        min(config["max_threads"], 24)
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60*attempt,
    log:
        "logs/rosella.log"
    benchmark:
        "benchmarks/rosella.benchmark.txt"
    shell:
        "rm -rf data/rosella_bins/; " + \
        pixi_run + " -e rosella rosella recover -r {input.fasta} -C {input.coverage} -t {threads} -o data/rosella_bins "
        "--min-contig-size {params.min_contig_size} --min-bin-size {params.min_bin_size} --n-neighbors 100 > {log} 2>&1 "
        "&& touch {output[0]} {params.really_done} {params.touch}"

rule semibin:
    input:
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta",
        bams_indexed = ancient("data/binning_bams/done")
    params:
        # Can't use premade model with multiple samples, so disregard if provided
        semibin_model = f"--environment {config['semibin_model']} " if get_num_samples() == 1 else "",
        semibin_sequencing_type = "--sequencing-type=long_read" if config["long_reads"] != "none" else "",
        pixi_env = "semibin-gpu" if config["request_gpu"] else "semibin",
        touch = "" if config["strict"] else "|| touch data/semibin_bins/done",
        really_done = "data/semibin_bins/really_done",
    output:
        "data/semibin_bins/done"
    threads:
        min(config["max_threads"], 24)
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60 + 48*60*(attempt-1),
        gpus = 1 if config["request_gpu"] else 0
    log:
        "logs/semibin.log"
    benchmark:
        "benchmarks/semibin.benchmark.txt"
    shell:
        pixi_run + " -e {params.pixi_env} bash -c '"
        "rm -rf data/semibin_bins/; "
        "mkdir -p data/semibin_bins/output_bins/ && "
        "SemiBin2 single_easy_bin "
        "-i {input.fasta} "
        "-b data/binning_bams/*.bam "
        "-o data/semibin_bins "
        "{params.semibin_model} "
        "-p {threads} "
        "--self-supervised "
        "--compression none "
        "{params.semibin_sequencing_type} "
        "> {log} 2>&1 "
        "&& touch {output[0]} {params.really_done} {params.touch}'"


rule comebin:
    input:
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta",
        bams_indexed = ancient("data/binning_bams/done")
    output:
        done = "data/comebin_bins/done"
    threads:
        config["max_threads"]
    params:
        pixi_env = "comebin-gpu" if config["request_gpu"] else "comebin",
        touch = "" if config["strict"] else "|| touch data/comebin_bins/done",
        really_done = "data/comebin_bins/really_done",
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60*attempt,
        gpus = 1 if config["request_gpu"] else 0
    params:
        pixi_env = "comebin-gpu" if config["request_gpu"] else "comebin"
    log:
        "logs/comebin.log"
    benchmark:
        "benchmarks/comebin.benchmark.txt"
    shell:
        "rm -rf data/comebin_bins/; " + \
        pixi_run + " -e {params.pixi_env} run_comebin.sh -a {input.fasta} -p data/binning_bams -t {threads} -o data/comebin_bins > {log} 2>&1 "
        "&& touch {output[0]} {params.really_done} {params.touch}"


rule completebin:
    input:
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta",
        bams_indexed = ancient("data/binning_bams/done")
    output:
        done = "data/completebin_bins/done"
    threads:
        config["max_threads"]
    params:
        pixi_env = "completebin-gpu" if config["request_gpu"] else "completebin",
        touch = "" if config["strict"] else "|| touch data/completebin_bins/done",
        really_done = "data/completebin_bins/really_done",
        device = "cuda:0" if config["request_gpu"] else "cpu",
        min_contig_length = max(900, config["min_contig_size"]),  # CompleteBin minimum is 900bp by default
        batch_size = 512 if config["request_gpu"] else 128,  # Reduce batch size for GPU memory management
        db_param = f"--db_files_path {config['completebin_db']}" if "completebin_db" in config and config["completebin_db"] else "",
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 256*1024*attempt),  # CompleteBin needs more memory
        runtime = lambda wildcards, attempt: 48*60*attempt,  # Allow longer runtime
        gpus = 1 if config["request_gpu"] else 0
    log:
        "logs/completebin.log"
    benchmark:
        "benchmarks/completebin.benchmark.txt"
    shell:
        "rm -rf data/completebin_bins/; " + \
        "mkdir -p data/completebin_bins/temp/ && " + \
        pixi_run + " -e {params.pixi_env} bash -c '"
        "completebin "
        "-c {input.fasta} "
        "-b data/binning_bams/*.bam "
        "-o data/completebin_bins/ "
        "-temp data/completebin_bins/temp/ "
        "--device {params.device} "
        "--min_contig_length {params.min_contig_length} "
        "--batch_size {params.batch_size} "
        "--num_workers {threads} "
        "{params.db_param} "
        "> {log} 2>&1 "
        "&& touch {output[0]} {params.really_done} {params.touch}'"


rule checkm_rosella:
    input:
        done = ancient("data/rosella_bins/done")
    params:
        pplacer_threads = lambda wildcards, threads: min(threads, config["pplacer_threads"]),
        checkm2_db_path = config["checkm2_db_folder"],
        bin_folder = "data/rosella_bins/",
        extension = "fna",
        refinery_max_iterations = config["refinery_max_iterations"],
    output:
        output_folder = directory("data/rosella_bins/checkm2_out/"),
        output_file = "data/rosella_bins/checkm.out"
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 8*60*attempt,
    log:
        "logs/checkm_rosella.log"
    shell:
        f'{pixi_run} -e checkm2 {BINNING_SCRIPTS_DIR}/'+\
        """run_checkm.py \
        --checkm2-db {params.checkm2_db_path} \
        --bin-folder {params.bin_folder} \
        --bin-ext {params.extension} \
        --refinery-max-iterations {params.refinery_max_iterations} \
        --output-folder {output.output_folder} \
        --output-file {output.output_file} \
        --threads {threads} \
        --log {log}
        """

rule checkm_metabat2:
    input:
        done = ancient("data/metabat_bins_2/done")
    params:
        pplacer_threads = lambda wildcards, threads: min(threads, config["pplacer_threads"]),
        checkm2_db_path = config["checkm2_db_folder"],
        bin_folder = "data/metabat_bins_2/",
        extension = "fa",
        refinery_max_iterations = config["refinery_max_iterations"],
    output:
        output_folder = directory("data/metabat_bins_2/checkm2_out/"),
        output_file = "data/metabat_bins_2/checkm.out"
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 8*60*attempt,
    log:
        "logs/checkm_metabat2.log"
    shell:
        f'{pixi_run} -e checkm2 {BINNING_SCRIPTS_DIR}/'+\
        """run_checkm.py \
        --checkm2-db {params.checkm2_db_path} \
        --bin-folder {params.bin_folder} \
        --bin-ext {params.extension} \
        --refinery-max-iterations {params.refinery_max_iterations} \
        --output-folder {output.output_folder} \
        --output-file {output.output_file} \
        --threads {threads} \
        --log {log}
        """

rule checkm_semibin:
    input:
        done = ancient("data/semibin_bins/done")
    params:
        pplacer_threads = lambda wildcards, threads: min(threads, config["pplacer_threads"]),
        checkm2_db_path = config["checkm2_db_folder"],
        bin_folder = "data/semibin_bins/output_bins/",
        extension = "fa",
        refinery_max_iterations = config["refinery_max_iterations"],
    output:
        output_folder = directory("data/semibin_bins/checkm2_out/"),
        output_file = "data/semibin_bins/checkm.out"
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 8*60*attempt,
    log:
        "logs/checkm_semibin.log"
    shell:
        f'{pixi_run} -e checkm2 {BINNING_SCRIPTS_DIR}/'+\
        """run_checkm.py \
        --checkm2-db {params.checkm2_db_path} \
        --bin-folder {params.bin_folder} \
        --bin-ext {params.extension} \
        --refinery-max-iterations {params.refinery_max_iterations} \
        --output-folder {output.output_folder} \
        --output-file {output.output_file} \
        --threads {threads} \
        --log {log}
        """

rule checkm_completebin:
    input:
        done = ancient("data/completebin_bins/done")
    params:
        pplacer_threads = lambda wildcards, threads: min(threads, config["pplacer_threads"]),
        checkm2_db_path = config["checkm2_db_folder"],
        bin_folder = "data/completebin_bins/",
        extension = "fasta",
        refinery_max_iterations = config["refinery_max_iterations"],
    output:
        output_folder = directory("data/completebin_bins/checkm2_out/"),
        output_file = "data/completebin_bins/checkm.out"
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 8*60*attempt,
    log:
        "logs/checkm_completebin.log"
    shell:
        f'{pixi_run} -e checkm2 {BINNING_SCRIPTS_DIR}/'+\
        """run_checkm.py \
        --checkm2-db {params.checkm2_db_path} \
        --bin-folder {params.bin_folder} \
        --bin-ext {params.extension} \
        --refinery-max-iterations {params.refinery_max_iterations} \
        --output-folder {output.output_folder} \
        --output-file {output.output_file} \
        --threads {threads} \
        --log {log}
        """

rule refine_rosella:
    input:
        checkm = ancient('data/rosella_bins/checkm.out'),
        rosella = ancient('data/rosella_bins/done'),
        coverage = ancient("data/coverm.cov"),
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta",
        # kmers = "data/rosella_bins/kmer_frequencies.tsv"
    output:
        'data/rosella_refined/done'
    benchmark:
        'benchmarks/refine_rosella.benchmark.txt'
    params:
        bin_folder = "data/rosella_bins/",
        extension = "fna",
        output_folder = "data/rosella_refined/",
        min_bin_size = config["min_bin_size"],
        max_iterations = config["refinery_max_iterations"],
        max_retries = config["refinery_max_retries"],
        pplacer_threads = lambda wildcards, threads: min(threads, config["pplacer_threads"]),
        max_contamination = 15,
        final_refining = False,
        bin_prefix = "rosella"
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 48*60 + 24*60*attempt,
    log:
        "logs/refine_rosella.log"
    shell:
        f'{pixi_run} -e checkm2 {BINNING_SCRIPTS_DIR}/'+\
        """rosella_refine.py \
        --checkm {input.checkm} \
        --coverage {input.coverage} \
        --fasta {input.fasta} \
        --output-folder {params.output_folder} \
        --final-refining {params.final_refining} \
        --min-bin-size {params.min_bin_size} \
        --max-iterations {params.max_iterations} \
        --max-retries {params.max_retries} \
        --pplacer-threads {params.pplacer_threads} \
        --threads {threads} \
        --max-contamination {params.max_contamination} \
        --bin-folder {params.bin_folder} \
        --extension {params.extension} \
        --bin-prefix {params.bin_prefix} \
        --log {log}
        """

rule refine_metabat2:
    input:
        checkm = ancient('data/metabat_bins_2/checkm.out'),
        rosella = ancient('data/metabat_bins_2/done'),
        coverage = ancient("data/coverm.cov"),
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta",
        # kmers = "data/rosella_bins/kmer_frequencies.tsv"
    output:
        'data/metabat2_refined/done'
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 48*60 + 24*60*attempt,
    benchmark:
        'benchmarks/refine_metabat2.benchmark.txt'
    params:
        bin_folder = "data/metabat_bins_2/",
        extension = "fa",
        output_folder = "data/metabat2_refined/",
        min_bin_size = config["min_bin_size"],
        max_iterations = config["refinery_max_iterations"],
        max_retries = config["refinery_max_retries"],
        pplacer_threads = lambda wildcards, threads: min(threads, config["pplacer_threads"]),
        max_contamination = 15,
        final_refining = False,
        bin_prefix = "metabat2"
    log:
        "logs/refine_metabat2.log"
    shell:
        f'{pixi_run} -e checkm2 {BINNING_SCRIPTS_DIR}/'+\
        """rosella_refine.py \
        --checkm {input.checkm} \
        --coverage {input.coverage} \
        --fasta {input.fasta} \
        --output-folder {params.output_folder} \
        --final-refining {params.final_refining} \
        --min-bin-size {params.min_bin_size} \
        --max-iterations {params.max_iterations} \
        --max-retries {params.max_retries} \
        --pplacer-threads {params.pplacer_threads} \
        --threads {threads} \
        --max-contamination {params.max_contamination} \
        --bin-folder {params.bin_folder} \
        --extension {params.extension} \
        --bin-prefix {params.bin_prefix} \
        --log {log}
        """

rule refine_semibin:
    input:
        checkm = ancient('data/semibin_bins/checkm.out'),
        rosella = ancient('data/semibin_bins/done'),
        coverage = ancient("data/coverm.cov"),
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta",
        # kmers = "data/rosella_bins/kmer_frequencies.tsv"
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 48*60 + 24*60*attempt,
    output:
        touch('data/semibin_refined/done')
    benchmark:
        'benchmarks/refine_semibin.benchmark.txt'
    params:
        bin_folder = "data/semibin_bins/output_bins/",
        extension = "fa",
        output_folder = "data/semibin_refined/",
        min_bin_size = config["min_bin_size"],
        max_iterations = config["refinery_max_iterations"],
        max_retries = config["refinery_max_retries"],
        pplacer_threads = lambda wildcards, threads: min(threads, config["pplacer_threads"]),
        max_contamination = 15,
        final_refining = False,
        bin_prefix = "semibin2"
    log:
        "logs/refine_semibin.log"
    shell:
        f'{pixi_run} -e checkm2 {BINNING_SCRIPTS_DIR}/'+\
        """rosella_refine.py \
        --checkm {input.checkm} \
        --coverage {input.coverage} \
        --fasta {input.fasta} \
        --output-folder {params.output_folder} \
        --final-refining {params.final_refining} \
        --min-bin-size {params.min_bin_size} \
        --max-iterations {params.max_iterations} \
        --max-retries {params.max_retries} \
        --pplacer-threads {params.pplacer_threads} \
        --threads {threads} \
        --max-contamination {params.max_contamination} \
        --bin-folder {params.bin_folder} \
        --extension {params.extension} \
        --bin-prefix {params.bin_prefix} \
        --log {log}
        """

rule refine_completebin:
    input:
        checkm = ancient('data/completebin_bins/checkm.out'),
        rosella = ancient('data/completebin_bins/done'),
        coverage = ancient("data/coverm.cov"),
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta",
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 48*60 + 24*60*attempt,
    output:
        touch('data/completebin_refined/done')
    benchmark:
        'benchmarks/refine_completebin.benchmark.txt'
    params:
        bin_folder = "data/completebin_bins/",
        extension = "fasta",
        output_folder = "data/completebin_refined/",
        min_bin_size = config["min_bin_size"],
        max_iterations = config["refinery_max_iterations"],
        max_retries = config["refinery_max_retries"],
        pplacer_threads = lambda wildcards, threads: min(threads, config["pplacer_threads"]),
        max_contamination = 15,
        final_refining = False,
        bin_prefix = "completebin"
    log:
        "logs/refine_completebin.log"
    shell:
        f'{pixi_run} -e checkm2 {BINNING_SCRIPTS_DIR}/'+\
        """rosella_refine.py \
        --checkm {input.checkm} \
        --coverage {input.coverage} \
        --fasta {input.fasta} \
        --output-folder {params.output_folder} \
        --final-refining {params.final_refining} \
        --min-bin-size {params.min_bin_size} \
        --max-iterations {params.max_iterations} \
        --max-retries {params.max_retries} \
        --pplacer-threads {params.pplacer_threads} \
        --threads {threads} \
        --max-contamination {params.max_contamination} \
        --bin-folder {params.bin_folder} \
        --extension {params.extension} \
        --bin-prefix {params.bin_prefix} \
        --log {log}
        """

rule amber_checkm_output:
    input:
        amber_done = "data/amber_refine/for_refine/index.html"
    output:
        metabat_checkm = 'data/metabat_bins_2/checkm.out',
        rosella_checkm = 'data/rosella_bins/checkm.out',
        semibin_checkm = 'data/semibin_bins/checkm.out',
        completebin_checkm = 'data/completebin_bins/checkm.out'
        # dastool_checkm = 'data/das_tool_bins_with_refine/checkm.out'
    run:
        import pandas as pd

        def amber_to_checkm_like(amber_table, prefix="data/rosella_bins/", output="data/rosella_bins/checkm.out", extension="fna"):
            amber_table["Bin Id"] = amber_table["Bin ID"].replace(prefix, "", regex=True)
            amber_table["Bin Id"] = amber_table["Bin Id"].replace(f".{extension}", "", regex=True)
            amber_table["Completeness"] = amber_table["Completeness (bp)"] * 100
            amber_table["Contamination"] = round((1 - amber_table["Purity (bp)"]) * 100, 2)
            amber_table.to_csv(output, sep='\t', index=False)

        try:
            rosella_amber = pd.read_csv("data/amber_refine/for_refine/genome/rosella_amber.tsv/metrics_per_bin.tsv", sep='\t')
            amber_to_checkm_like(rosella_amber, "data/rosella_bins/", "data/rosella_bins/checkm.out", "fna")
        except FileNotFoundError:
            pass

        try:
            metabat_amber = pd.read_csv("data/amber_refine/for_refine/genome/m2_amber.tsv/metrics_per_bin.tsv", sep='\t')
            amber_to_checkm_like(metabat_amber, "data/metabat_bins_2/", "data/metabat_bins_2/checkm.out", "fa")
        except FileNotFoundError:
            pass

        try:
            semibin_amber = pd.read_csv("data/amber_refine/for_refine/genome/semibin_amber.tsv/metrics_per_bin.tsv", sep='\t')
            amber_to_checkm_like(semibin_amber, "data/semibin_bins/output_bins/", "data/semibin_bins/checkm.out", "fa")
        except FileNotFoundError:
            pass

        try:
            completebin_amber = pd.read_csv("data/amber_refine/for_refine/genome/completebin_amber.tsv/metrics_per_bin.tsv", sep='\t')
            amber_to_checkm_like(completebin_amber, "data/completebin_bins/", "data/completebin_bins/checkm.out", "fasta")
        except FileNotFoundError:
            pass


rule das_tool:
    """
    Runs dasTool on the output of all binning algorithms. If a binner failed to produce bins then their output is ignored
    """
    input:
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta",
        metabat2_done = [] if "metabat2" in config["skip_binners"] else "data/metabat2_refined/done",
        concoct_done = [] if "concoct" in config["skip_binners"] else "data/concoct_bins/done",
        maxbin_done = [] if "maxbin2" in config["skip_binners"] else "data/maxbin2_bins/done",
        metabat_sspec = [] if "metabat_sspec" in config["skip_binners"] else "data/metabat_bins_sspec/done",
        metabat_spec = [] if "metabat_spec" in config["skip_binners"] else "data/metabat_bins_spec/done",
        metabat_ssens = [] if "metabat_ssens" in config["skip_binners"] else "data/metabat_bins_ssens/done",
        metabat_sense = [] if "metabat_sens" in config["skip_binners"] else "data/metabat_bins_sens/done",
        rosella_done = [] if "rosella" in config["skip_binners"] else "data/rosella_refined/done",
        semibin_done = [] if "semibin" in config["skip_binners"] else "data/semibin_refined/done",
        comebin_done = [] if "comebin" in config["skip_binners"] else "data/comebin_bins/done",
        completebin_done = [] if "completebin" in config["skip_binners"] else "data/completebin_refined/done",
        vamb_done = [] if "vamb" in config["skip_binners"] else "data/vamb_bins/done",
        taxvamb_done = [] if "taxvamb" in config["skip_binners"] else "data/taxvamb_bins/done",
    threads:
        min(config["max_threads"], 10)
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 12*60*attempt,
    output:
        touch("data/das_tool_bins_pre_refine/done")
    log:
        "logs/das_tool.log"
    benchmark:
        "benchmarks/das_tool.benchmark.txt"
    shell:
        f'{pixi_run} -e das-tool {BINNING_SCRIPTS_DIR}/'+\
        """das_tool.py \
        --skip-binners {config[skip_binners]} \
        --fasta {input.fasta} \
        --threads {threads} \
        --log {log}
        """

rule refine_dastool:
    input:
        checkm = 'data/das_tool_bins_pre_refine/checkm.out',
        das_tool = 'data/das_tool_bins_pre_refine/done',
        coverage = ancient("data/coverm.cov"),
        large_contigs_done = "data/done/filter_contigs_by_size.done",
        fasta = "data/large_contigs.fasta",
        # kmers = "data/rosella_bins/kmer_frequencies.tsv"
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 48*60 + 24*60*attempt,
    output:
        temporary('bins/checkm.out'),
        directory('bins/final_bins')
    benchmark:
        'benchmarks/refine_dastool.benchmark.txt'
    params:
        bin_folder = "data/das_tool_bins_pre_refine/das_tool_DASTool_bins/",
        extension = "fa",
        output_folder = "data/refined_bins/",
        min_bin_size = config["min_bin_size"],
        max_iterations = config["refinery_max_iterations"],
        max_retries = config["refinery_max_retries"],
        pplacer_threads = lambda wildcards, threads: min(threads, config["pplacer_threads"]),
        max_contamination = 15,
        final_refining = True,
        bin_prefix = "dastool"
    log:
        "logs/refine_dastool.log"
    shell:
        f'{pixi_run} -e checkm2 {BINNING_SCRIPTS_DIR}/'+\
        """rosella_refine.py \
        --checkm {input.checkm} \
        --coverage {input.coverage} \
        --fasta {input.fasta} \
        --output-folder {params.output_folder} \
        --final-refining {params.final_refining} \
        --min-bin-size {params.min_bin_size} \
        --max-iterations {params.max_iterations} \
        --max-retries {params.max_retries} \
        --pplacer-threads {params.pplacer_threads} \
        --threads {threads} \
        --max-contamination {params.max_contamination} \
        --bin-folder {params.bin_folder} \
        --extension {params.extension} \
        --bin-prefix {params.bin_prefix} \
        --log {log}
        """

rule get_abundances:
    input:
        "bins/checkm.out"
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 512*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60 + 24*60*attempt,
    output:
        "data/coverm_abundances.tsv"
    log:
        "logs/coverm_abundances.log"
    shell:
        f'{pixi_run} -e coverm {BINNING_SCRIPTS_DIR}/'+\
        """get_abundances.py \
        --long-reads {config[long_reads]} \
        --short-reads-1 {config[short_reads_1]} \
        --short-reads-2 {config[short_reads_2]} \
        --long-read-type {config[long_read_type]} \
        --threads {threads} \
        --strain-analysis {config[strain_analysis]} \
        --log {log}
        """

rule finalise_stats:
    input:
        checkm1_done = "bins/checkm.out",
        checkm2_done = "bins/checkm2_output/quality_report.tsv",
        coverage_file = "data/coverm_abundances.tsv" if not config["skip_abundances"] else [],
        gtdbtk_done = "data/gtdbtk/done" if not config["skip_taxonomy"] else []
    output:
        bin_stats = "bins/bin_info.tsv",
        checkm_minimal = "bins/checkm_minimal.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 16*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60*attempt,
    script:
        "scripts/finalise_stats.py"



rule checkm_das_tool:
    input:
        done = "data/das_tool_bins_pre_refine/done"
    params:
        pplacer_threads = config["pplacer_threads"]
    output:
        "data/das_tool_bins_pre_refine/checkm.out"
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 8*60*attempt,
    log:
        "logs/checkm_das_tool.log"
    shell:
        f"{pixi_run} -e checkm bash -c '"
        "checkm lineage_wf -t {threads} --pplacer_threads {params.pplacer_threads} "
        "-x fa data/das_tool_bins_pre_refine/das_tool_DASTool_bins data/das_tool_bins_pre_refine/checkm --tab_table "
        "-f data/das_tool_bins_pre_refine/checkm.out > {log} 2>&1 && "
        "checkm qa -o 2 --tab_table -f data/das_tool_bins_pre_refine/checkm.out "
        "data/das_tool_bins_pre_refine/checkm/lineage.ms data/das_tool_bins_pre_refine/checkm/ >> {log} 2>&1; "
        "'"


rule singlem_pipe_reads:
    output:
        "data/singlem_out/metagenome.combined_otu_table.csv"
    params:
        package_path = config['singlem_metapackage']
    threads: min(config["max_threads"], 48)
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 8*1024*attempt),
        runtime = lambda wildcards, attempt: 12*60*attempt,
    log:
        "logs/singlem_pipe_reads_log.txt"
    shell:
        f'{pixi_run} -e singlem {BASE_SCRIPTS_DIR}/'+\
        """singlem_reads.py \
        --long-reads {config[long_reads]} \
        --short-reads-1 {config[short_reads_1]} \
        --short-reads-2 {config[short_reads_2]} \
        --threads {threads} \
        --log {log} \
        --package-path {params.package_path}
        """

rule singlem_appraise:
    input:
        pipe_results = "data/singlem_out/metagenome.combined_otu_table.csv",
        assembly = config["fasta"],
        # gtdbtk_done = "data/gtdbtk/done",
        bins_complete = "bins/checkm.out"
    output:
        binned = "data/singlem_out/binned.otu_table.csv",
        unbinned = "data/singlem_out/unbinned.otu_table.csv",
        plot = "data/singlem_out/singlem_appraise.svg",
        assembled = "data/singlem_out/assembled.otu_table.csv",
        singlem = "data/singlem_out/singlem_appraisal.tsv"
    params:
        package_path = config['singlem_metapackage'],
        genomes_folder = "data/refined_bins/final_bins/"
    threads: min(config["max_threads"], 48)
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 8*1024*attempt),
        runtime = lambda wildcards, attempt: 12*60*attempt,
    log:
        "logs/singlem_appraise_log.txt"
    shell:
        f'{pixi_run} -e singlem {BASE_SCRIPTS_DIR}/'+\
        """singlem_appraise.py \
        --assembly {input.assembly} \
        --genomes-folder {params.genomes_folder} \
        --pipe-results {input.pipe_results} \
        --threads {threads} \
        --package-path {params.package_path} \
        --log {log}
        """

rule recover_mags:
    input:
        final_bins = "bins/bin_info.tsv",
        gtdbtk = "data/gtdbtk/done" if not config["skip_taxonomy"] else [],
        coverm = "data/coverm_abundances.tsv" if not config["skip_abundances"] else [],
        contig_coverage = "data/coverm.cov",
        singlem = "data/singlem_out/singlem_appraisal.tsv" if not config["skip_singlem"] else [],
    output:
        bins = "bins/done"
    threads:
        config["max_threads"]
    script:
        "scripts/finalise_recovery.py"

rule recover_mags_no_singlem:
    input:
        final_bins = "bins/bin_info.tsv",
        gtdbtk = [],
        coverm = "data/coverm_abundances.tsv" if not config["skip_abundances"] else [],
        contig_coverage = "data/coverm.cov",
        singlem = [],
    output:
        bins = "bins/done",
    threads:
        config["max_threads"]
    script:
        "scripts/finalise_recovery.py"

# Special rule to help out with a buggy output
rule dereplicate_and_get_abundances_paired:
    input:
        pe_1 = config["short_reads_1"],
        pe_2 = config["short_reads_2"]
    output:
        output_abundances = 'bins/coverm_abundances.tsv'
    params:
        final_bins = 'bins/final_bins',
        derep_ani = 0.97
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 512*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60 + 24*60*attempt,
    log:
        "logs/coverm_abundances_paired.log"
    shell:
        pixi_run + " -e coverm "
        "coverm genome -t {threads} -d bins/final_bins/ -1 {input.pe_1} -2 {input.pe_2} --min-covered-fraction 0.0 -x fna > bins/coverm_abundances.tsv 2> {log}; "

# Special rule to help out with a buggy output
rule dereplicate_and_get_abundances_interleaved:
    input:
        pe_1 = config["short_reads_1"],
    output:
        output_abundances = 'bins/coverm_abundances.tsv'
    params:
        final_bins = 'bins/final_bins',
        derep_ani = 0.97
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 512*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60 + 24*60*attempt,
    log:
        "logs/coverm_abundances_interleaved.log"
    shell:
        pixi_run + " -e coverm coverm genome -t {threads} -d bins/final_bins/ --interleaved {input.pe_1} --min-covered-fraction 0.0 -x fna > bins/coverm_abundances.tsv 2> {log}; "
