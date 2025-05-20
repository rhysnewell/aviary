
QC_SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), 'scripts')
from aviary.modules.common import pixi_run

localrules: assembly_size, assembly_quality, complete_qc_short, complete_qc_long, complete_qc_all

ruleorder: complete_qc_all > complete_qc_long > complete_qc_short

if config['fasta'] == 'none':
    config['fasta'] = 'assembly/final_contigs.fasta'

### Filter illumina reads against provided reference
rule qc_short_reads:
    input:
        short_reads_1 = config["short_reads_1"],
        short_reads_2 = config["short_reads_2"],
    output:
        bam = temp("data/short_unmapped_ref.bam"),
        fastq = "data/short_reads.fastq.gz",
        filtered = "data/short_filter.done"
    params:
        disable_adapter_trimming = config["disable_adapter_trimming"],
        min_length = config['min_short_read_size'],
        max_length = config['max_short_read_size'],
        quality_cutoff = config['quality_cutoff'],
        unqualified_percent_limit = config['unqualified_percent_limit'],
        extra_fastp_params = config['extra_fastp_params'],
        coassemble = config["coassemble"],
        reference_filter = [] if "none" in config["reference_filter"] else config["reference_filter"],
        skip_qc = config["skip_qc"]
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 512*1024*attempt),
        runtime = lambda wildcards, attempt: 8*60*attempt,
    log:
        "logs/qc_short_reads.log"
    benchmark:
        "benchmarks/qc_short_reads.benchmark.txt"
    shell:
        f'{pixi_run} -e minimap2 {QC_SCRIPTS_DIR}/'+\
        """qc_short_reads.py \
          --disable-adapter-trimming {params.disable_adapter_trimming} \
          --skip-qc {params.skip_qc} \
          --coassemble {params.coassemble} \
          --short-reads-1 {input.short_reads_1} \
          --short-reads-2 {input.short_reads_2} \
          --output-bam {output.bam} \
          --output-fastq {output.fastq} \
          --output-filtered {output.filtered} \
          --quality-cutoff {params.quality_cutoff} \
          --unqualified-percent-limit {params.unqualified_percent_limit} \
          --min-length {params.min_length} \
          --max-length {params.max_length} \
          --extra-fastp-params "{params.extra_fastp_params}" \
          --reference-filter {params.reference_filter} \
          --threads {threads} \
          --log {log}
        """

rule qc_long_reads:
    input:
        long_reads = config["long_reads"]
    output:
        output_long = temp("data/long_reads.fastq.gz")
    params:
        coassemble = config["coassemble"],
        min_length = config['min_read_size'],
        keep_percent = config['keep_percent'],
        min_mean_q = config['min_mean_q'],
        reference_filter = [] if config["reference_filter"] == "none" else config["reference_filter"],
        skip_qc = config["skip_qc"]
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 12*60*attempt,
    log:
        "logs/qc_long_reads.log"
    benchmark:
        "benchmarks/qc_long_reads.benchmark.txt"
    shell:
        f'{pixi_run} -e chopper {QC_SCRIPTS_DIR}/'+\
        """qc_long_reads.py \
            --long-reads {input.long_reads} \
            --output-long {output.output_long} \
            --coassemble {params.coassemble} \
            --min-length {params.min_length} \
            --min-quality {params.min_mean_q} \
            --keep-percent {params.keep_percent} \
            --reference-filter {params.reference_filter} \
            --skip-qc {params.skip_qc} \
            --threads {threads} \
            --log {log}
        """

rule fastqc:
    input:
        config['short_reads_1']
    output:
        output = "www/fastqc/done"
    benchmark:
        "benchmarks/fastqc_short.benchmark.txt"
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 12*60*attempt,
    log:
        "logs/fastqc.log"
    shell:
        f'{pixi_run} -e fastqc {QC_SCRIPTS_DIR}/'+\
        """run_fastqc.py \
        --short-reads-1 {config[short_reads_1]} \
        --short-reads-2 {config[short_reads_2]} \
        --threads {threads} \
        --log {log[0]}
        """

rule fastqc_long:
    input:
        'data/long_reads.fastq.gz'
    output:
        output = "www/fastqc_long/done"
    benchmark:
        "benchmarks/fastqc_long.benchmark.txt"
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 12*60*attempt,
    log:
        "logs/fastqc_long.log"
    shell:
        f"{pixi_run} -e fastqc "
        "fastqc -o www/fastqc_long/ -t {threads} {input[0]} > {log} 2>&1; touch {output.output}"

rule nanoplot:
    input:
        long = config['long_reads']
         # "data/long_reads.fastq.gz"
    output:
        output = "www/nanoplot/done"
    benchmark:
        "benchmarks/nanoplot.benchmark.txt"
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 12*60*attempt,
    log:
        "logs/nanoplot.log"
    shell:
        f"{pixi_run} -e nanoplot "
        "NanoPlot -o www/nanoplot -p longReads -t {threads} --fastq {input.long} > {log} 2>&1; touch {output.output}"

rule metaquast:
    """
    MetaQuast on input assembly (one or more). Compare against GSA if one is provided
    """
    input:
        assembly = config['fasta']
    params:
        gsa = f" -r {','.join(config['gsa'])} " if config['gsa'][0] != 'none' else ""
    output:
        "www/metaquast/report.html"
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 128*1024*attempt),
        runtime = lambda wildcards, attempt: 12*60*attempt,
    log:
        "logs/metaquast.log"
    shell:
        f"{pixi_run} -e quast "
        "metaquast.py {input.assembly} -t {threads} -o www/metaquast {params.gsa} --min-identity 80 --extensive-mis-size 20000 > {log} 2>&1"

rule read_fraction_recovered:
    """
    CoverM genome on assembly to get percentage of reads mapping
    """
    input:
        fasta = config["fasta"]
    output:
        "www/fraction_recovered/short_fraction_recovered" if config['short_reads_1'] != 'none' else "www/fraction_recovered/long_fraction_recovered"
    threads:
        config["max_threads"]
    resources:
        mem_mb = lambda wildcards, attempt: min(int(config["max_memory"])*1024, 512*1024*attempt),
        runtime = lambda wildcards, attempt: 24*60 + 24*60*attempt,
    log:
        "logs/fraction_recovered.log"
    shell:
        f'{pixi_run} -e coverm {QC_SCRIPTS_DIR}/'+\
        """fraction_recovered.py \
        --input-fasta {input.fasta} \
        --long-reads {config[long_reads]} \
        --short-reads-1 {config[short_reads_1]} \
        --short-reads-2 {config[short_reads_2]} \
        --long-read-type {config[long_read_type][0]} \
        --threads {threads} \
        --log {log}
        """

rule assembly_size:
    """
    Uses the bbmap stats.sh script to retun N/L50 values and the distribution of contig sizes
    """
    input:
        fasta = config["fasta"]
    output:
        sizes = "www/assembly_stats.txt"
    log:
        "logs/assembly_stats.log"
    shell:
        f"{pixi_run} -e bbmap "
        "stats.sh {input.fasta} > {output.sizes} 2> {log}"


rule assembly_quality:
    """
    Final return point for assembly quality statistics
    """
    input:
        # "www/metaquast/report.txt",
        "www/fraction_recovered/short_fraction_recovered" if config['short_reads_1'] != 'none' else "www/fraction_recovered/long_fraction_recovered",
        "www/assembly_stats.txt"
    log:
        temp('www/assembly_quality_stats')

rule complete_qc_short:
    input:
        'www/fastqc/done',
    output:
        temp('data/qc_done')
    shell:
        'touch data/qc_done'

rule complete_qc_long:
    input:
        'www/nanoplot/done',
        'www/fastqc_long/done',
        'data/long_reads.fastq.gz'
    output:
        temp('data/qc_done')
    shell:
        'touch data/qc_done'

rule complete_qc_all:
    input:
        'www/nanoplot/done',
        'www/fastqc/done',
        'www/fastqc_long/done',
        'data/long_reads.fastq.gz'
    output:
        temp('data/qc_done')
    shell:
        'touch data/qc_done'