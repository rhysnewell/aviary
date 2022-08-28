ruleorder: filtlong_no_reference > link_reads
ruleorder: complete_qc_all > complete_qc_long > complete_qc_short
# ruleorder: filtlong_paired > filtlong_single
# ruleorder: filtlong_paired > filtlong_no_reference

if config['fasta'] == 'none':
    config['fasta'] = 'assembly/final_contigs.fasta'

# If you don't want to filter the reads using a genome just link them into the folder
rule link_reads:
    input:
        fastq = config["long_reads"],
    params:
        coassemble = config["coassemble"]
    output:
        "data/long_reads.fastq.gz"
    threads:
        config['max_threads']
    group: 'qc'
    run:
        import subprocess
        import os
        import sys
        if len(input.fastq) == 1 or isinstance(input.fastq, str): # Check if only one longread sample
            shell("ln -s {input.fastq} {output}")
        elif params.coassemble and len(input.fastq) >= 1: # Check if only one longread sample
            shell("ln -s {input.fastq[0]} {output}")
        elif len(input.fastq) > 1 and not isinstance(input.fastq, str):
            for reads in input.fastq:
                shell(f"cat {reads} >> data/long_reads.fastq.gz")
            # shell("pigz -p {threads} data/long_reads.fastq")
        else:
            shell("touch {output}")

rule filtlong_no_reference:
    input:
        long = config['long_reads']
    output:
        long = "data/long_reads.fastq.gz",
    params:
        min_length = config['min_long_read_length'],
        keep_percent = config['keep_percent'],
        min_mean_q = config['min_mean_q'],
        coassemble = config["coassemble"]
    threads:
        config['max_threads']
    benchmark:
        "benchmarks/filtlong.benchmark.txt"
    group: 'qc'
    conda:
        'envs/filtlong.yaml'
    shell:
        '''
        for long_reads in {input.long}
        do
            filtlong --min_length {params.min_length} --min_mean_q {params.min_mean_q} $long_reads | pigz -p {threads} >> {output.long}
            if ! [ {params.coassemble} ]
            then
                break
            fi
        done
        '''

# rule filtlong_reference:
#     input:
#         reference = config['reference_filter'],
#         long = config['long_reads'],
#     output:
#         long = "data/long_reads.fastq.gz",
#     params:
#         min_length = config['min_long_read_length'],
#         keep_percent = config['keep_percent'],
#         min_mean_q = config['min_mean_q']
#     threads:
#         config['max_threads']
#     conda:
#         'envs/filtlong.yaml'
#     shell:
#         'filtlong -a {input.reference} --min_length {params.min_length} --min_mean_q {params.min_mean_q} -p {params.keep_percent} {input.long} '
#         '| pigz -p {threads} > {output.long}'

# rule filtlong_single:
#     input:
#         pe1 = config['short_reads_1'],
#         long = config['long_reads'],
#     output:
#         long = "data/long_reads.fastq.gz",
#     params:
#         min_length = config['min_long_read_length'],
#         keep_percent = config['keep_percent'],
#         min_mean_q = config['min_mean_q']
#     threads:
#         config['max_threads']
#     conda:
#         'envs/filtlong.yaml'
#     shell:
#         'filtlong -1 {input.pe1} --min_length {params.min_length} --min_mean_q {params.min_mean_q} -p {params.keep_percent} {input.long} \ '
#         '| pigz -p {threads} > {output.long}'

# rule filtlong_paired:
#     input:
#         pe1 = config['short_reads_1'],
#         pe2 = config['short_reads_2'],
#         long = config['long_reads'],
#     output:
#         long = "data/long_reads.fastq.gz",
#     params:
#         min_length = config['min_long_read_length'],
#         keep_percent = config['keep_percent'],
#         min_mean_q = config['min_mean_q']
#     threads:
#         config['max_threads']
#     conda:
#         'envs/filtlong.yaml'
#     shell:
#         'filtlong -1 {input.pe1} -2 {input.pe2} --min_length {params.min_length} --min_mean_q {params.min_mean_q} -p {params.keep_percent} {input.long} '
#         '| pigz -p {threads} > {output.long}'


rule fastqc:
    input:
        config['short_reads_1']
    output:
        directory("www/fastqc/")
    conda:
        "envs/fastqc.yaml"
    benchmark:
        "benchmarks/fastqc_short.benchmark.txt"
    group: 'qc'
    threads:
        config["max_threads"]
    script:
        "scripts/run_fastqc.py"

rule fastqc_long:
    input:
        'data/long_reads.fastq.gz'
    output:
        directory("www/fastqc_long/")
    conda:
        "envs/fastqc.yaml"
    benchmark:
        "benchmarks/fastqc_long.benchmark.txt"
    group: 'qc'
    threads:
        config["max_threads"]
    shell:
        "fastqc -o www/fastqc_long/ -t {threads} {input[0]}; touch www/fastqc_long/done"

rule nanoplot:
    input:
        long = config['long_reads']
         # "data/long_reads.fastq.gz"
    output:
        directory("www/nanoplot/")
    conda:
        "envs/nanoplot.yaml"
    benchmark:
        "benchmarks/nanoplot.benchmark.txt"
    group: 'qc'
    threads:
        config["max_threads"]
    shell:
        "NanoPlot -o www/nanoplot -p longReads -t {threads} --fastq {input.long}; touch www/nanoplot/done"

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
    conda:
        "envs/quast.yaml"
    shell:
        "metaquast.py {input.assembly} -t {threads} -o www/metaquast {params.gsa} --min-identity 80 --extensive-mis-size 20000"

rule read_fraction_recovered:
    """
    CoverM genome on assembly to get percentage of reads mapping
    """
    input:
        fasta = config["fasta"]
    group: 'binning'
    output:
        "www/fraction_recovered/short_fraction_recovered" if config['short_reads_1'] != 'none' else "www/fraction_recovered/long_fraction_recovered"
    conda:
        "../../envs/coverm.yaml"
    threads:
        config["max_threads"]
    script:
        "scripts/fraction_recovered.py"

rule assembly_size:
    """
    Uses the bbmap stats.sh script to retun N/L50 values and the distribution of contig sizes
    """
    input:
        fasta = config["fasta"]
    output:
        sizes = "www/assembly_stats.txt"
    shell:
        "stats.sh {input.fasta} > {output.sizes}"


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
    group: 'qc'
    shell:
        'touch data/qc_done'

rule complete_qc_long:
    input:
        'www/nanoplot/done',
        'www/fastqc_long/done',
        'data/long_reads.fastq.gz'
    output:
        temp('data/qc_done')
    group: 'qc'
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
    group: 'qc'
    shell:
        'touch data/qc_done'