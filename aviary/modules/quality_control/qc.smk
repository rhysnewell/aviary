ruleorder: filtlong_no_reference > link_reads
ruleorder: complete_qc_all > complete_qc_long > complete_qc_short
# ruleorder: filtlong_paired > filtlong_single
# ruleorder: filtlong_paired > filtlong_no_reference

# If you don't want to filter the reads using a genome just link them into the folder
rule link_reads:
    input:
        fastq = config["long_reads"],
    output:
        temp("data/long_reads.fastq.gz")
    threads:
        config['max_threads']
    run:
        import subprocess
        if len(input.fastq) == 1: # Check if only one longread sample
            shell("ln -s {input.fastq} {output}")
        elif len(input.fastq) > 1:
            subprocess.Popen("ln -s %s %s" % (input.fastq[0], output))
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
        min_mean_q = config['min_mean_q']
    threads:
        config['max_threads']
    benchmark:
        "benchmarks/filtlong.benchmark.txt"
    conda:
        'envs/filtlong.yaml'
    shell:
        'filtlong --min_length {params.min_length} --min_mean_q {params.min_mean_q} {input.long} '
        '| pigz -p {threads} > {output.long}'

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
        directory("www/fastqc/"),
        temp('www/fastqc/done')
    conda:
        "envs/fastqc.yaml"
    benchmark:
        "benchmarks/fastqc_short.benchmark.txt"
    threads:
        config["max_threads"]
    script:
        "scripts/run_fastqc.py"

rule fastqc_long:
    input:
        'data/long_reads.fastq.gz'
    output:
        directory("www/fastqc_long/"),
        temp('www/fastqc_long/done')
    conda:
        "envs/fastqc.yaml"
    benchmark:
        "benchmarks/fastqc_long.benchmark.txt"
    threads:
        config["max_threads"]
    shell:
        "fastqc -o www/fastqc_long/ -t {threads} {input[0]}; touch www/fastqc_long/done"

rule nanoplot:
    input:
        long = config['long_reads']
         # "data/long_reads.fastq.gz"
    output:
        directory("www/nanoplot/"),
        temp('www/nanoplot/done')
    conda:
        "envs/nanoplot.yaml"
    benchmark:
        "benchmarks/nanoplot.benchmark.txt"
    threads:
        config["max_threads"]
    shell:
        "NanoPlot -o www/nanoplot -p longReads -t {threads} --fastq {input.long}; touch www/nanoplot/done"


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