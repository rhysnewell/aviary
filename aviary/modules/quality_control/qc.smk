ruleorder: filtlong_paired > filtlong_single > filtlong_reference > filtlong_no_reference


rule filtlong_no_reference:
    input:
        long = config['long_reads']
    output:
        long = "data/long_reads.fastq.gz",
    params:
        min_length = config['min_long_read_length'],
        keep_percent = config['keep_percent']
    threads:
        config['max_threads']
    conda:
        'envs/filtlong.yaml'
    shell:
        'filtlong --min_length {params.min_length} -p {params.keep_percent} {input.long} \ '
        '| pigz -p {threads} > {output.long}'

rule filtlong_reference:
    input:
        reference = config['reference_filter'],
        long = config['long_reads'],
    output:
        long = "data/long_reads.fastq.gz",
    params:
        min_length = config['min_long_read_length'],
        keep_percent = config['keep_percent']
    threads:
        config['max_threads']
    conda:
        'envs/filtlong.yaml'
    shell:
        'filtlong -a {input.reference} --min_length {params.min_length} -p {params.keep_percent} {input.long} \ '
        '| pigz -p {threads} > {output.long}'

rule filtlong_single:
    input:
        pe1 = config['short_reads_1'],
        long = config['long_reads'],
    output:
        long = "data/long_reads.fastq.gz",
    params:
        min_length = config['min_long_read_length'],
        keep_percent = config['keep_percent']
    threads:
        config['max_threads']
    conda:
        'envs/filtlong.yaml'
    shell:
        'filtlong -1 {input.pe1} --min_length {params.min_length} -p {params.keep_percent} {input.long} \ '
        '| pigz -p {threads} > {output.long}'

rule filtlong_paired:
    input:
        pe1 = config['short_reads_1'],
        pe2 = config['short_reads_2'],
        long = config['long_reads'],
    output:
        long = "data/long_reads.fastq.gz",
    params:
        min_length = config['min_long_read_length'],
        keep_percent = config['keep_percent']
    threads:
        config['max_threads']
    conda:
        'envs/filtlong.yaml'
    shell:
        'filtlong -1 {input.pe1} -2 {input.pe2} --min_length {params.min_length} -p {params.keep_percent} {input.long} \ '
        '| pigz -p {threads} > {output.long}'


rule fastqc:
    input:
        config['short_reads_1']
    output:
        directory("www/fastqc/")
    conda:
        "envs/fastqc.yaml"
    threads:
        config["max_threads"]
    script:
        "../../script/run_fastqc.py"


rule nanoplot:
    input:
        "data/long_reads.fastq.gz"
    output:
        directory("www/nanoplot/")
    conda:
        "envs/nanoplot.yaml"
    threads:
        config["max_threads"]
    shell:
        "NanoPlot -o www/nanoplot -p longReads --fastq {input}"


rule complete_qc:
    input:
        directory('www/nanoplot'),
        directory('www/fastqc'),
        'data/long_reads.fastq.gz'
    output:
        'data/qc_done'
    shell:
        'touch data/qc_done'