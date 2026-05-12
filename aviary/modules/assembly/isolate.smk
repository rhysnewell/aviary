import os
from aviary.modules.common import pixi_run
ISOLATE_SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)), 'scripts')

ruleorder: polish_isolate_racon_ill > skip_illumina_polish

####################
# isolate assembly #
####################
# rule assemble_reads_redbean:
#     input:
#         reads = config["long_reads"]
#     output:
#         contigs = "isolate/redbean/{assembly}"

rule assemble_reads_flye:
    input:
        reads = config["long_reads"]
    output:
        contigs = "isolate/flye/assembly.fasta"
    threads:
        config["max_threads"]
    benchmark:
        "benchmarks/assemble_reads_flye.benchmark.txt"
    shell:
        f'{pixi_run} -e flye flye --nano-raw {{input.reads}} --threads {{threads}} -o isolate/flye'


rule polish_isolate_racon:
    input:
        fastq = config["long_reads"],
        fasta = "isolate/flye/assembly.fasta"
    threads:
        config["max_threads"]
    params:
        prefix = "second",
        maxcov = 1000,
        rounds = 4,
        illumina = False
    output:
        fasta = "isolate/isolate.pol.rac.fasta"
    benchmark:
        "benchmarks/polish_isolate_racon.benchmark.txt"
    shell:
        f'{pixi_run} -e polishing {ISOLATE_SCRIPTS_DIR}/'+\
        """polish.py \
        --input-fastq {input.fastq} \
        --reference {input.fasta} \
        --output-dir data/polishing \
        --output-prefix {params.prefix} \
        --output-fasta {output.fasta} \
        --rounds {params.rounds} \
        --long-read-type {config[long_read_type]} \
        --medaka-model {config[medaka_model]} \
        --illumina {params.illumina} \
        --max-cov {params.maxcov} \
        --threads {threads} \
        --coassemble False \
        --log logs/polish_isolate_racon.log
        """


rule polish_isolate_medaka:
    input:
        reads = config["long_reads"],
        contigs = "isolate/isolate.pol.rac.fasta"
    threads:
        config["max_threads"]
    params:
        model = config["guppy_model"]
    output:
        fasta = "isolate/isolate.pol.med.fasta"
    benchmark:
        "benchmarks/polish_isolate_medaka.benchmark.txt"
    shell:
        f'{pixi_run} -e polishing medaka_consensus '
        '-i {input.reads} -d {input.contigs} -o isolate/medaka/ -t {threads} -m {params.model} && '
        'cp isolate/medaka/consensus.fasta {output.fasta}'


rule polish_isolate_pilon:
    input:
        reads = "data/short_reads.fastq.gz",
        fasta = "isolate/isolate.pol.med.fasta"
    output:
        fasta = "isolate/isolate.pol.pil.fasta",
    threads:
        config["max_threads"]
    benchmark:
        "benchmarks/polish_isolate_pilon.benchmark.txt"
    shell:
        f'{pixi_run} -e pilon bash -c "'
        'minimap2 -ax sr -t {threads} {input.fasta} {input.reads} | samtools view -b | '
        'samtools sort -o isolate/pilon.sort.bam - && '
        'samtools index isolate/pilon.sort.bam && '
        'pilon -Xmx64000m --genome {input.fasta} --frags isolate/pilon.sort.bam --threads {threads} '
        '--output isolate/isolate.pol.pil --fix bases"'


rule polish_isolate_racon_ill:
    input:
        fastq = "data/short_reads.fastq.gz",
        fasta = "isolate/isolate.pol.pil.fasta"
    output:
        fasta = "isolate/isolate.pol.fin.fasta"
    threads:
        config["max_threads"]
    params:
        prefix = "racon_ill",
        maxcov = 1000,
        rounds = 1,
        illumina = True
    benchmark:
        "benchmarks/polish_isolate_racon_ill.benchmark.txt"
    shell:
        f'{pixi_run} -e polishing {ISOLATE_SCRIPTS_DIR}/'+\
        """polish.py \
        --short-reads-1 {input.fastq} \
        --short-reads-2 none \
        --reference {input.fasta} \
        --output-dir data/polishing_ill \
        --output-prefix {params.prefix} \
        --output-fasta {output.fasta} \
        --rounds {params.rounds} \
        --long-read-type {config[long_read_type]} \
        --medaka-model {config[medaka_model]} \
        --illumina {params.illumina} \
        --max-cov {params.maxcov} \
        --threads {threads} \
        --coassemble False \
        --log logs/polish_isolate_racon_ill.log
        """


rule skip_illumina_polish:
    input:
        fasta = "isolate/isolate.pol.med.fasta"
    output:
        fasta = "isolate/isolate.pol.fin.fasta"
    shell:
        "cp {input.fasta} {output.fasta}"


rule dnaapler:
    input:
        fasta = "isolate/isolate.pol.fin.fasta"
    output:
        fasta = "isolate/completed_assembly.fasta"
    threads:
        config["max_threads"]
    benchmark:
        "benchmarks/dnaapler.benchmark.txt"
    shell:
        f'{pixi_run} -e dnaapler dnaapler all '
        '-i {input.fasta} -o isolate/dnaapler -t {threads} && '
        'cp isolate/dnaapler/dnaapler_reoriented.fasta {output.fasta}'
