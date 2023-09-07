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
    conda:
        "envs/flye.yaml"
    threads:
        config["max_threads"]
    benchmark:
        "benchmarks/assemble_reads_flye.benchmark.txt"
    shell:
        "flye --nano-raw {input.reads} --threads {threads} -o isolate/flye"


rule polish_isolate_racon:
    input:
        fastq = config["long_reads"],
        fasta = "isolate/flye/assembly.fasta"
    conda:
        "envs/racon.yaml"
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
    script:
        "scripts/racon_polish.py"


rule polish_isolate_medaka:
    input:
        reads = config["long_reads"],
        contigs = "isolate/isolate.pol.rac.fasta"
    conda:
        "envs/medaka.yaml"
    threads:
        config["max_threads"]
    params:
        model = config["guppy_model"]
    output:
        fasta = "isolate/isolate.pol.med.fasta"
    benchmark:
        "benchmarks/polish_isolate_medaka.benchmark.txt"
    shell:
        """
        medaka_consensus -i {input.reads} -d {input.contigs} -o isolate/medaka/ -t {threads} -m {params.model} && \
        cp isolate/medaka/consensus.fasta {output.fasta}
        """


rule polish_isolate_pilon:
    input:
        reads = "data/short_reads.fastq.gz",
        fasta = "isolate/isolate.pol.med.fasta"
    output:
        fasta = "isolate/isolate.pol.pil.fasta",
    threads:
        config["max_threads"]
    conda:
        "envs/pilon.yaml"
    benchmark:
        "benchmarks/polish_isolate_pilon.benchmark.txt"
    shell:
        """
        minimap2 -ax sr -t {threads} {input.fasta} {input.reads} | samtools view -b | 
        samtools sort -o isolate/pilon.sort.bam - && \
        samtools index isolate/pilon.sort.bam && \
        pilon -Xmx64000m --genome {input.fasta} --frags isolate/pilon.sort.bam --threads {threads} \
        --output isolate/isolate.pol.pil --fix bases
        """


rule polish_isolate_racon_ill:
    input:
        fastq = "data/short_reads.fastq.gz",
        fasta = "isolate/isolate.pol.pil.fasta"
    output:
        fasta = "isolate/isolate.pol.fin.fasta"
    threads:
        config["max_threads"]
    conda:
        "envs/racon.yaml"
    params:
        prefix = "racon_ill",
        maxcov = 1000,
        rounds = 1,
        illumina = True
    benchmark:
        "benchmarks/polish_isolate_racon_ill.benchmark.txt"
    script:
        "scripts/racon_polish.py"


rule skip_illumina_polish:
    input:
        fasta = "isolate/isolate.pol.med.fasta"
    output:
        fasta = "isolate/isolate.pol.fin.fasta"
    shell:
        "cp {input.fasta} {output.fasta}"


rule circlator:
    input:
        fasta = "isolate/isolate.pol.fin.fasta",
        reads = config["long_reads"]
    output:
        fasta = "isolate/completed_assembly.fasta"
    threads:
        config["max_threads"]
    conda:
        "envs/circlator.yaml"
    benchmark:
        "benchmarks/circlator.benchmark.txt"
    shell:
        """
        circlator all --threads {threads} {input.fasta} {input.reads} isolate/circlator && \
        cp isolate/circlator/06.fixstart.fasta {output.fasta}
        """
