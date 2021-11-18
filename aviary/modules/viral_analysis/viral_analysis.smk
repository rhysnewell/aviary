rule rename_contigs:
    input:
        fasta = config["fasta"]
    output:
        "data/renamed_contigs.fasta"
    shell:
        "sed -i 's/>/>${{input.fasta}%%_*}_/' {input.fasta}"


rule run_virsorter:
    input:
        fasta = "data/renamed_contigs.fasta",
        virsorter_data = config["virsorter_data"]
    output:
        "data/virsorter/done"
    conda:
        "envs/virsorter.yaml"
    threads:
        config["max_threads"]
    shell:
        "virsorter -f {input.fasta} --wdir data/virsorter --data-dir {input.virsorter_data} --ncpu {threads} &&" \
        "touch data/virsorter/done"