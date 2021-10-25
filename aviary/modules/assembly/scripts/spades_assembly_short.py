import subprocess

if snakemake.config['short_reads_2'] != 'none':
    if len(snakemake.config['short_reads_2']) == 1:
        subprocess.Popen(
            "spades.py --memory %s --meta -t %d -o data/spades_assembly -k 21,33,55,81,99,127 -1 %s -2 %s" %
            (snakemake.config["max_memory"], snakemake.threads, snakemake.config["short_reads_1"][0],
             snakemake.config["short_reads_2"][0]), shell=True).wait()
    else: # co-assembly so use megahit
        subprocess.Popen(
            "megahit -1 %s -2 %s -t %d -m %d -o data/megahit_assembly" %
            (
                ",".join(snakemake.config["short_reads_1"]),
                ",".join(snakemake.config["short_reads_2"]),
                snakemake.threads,
                snakemake.config["max_memory"]
             ),
            shell=True
        ).wait()
        subprocess.Popen(
            "mkdir -p data/spades_assembly; cp data/megahit_assembly/final.contigs.fa data/spades_assembly/scaffolds.fasta", shell=True
        ).wait()

elif snakemake.config['short_reads_1']  != 'none':
    if len(snakemake.config["short_reads_1"]) == 1:
        subprocess.Popen(
            "spades.py --memory %s --meta -t %d -o data/spades_assembly -k 21,33,55,81,99,127 -1 %s" %
            (snakemake.config["max_memory"], snakemake.threads, snakemake.config["short_reads_1"][0]),
            shell=True).wait()
    else: # co-assembly so use megahit
        subprocess.Popen(
            "megahit --12 %s -t %d -m %d -o data/megahit_assembly" %
            (
                ",".join(snakemake.config["short_reads_1"]),
                snakemake.threads,
                snakemake.config["max_memory"]
             ),
            shell=True
        ).wait()
        subprocess.Popen(
            "mkdir -p data/spades_assembly; cp data/megahit_assembly/final.contigs.fa data/spades_assembly/scaffolds.fasta", shell=True
        ).wait()


# subprocess.Popen(
#     "ln -s data/spades_assembly/scaffolds.fasta data/final_contigs.fasta", shell=True
# ).wait()
