import subprocess

if snakemake.config['short_reads_2'] != 'none' and snakemake.config['short_reads_1'] != 'none':
    if len(snakemake.config['short_reads_2']) == 1:
        subprocess.Popen(
            "spades.py --memory %s --meta -t %d -o data/short_read_assembly -1 %s -2 %s -k %s" %
            (snakemake.config["max_memory"], snakemake.threads, snakemake.config["short_reads_1"][0],
             snakemake.config["short_reads_2"][0], " ".join(snakemake.params.kmer_sizes)), shell=True).wait()
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
            "mkdir -p data/short_read_assembly; cp data/megahit_assembly/final.contigs.fa data/short_read_assembly/scaffolds.fasta", shell=True
        ).wait()

elif snakemake.config['short_reads_1']  != 'none':
    if len(snakemake.config["short_reads_1"]) == 1:
        subprocess.Popen(
            "spades.py --memory %s --meta -t %d -o data/short_read_assembly --12 %s -k %s" %
            (snakemake.config["max_memory"], snakemake.threads, snakemake.config["short_reads_1"][0], " ".join(snakemake.params.kmer_sizes)),
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
            "mkdir -p data/short_read_assembly; cp data/megahit_assembly/final.contigs.fa data/short_read_assembly/scaffolds.fasta", shell=True
        ).wait()


# subprocess.Popen(
#     "ln -s data/spades_assembly/scaffolds.fasta data/final_contigs.fasta", shell=True
# ).wait()
