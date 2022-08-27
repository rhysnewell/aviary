import subprocess
from pathlib import Path
import os

if snakemake.config['short_reads_2'] != 'none' and snakemake.config['short_reads_1'] != 'none':
    if len(snakemake.config['short_reads_2']) == 1 or not snakemake.params.coassemble:
        subprocess.Popen(
            "spades.py --memory %s --meta -t %d -o data/short_read_assembly -1 %s -2 %s -k %s --tmp-dir %s" %
            (snakemake.config["max_memory"], snakemake.threads, snakemake.config["short_reads_1"][0],
             snakemake.config["short_reads_2"][0], " ".join(snakemake.params.kmer_sizes), snakemake.params.tmpdir), shell=True).wait()
    elif snakemake.params.use_megahit: # co-assembly so use megahit
        subprocess.Popen(
            "megahit -1 %s -2 %s -t %d -m %d -o data/megahit_assembly --tmp-dir %s" %
            (
                ",".join(snakemake.config["short_reads_1"]),
                ",".join(snakemake.config["short_reads_2"]),
                snakemake.threads,
                snakemake.config["max_memory"],
                snakemake.params.tmpdir
             ),
            shell=True
        ).wait()
        subprocess.Popen(
            "mkdir -p data/short_read_assembly; cp data/megahit_assembly/final.contigs.fa data/short_read_assembly/scaffolds.fasta", shell=True
        ).wait()
    else:
        for reads1, reads2 in zip(snakemake.config['short_reads_1'], snakemake.config['short_reads_2']):
            subprocess.Popen(f"cat {reads1} >> data/short_reads.1.fastq.gz", shell=True).wait()
            subprocess.Popen(f"cat {reads2} >> data/short_reads.2.fastq.gz", shell=True).wait()

        subprocess.Popen(
            "spades.py --memory %s --meta -t %d -o data/short_read_assembly "
            "-1 data/short_reads.1.fastq.gz -2 data/short_reads.2.fastq.gz -k %s --tmp-dir %s" %
            (snakemake.config["max_memory"], snakemake.threads, " ".join(snakemake.params.kmer_sizes), snakemake.params.tmpdir),
            shell=True).wait()

        os.remove("data/short_reads.1.fastq.gz")
        os.remove("data/short_reads.2.fastq.gz")


elif snakemake.config['short_reads_1']  != 'none':
    if len(snakemake.config["short_reads_1"]) == 1 or not snakemake.params.coassemble:
        subprocess.Popen(
            "spades.py --memory %s --meta -t %d -o data/short_read_assembly --12 %s -k %s --tmp-dir %s" %
            (snakemake.config["max_memory"], snakemake.threads, snakemake.config["short_reads_1"][0], " ".join(snakemake.params.kmer_sizes), snakemake.params.tmpdir),
            shell=True).wait()
    elif snakemake.params.use_megahit: # co-assembly so use megahit
        subprocess.Popen(
            "megahit --12 %s -t %d -m %d -o data/megahit_assembly --tmp-dir %s" %
            (
                ",".join(snakemake.config["short_reads_1"]),
                snakemake.threads,
                snakemake.config["max_memory"],
                snakemake.params.tmpdir
             ),
            shell=True
        ).wait()
        subprocess.Popen(
            "mkdir -p data/short_read_assembly; cp data/megahit_assembly/final.contigs.fa data/short_read_assembly/scaffolds.fasta", shell=True
        ).wait()
    else:
        for reads1 in snakemake.config['short_reads_1']:
            subprocess.Popen(f"cat {reads1} >> data/short_reads.1.fastq.gz", shell=True).wait()

        subprocess.Popen(
            "spades.py --memory %s --meta -t %d -o data/short_read_assembly "
            "--12 data/short_reads.1.fastq.gz -k %s --tmp-dir %s" %
            (snakemake.config["max_memory"], snakemake.threads, " ".join(snakemake.params.kmer_sizes),
             snakemake.params.tmpdir),
            shell=True).wait()

        os.remove("data/short_reads.1.fastq.gz")



# subprocess.Popen(
#     "ln -s data/spades_assembly/scaffolds.fasta data/final_contigs.fasta", shell=True
# ).wait()
