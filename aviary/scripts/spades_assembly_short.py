import subprocess
import os

if snakemake.config['short_reads_2'] != 'none':

    if len(snakemake.config['short_reads_2']) > 1:
        subprocess.Popen(
            "spades.py --memory %s -t %d -o data/spades_assembly -k 21,33,55,81,99,127 %s %s" %
            (snakemake.config["max_memory"], snakemake.threads,
             [" ".join(['-pe-1 ' + str(tup[0] + 1), tup[1]]) for tup in enumerate(snakemake.config['short_reads_1'])],
             [" ".join(['-pe-2 ' + str(tup[0] + 1), tup[1]]) for tup in enumerate(snakemake.config['short_reads_2'])]),
            shell=True).wait()
    else:
        subprocess.Popen(
            "spades.py --memory %s --meta -t %d -o data/spades_assembly -k 21,33,55,81,99,127 -1 %s -2 %s" %
            (snakemake.config["max_memory"], snakemake.threads, " ".join(snakemake.config["short_reads_1"]),
             " ".join(snakemake.config["short_reads_2"])), shell=True).wait()

elif snakemake.config['short_reads_1']  != 'none':
    if len(snakemake.config['short_reads_2']) > 1:
        subprocess.Popen(
            "spades.py --memory %s -t %d -o data/spades_assembly -k 21,33,55,81,99,127 %s" %
            (snakemake.config["max_memory"], snakemake.threads,
             [" ".join(['-pe-12 ' + str(tup[0] + 1), tup[1]]) for tup in enumerate(snakemake.config['short_reads_1'])]),
            shell=True).wait()
    else:
        subprocess.Popen(
            "spades.py --memory %s --meta -t %d -o data/spades_assembly -k 21,33,55,81,99,127 -12 %s" %
            (snakemake.config["max_memory"], snakemake.threads, " ".join(snakemake.config["short_reads_1"])),
            shell=True).wait()

subprocess.Popen(
    "ln -s data/spades_assembly/scaffolds.fasta data/final_contigs.fasta", shell=True
).wait()
