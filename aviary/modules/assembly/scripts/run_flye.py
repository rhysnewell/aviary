import os
import subprocess

if snakemake.params.long_read_type == 'ont':
    subprocess.Popen(
        "flye --nano-raw %s --meta -o data/flye -t %d" %
        (snakemake.input.fastq, snakemake.threads), shell=True).wait()
elif snakemake.params.long_read_type == 'ont_hq':
    subprocess.Popen(
        "flye --nano-hq %s --meta -o data/flye -t %d" %
        (snakemake.input.fastq, snakemake.threads), shell=True).wait()
elif snakemake.params.long_read_type == 'ccs':
    subprocess.Popen(
        "flye --pacbio-hifi %s --meta -o data/flye -t %d" %
        (snakemake.input.fastq, snakemake.threads), shell=True).wait()
else:
    subprocess.Popen(
        "flye --pacbio-raw %s --meta -o data/flye -t %d" %
        (snakemake.input.fastq, snakemake.threads), shell=True).wait()