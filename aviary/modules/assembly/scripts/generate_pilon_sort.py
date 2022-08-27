import subprocess
import os
import logging

logging.info("Generating BAM files for pilon...")
if os.path.exists('data/short_reads.fastq.gz'):
    subprocess.Popen("minimap2 -ax sr -t %d %s data/short_reads.fastq.gz | samtools view -b -F 4 -@ %d | "
                     "samtools sort -@ %d -o %s -; samtools index -@ %d %s; " %
                     (snakemake.threads, snakemake.input.fasta,
                      snakemake.threads, snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.bam,
                      ), shell=True).wait()

elif snakemake.config['short_reads_2'] != 'none':
    if len(snakemake.config['short_reads_2']) == 1 or not snakemake.params.coassemble:
        pe1 = snakemake.config['short_reads_1'][0]
        pe2 = snakemake.config['short_reads_2'][0]
    else:
        if not os.path.exists("data/short_reads.1.fastq.gz"):
            for reads1, reads2 in zip(snakemake.config['short_reads_1'], snakemake.config['short_reads_2']):
                subprocess.Popen(f"cat {reads1} >> data/short_reads.1.fastq.gz", shell=True).wait()
                subprocess.Popen(f"cat {reads2} >> data/short_reads.2.fastq.gz", shell=True).wait()
        pe1 = "data/short_reads.1.fastq.gz"
        pe2 = "data/short_reads.2.fastq.gz"

    subprocess.Popen("minimap2 -ax sr -t %d %s %s %s | samtools view -b -F 4 -@ %d | "
                     "samtools sort -@ %d -o %s -; samtools index -@ %d %s; " %
                     (snakemake.threads, snakemake.input.fasta, pe1, pe2,
                      snakemake.threads, snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.bam,
                      ), shell=True).wait()

    # os.remove("data/short_reads.1.fastq.gz")
    # os.remove("data/short_reads.2.fastq.gz")

elif snakemake.config['short_reads_1']  != 'none':
    if len(snakemake.config['short_reads_1']) == 1:
        pe1 = snakemake.config['short_reads_1'][0]
    else:
        if not os.path.exists("data/short_reads.1.fastq.gz") or not snakemake.params.coassemble:
            for reads1 in snakemake.config['short_reads_1']:
                subprocess.Popen(f"cat {reads1} >> data/short_reads.1.fastq.gz", shell=True).wait()
        pe1 = "data/short_reads.1.fastq.gz"

    subprocess.Popen("minimap2 -ax sr -t %d %s %s | samtools view -b -F 4 -@ %d | "
                     "samtools sort -@ %d -o %s -; samtools index -@ %d %s; " %
                     (snakemake.threads, snakemake.input.fasta, pe1,
                      snakemake.threads, snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.bam,
                      ), shell=True).wait()

    # os.remove("data/short_reads.1.fastq.gz")
