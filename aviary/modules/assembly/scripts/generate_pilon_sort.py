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
    for reads1, reads2 in zip(snakemake.config['short_reads_1'], snakemake.config['short_reads_2']):
        subprocess.Popen(f"cat {reads1} >> data/short_reads.1.fastq.gz", shell=True).wait()
        subprocess.Popen(f"cat {reads2} >> data/short_reads.2.fastq.gz", shell=True).wait()

    subprocess.Popen("minimap2 -ax sr -t %d %s data/short_reads.1.fastq.gz data/short_reads.2.fastq.gz | samtools view -b -F 4 -@ %d | "
                     "samtools sort -@ %d -o %s -; samtools index -@ %d %s; " %
                     (snakemake.threads, snakemake.input.fasta,
                      snakemake.threads, snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.bam,
                      ), shell=True).wait()

    os.remove("data/short_reads.1.fastq.gz")
    os.remove("data/short_reads.2.fastq.gz")

elif snakemake.config['short_reads_1']  != 'none':
    for reads1 in snakemake.config['short_reads_1']:
        subprocess.Popen(f"cat {reads1} >> data/short_reads.1.fastq.gz", shell=True).wait()

    subprocess.Popen("minimap2 -ax sr -t %d %s data/short_reads.1.fastq.gz | samtools view -b -F 4 -@ %d | "
                     "samtools sort -@ %d -o %s -; samtools index -@ %d %s; " %
                     (snakemake.threads, snakemake.input.fasta,
                      snakemake.threads, snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.bam,
                      ), shell=True).wait()

    os.remove("data/short_reads.1.fastq.gz")
