import subprocess
import os

if os.path.exists('data/short_reads.fastq.gz'):
    subprocess.Popen("minimap2 -ax sr -t %d %s data/short_reads.fastq.gz | samtools view -b -@ %d | "
                     "samtools sort -@ %d -o %s -; samtools index -@ %d %s; "
                     "samtools bam2fq -@ %d -f 12 %s | pigz -p %d > %s" %
                     (snakemake.threads, snakemake.input.fasta,
                      snakemake.threads, snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.fastq,
                      ), shell=True).wait()
elif snakemake.config['short_reads_2'] != 'none':
    subprocess.Popen("minimap2 -ax sr -t %d %s %s %s | samtools view -b -@ %d | "
                     "samtools sort -@ %d -o %s -; samtools index -@ %d %s; "
                     "samtools bam2fq -@ %d -f 12 %s | pigz -p %d > %s" %
                     (snakemake.threads, snakemake.input.fasta,
                      " ".join(snakemake.config["short_reads_1"]),
                      " ".join(snakemake.config["short_reads_2"]),
                      snakemake.threads, snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.fastq,
                      ), shell=True).wait()

elif snakemake.config['short_reads_1']  != 'none':
    subprocess.Popen("minimap2 -ax sr -t %d %s %s | samtools view -b -@ %d | "
                     "samtools sort -@ %d -o %s -; "
                     "samtools index -@ %d %s; "
                     "samtools bam2fq -@ %d -f 12 %s | "
                     "pigz -p %d > %s" %
                     (snakemake.threads, snakemake.input.fasta,
                      " ".join(snakemake.config["short_reads_1"]),
                      snakemake.threads, snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.fastq,
                      ), shell=True).wait()