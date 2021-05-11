import subprocess
from pathlib import Path

if snakemake.config['short_reads_2'] != 'none':
    subprocess.Popen("minimap2 -ax sr -t %d %s %s %s | samtools view -b -f 12 -@ %d > %s; "
                     "samtools index -@ %d %s; "
                     "samtools bam2fq -@ %d -f 12 %s | pigz -p %d > %s" %
                     (snakemake.threads, snakemake.input.reference_filter,
                      " ".join(snakemake.config["short_reads_1"]),
                      " ".join(snakemake.config["short_reads_2"]),
                      snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.fastq,
                      ), shell=True).wait()

elif snakemake.config['short_reads_1']  != 'none':
    subprocess.Popen("minimap2 -ax sr -t %d %s %s | samtools view -b -f 12 -@ %d > %s; "
                     "samtools index -@ %d %s; "
                     "samtools bam2fq -@ %d -f 12 %s | pigz -p %d > %s" %
                     (snakemake.threads, snakemake.input.reference_filter,
                      " ".join(snakemake.config["short_reads_1"]),
                      snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.bam,
                      snakemake.threads, snakemake.output.fastq,
                      ), shell=True).wait()

Path(snakemake.output.filtered).touch()