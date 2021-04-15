import subprocess
import os
import logging


if not os.path.exists(snakemake.output.fasta):
    try:
        subprocess.Popen("cat %s %s > %s" % (snakemake.input.flye_fasta,
                                             snakemake.input.unicyc_fasta,
                                             snakemake.output.fasta))
        logging.info("Flye and Unicycler assemblies combined...")
    except AttributeError:
        subprocess.Popen("ln %s %s" % (snakemake.input.fasta,
                                       snakemake.output.fasta))
        logging.info("Treating Flye assembly as final assembly...")


#
# if os.path.exists(snakemake.output.fasta) \
#         and not os.path.exists(snakemake.output.long_bam):
#     logging.info("Mapping long reads to final assembly...")
#     subprocess.Popen("minimap2 -t %d -ax %s -a %s %s | "
#                      "samtools view -@ %d -b |"
#                      "samtools sort -@ %d -o %s -; "
#                      "samtools index -@ %d %s" %
#                      (snakemake.threads, 'map-ont' if snakemake.config["long_read_type"] == 'ont' else 'map-pb', snakemake.output.fasta, snakemake.input.long_reads,
#                       snakemake.threads,
#                       snakemake.threads, snakemake.output.long_bam,
#                       snakemake.threads, snakemake.output.long_bam))
#
# try:
#     logging.info("Mapping short reads to final assembly...")
#     if not os.path.exists(snakemake.output.short_bam):
#         if os.path.exists('data/short_reads.fastq.gz'):
#             subprocess.Popen("minimap2 -ax sr -t %d %s data/short_reads.fastq.gz | samtools view -b -@ %d | "
#                              "samtools sort -@ %d -o %s -; samtools index -@ %d %s; " %
#                              (snakemake.threads, snakemake.output.fasta,
#                               snakemake.threads, snakemake.threads, snakemake.output.short_bam,
#                               snakemake.threads, snakemake.output.short_bam,
#                               ), shell=True).wait()
#         if snakemake.config['short_reads_2'] != 'none':
#             subprocess.Popen("minimap2 -ax sr -t %d %s %s %s | samtools view -b -F 4 -@ %d | "
#                              "samtools sort -@ %d -o %s -; samtools index -@ %d %s; " %
#                              (snakemake.threads, snakemake.output.fasta,
#                               " ".join(snakemake.config["short_reads_1"]),
#                               " ".join(snakemake.config["short_reads_2"]),
#                               snakemake.threads, snakemake.threads, snakemake.output.short_bam,
#                               snakemake.threads, snakemake.output.short_bam,
#                               ), shell=True).wait()
#
#         elif snakemake.config['short_reads_1'] != 'none':
#             subprocess.Popen("minimap2 -ax sr -t %d %s %s | samtools view -b -@ %d | "
#                              "samtools sort -@ %d -o %s -; samtools index -@ %d %s; " %
#                              (snakemake.threads, snakemake.output.fasta,
#                               " ".join(snakemake.config["short_reads_1"]),
#                               snakemake.threads, snakemake.threads, snakemake.output.short_bam,
#                               snakemake.threads, snakemake.output.short_bam,
#                               ), shell=True).wait()
# except AttributeError:
#     logging.info("No short reads to map...")
#     pass
