import subprocess
import os
import sys

try:
    os.makedirs('data/binning_bams')
except OSError:
    pass
already_done = set()


if snakemake.config["long_reads"] != "none":
    subprocess.Popen("coverm contig -t %d -r %s --single %s -p minimap2-ont -m metabat --bam-file-cache-directory data/binning_bams/ > data/long_cov.tsv" %
                     (snakemake.threads, snakemake.input.fasta, snakemake.config["long_reads"]), shell=True).wait()

if snakemake.config['short_reads_2'] != 'none':
    subprocess.Popen("coverm contig -t %d -r %s -1 %s -2 %s -m metabat --bam-file-cache-directory data/binning_bams/  > data/short_cov.tsv" %
                     (snakemake.threads, snakemake.input.fasta, snakemake.config["short_reads_1"], snakemake.config["short_reads_2"]), shell=True).wait()

elif snakemake.config['short_reads_1']  != 'none':
    subprocess.Popen("coverm contig -t %d -r %s --interleaved %s -m metabat --bam-file-cache-directory data/binning_bams/  > data/short_cov.tsv" %
                     (snakemake.threads, snakemake.input.fasta, snakemake.config["short_reads_1"]), shell=True).wait()

subprocess.Popen("ls data/binning_bams/*.bam | parallel -j1 samtools -@ %d index {} {}.bai" %
                     (snakemake.threads-1), shell=True).wait()
# Concatenate the two coverage files if both long and short exist
if snakemake.config["long_reads"] != "none" and (snakemake.config["short_reads_1"] != "none"):
    with open('data/coverm.cov', 'w') as file3:
        with open('data/short_cov.tsv', 'r') as file1:
            with open('data/long_cov.tsv', 'r') as file2:
                for line1, line2 in zip(file1, file2):
                    print(line1.strip(), "\t".join(line2.strip().split()[1::]), file=file3)
elif snakemake.config["long_reads"] != "none":  # rename long reads cov if only it exists
    os.rename("data/long_cov.tsv", "data/coverm.cov")
elif snakemake.config["short_reads_1"] != "none":  # rename shrot reads cov if only they exist
    os.rename("data/short_cov.tsv", "data/coverm.cov")

try:
    os.makedirs("data/maxbin_cov/")
except OSError:
    pass
with open("data/coverm.cov") as f, open("data/maxbin.cov.list", "w") as o:
    cov_list = []
    contig_list = []
    f.readline()
    for line in f:
        contig = line.split()[0]
        contig_list.append(contig)
        depths = line.split()[3:]
        cov_list.append([])
        for i in range(0, len(depths), 2):
            cov_list[-1].append(depths[i])
    for i in range(len(cov_list[0])):
        with open("data/maxbin_cov/p%d.cov" % i, 'w') as oo:
            for j in range(len(contig_list)):
                oo.write(contig_list[j] + '\t' + cov_list[j][i] + '\n')
        o.write("data/maxbin_cov/p%d.cov\n" % i)
with open("data/binning_bams/done", 'w') as o:
    o.write('done')




