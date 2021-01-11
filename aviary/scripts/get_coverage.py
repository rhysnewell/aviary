import subprocess
import os
import sys

if snakemake.config["long_reads"] != "none":
    if snakemake.config["long_read_type"] == "nanopore":
        subprocess.Popen("coverm contig -t %d -r %s --single %s -p minimap2-ont -m length trimmed_mean variance --bam-file-cache-directory data/binning_bams/ > data/long_cov.tsv" %
                         (snakemake.threads, snakemake.input.fasta, " ".join(snakemake.config["long_reads"])), shell=True).wait()
    elif snakemake.config["long_read_type"] == "pacbio":
        subprocess.Popen("coverm contig -t %d -r %s --single %s -p minimap2-pb -m length trimmed_mean variance --bam-file-cache-directory data/binning_bams/ > data/long_cov.tsv" %
                         (snakemake.threads, snakemake.input.fasta, " ".join(snakemake.config["long_reads"])), shell=True).wait()
    else:
        subprocess.Popen(
            "coverm contig -t %d -r %s --single %s -p minimap2-ont -m length trimmed_mean variance --bam-file-cache-directory data/binning_bams/ > data/long_cov.tsv" %
            (snakemake.threads, snakemake.input.fasta, " ".join(snakemake.config["long_reads"])), shell=True).wait()

if snakemake.config['short_reads_2'] != 'none':
    subprocess.Popen("coverm contig -t %d -r %s -1 %s -2 %s -m metabat --bam-file-cache-directory data/binning_bams/  > data/short_cov.tsv" %
                     (snakemake.threads, snakemake.input.fasta, " ".join(snakemake.config["short_reads_1"]), " ".join(snakemake.config["short_reads_2"])), shell=True).wait()

elif snakemake.config['short_reads_1']  != 'none':
    subprocess.Popen("coverm contig -t %d -r %s --interleaved %s -m metabat --bam-file-cache-directory data/binning_bams/  > data/short_cov.tsv" %
                     (snakemake.threads, snakemake.input.fasta, " ".join(snakemake.config["short_reads_1"])), shell=True).wait()

# subprocess.Popen("ls data/binning_bams/*.bam | parallel -j1 samtools index -@ %d {} {}.bai" % (snakemake.threads-1), shell=True).wait()
# Concatenate the two coverage files if both long and short exist
if snakemake.config["long_reads"] != "none" and (snakemake.config["short_reads_1"] != "none"):
    short_count = len(snakemake.config["short_reads_1"])
    long_count = len(snakemake.config["long_reads"])
    with open('data/coverm.cov', 'w') as file3:
        with open('data/short_cov.tsv', 'r') as file1:
            with open('data/long_cov.tsv', 'r') as file2:
                for idx, line1, line2 in enumerate(zip(file1, file2)):
                    long_values = line2.strip().split()[2::]
                    short_values = line1.strip().split()
                    if idx != 0:
                        long_cov = sum([float(x) for x in long_values[0::2]])
                        short_cov = sum([float(x) for x in short_values[3::2]])
                        tot_depth = (long_cov + short_cov) / (short_count + long_count)
                        line = (
                            "{contig_info}\t"
                            "{total_avg_depth}\t"
                            "{short}\t"
                            "{long}"
                        ).format(
                            contig_info = short_values[0:2],
                            total_avg_depth = tot_depth,
                            short = "\t".join(short_values[3::]),
                            long = "\t".join(long_values)
                        )
                        print(line, file=file3)
                    else:
                        line = (
                            "{short}\t{long}"
                        ).format(
                            short = short_values,
                            long = long_values
                        )
                        print(line, file=file3)
                        
elif snakemake.config["long_reads"] != "none":
    long_count = len(snakemake.config["long_reads"])
    with open('data/coverm.cov', 'w') as file3:
        with open('data/long_cov.tsv', 'r') as file1:
                for idx, line1, line2 in enumerate(zip(file1, file2)):
                    long_values = line2.strip().split()
                    if idx != 0:
                        long_cov = sum([float(x) for x in long_values[2::2]])
                        tot_depth = long_cov / long_count
                        line = (
                            "{contig_info}\t"
                            "{total_avg_depth}\t"
                            "{long}"
                        ).format(
                            contig_info = long_values[0:2],
                            total_avg_depth = tot_depth,
                            long = "\t".join(long_values)
                        )
                        print(line, file=file3)
                    else:
                        line = (
                            "{contig_info}\ttotalAvgDepth\t{samples}"
                        ).format(
                            contig_info = long_values[0:2],
                            samples = long_values[2::]
                        )
                        print(line, file=file3)
                        
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
# with open("data/binning_bams/done", 'w') as o:
#     o.write('done')




