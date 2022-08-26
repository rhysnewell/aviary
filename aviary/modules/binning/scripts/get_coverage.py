import subprocess
import os

if snakemake.config["long_reads"] != "none" and not os.path.exists("data/long_cov.tsv"):
    if snakemake.config["long_read_type"][0] in ["ont", "ont_hq"]:
        subprocess.Popen("TMPDIR=%s coverm contig -t %d -r %s --single %s -p minimap2-ont -m length trimmed_mean variance --bam-file-cache-directory data/binning_bams/ --discard-unmapped --min-read-percent-identity 0.85 > data/long_cov.tsv" %
                         (snakemake.params.tmpdir, snakemake.threads, snakemake.input.fasta, " ".join(snakemake.config["long_reads"])), shell=True).wait()
    elif snakemake.config["long_read_type"][0] in ["rs", "sq", "ccs"]:
        subprocess.Popen("TMPDIR=%s coverm contig -t %d -r %s --single %s -p minimap2-pb -m length trimmed_mean variance --bam-file-cache-directory data/binning_bams/ --discard-unmapped --min-read-percent-identity 0.9 > data/long_cov.tsv" %
                         (snakemake.params.tmpdir, snakemake.threads, snakemake.input.fasta, " ".join(snakemake.config["long_reads"])), shell=True).wait()
    else:
        subprocess.Popen(
            "TMPDIR=%s coverm contig -t %d -r %s --single %s -p minimap2-ont -m length trimmed_mean variance --bam-file-cache-directory data/binning_bams/ --discard-unmapped > data/long_cov.tsv" %
            (snakemake.params.tmpdir, snakemake.threads, snakemake.input.fasta, " ".join(snakemake.config["long_reads"])), shell=True).wait()

if snakemake.config['short_reads_2'] != 'none' and not os.path.exists("data/short_cov.tsv"):
    subprocess.Popen("TMPDIR=%s coverm contig -t %d -r %s -1 %s -2 %s -m metabat --bam-file-cache-directory data/binning_bams/ --discard-unmapped > data/short_cov.tsv" %
                     (snakemake.params.tmpdir, snakemake.threads, snakemake.input.fasta, " ".join(snakemake.config["short_reads_1"]), " ".join(snakemake.config["short_reads_2"])), shell=True).wait()

elif snakemake.config['short_reads_1']  != 'none' and not os.path.exists("data/short_cov.tsv"):
    subprocess.Popen("TMPDIR=%s coverm contig -t %d -r %s --interleaved %s -m metabat --bam-file-cache-directory data/binning_bams/ --discard-unmapped > data/short_cov.tsv" %
                     (snakemake.params.tmpdir, snakemake.threads, snakemake.input.fasta, " ".join(snakemake.config["short_reads_1"])), shell=True).wait()

# subprocess.Popen("ls data/binning_bams/*.bam | parallel -j1 samtools index -@ %d {} {}.bai" % (snakemake.threads-1), shell=True).wait()
# Concatenate the two coverage files if both long and short exist
if snakemake.config["long_reads"] != "none" and (snakemake.config["short_reads_1"] != "none"):
    short_count = len(snakemake.config["short_reads_1"])
    long_count = len(snakemake.config["long_reads"])
    with open('data/coverm.cov', 'w') as file3:
        with open('data/short_cov.tsv', 'r') as file1:
            with open('data/long_cov.tsv', 'r') as file2:
                for idx, (line1, line2) in enumerate(zip(file1, file2)):
                    long_values = line2.strip().split("\t")[1::]
                    del long_values[0::3] # delete length value
                    short_values = line1.strip().split("\t")
                    if idx != 0:
                        long_cov = sum([float(x) for x in long_values[0::2]])
                        short_cov = sum([float(x) for x in short_values[3::2]])
                        tot_depth = (long_cov + short_cov) / (short_count + long_count)
                        line = (
                            "{contig_name}\t"
                            "{contig_length}\t"
                            "{total_avg_depth}\t"
                            "{short}\t"
                            "{long}"
                        ).format(
                            contig_name = short_values[0],
                            contig_length = short_values[1],
                            total_avg_depth = tot_depth,
                            short = "\t".join(short_values[3::]),
                            long = "\t".join(long_values)
                        )
                        print(line, file=file3)
                    else:
                        line = (
                            "{short}\t{long}"
                        ).format(
                            short = "\t".join(short_values),
                            long = "\t".join(long_values)
                        )
                        print(line, file=file3)
                        
elif snakemake.config["long_reads"] != "none":
    with open('data/coverm.cov', 'w') as file3:
        with open('data/long_cov.tsv', 'r') as file1:
            for idx, line1 in enumerate(file1):
                long_values = line1.strip().split("\t")
                del long_values[4::3] # delete extra length values
                long_count = len(long_values[2::2])
                if idx != 0:
                    long_cov = sum([float(x) for x in long_values[2::2]])
                    tot_depth = long_cov / long_count
                    line = (
                        "{contig_name}\t"
                        "{contig_length}\t"
                        "{total_avg_depth}\t"
                        "{long}"
                    ).format(
                        contig_name = long_values[0],
                        contig_length = long_values[1],
                        total_avg_depth = tot_depth,
                        long = "\t".join(long_values[2::])
                    )
                    print(line, file=file3)
                else:
                    line = (
                        "contigName\tcontigLen\ttotalAvgDepth\t{samples}"
                    ).format(
                        samples = "\t".join(long_values[2::])
                    )
                    print(line, file=file3)
                        
elif snakemake.config["short_reads_1"] != "none":  # rename shrot reads cov if only they exist
    os.rename("data/short_cov.tsv", "data/coverm.cov")

if snakemake.config["long_reads"] != "none":
    with open('data/long.cov', 'w') as file3:
        with open('data/long_cov.tsv', 'r') as file:
            for idx, line1 in enumerate(file):
                long_values = line1.strip().split("\t")
                del long_values[4::3] # delete extra length values
                long_count = len(long_values[2::2])
                if idx != 0:
                    long_cov = sum([float(x) for x in long_values[2::2]])
                    tot_depth = long_cov / long_count
                    line = (
                        "{contig_name}\t"
                        "{contig_length}\t"
                        "{total_avg_depth}\t"
                        "{long}"
                    ).format(
                        contig_name = long_values[0],
                        contig_length = long_values[1],
                        total_avg_depth = tot_depth,
                        long = "\t".join(long_values[2::])
                    )
                    print(line, file=file3)
                else:
                    line = (
                        "contigName\tcontigLen\ttotalAvgDepth\t{samples}"
                    ).format(
                        samples = "\t".join(long_values[2::])
                    )
                    print(line, file=file3)
        # os.remove("data/long_cov.tsv")

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




