import subprocess
import os
import sys


if snakemake.config["long_reads"] != "none":
    if snakemake.config["long_read_type"][0] == "ont":
        subprocess.Popen("coverm genome -t %d -d data/galah_bins/ --single %s -p minimap2-ont --min-covered-fraction 0.0 -x fa %s --discard-unmapped > data/long_abundances.tsv" %
                         (snakemake.threads, " ".join(snakemake.config["long_reads"]),
                          '--bam-file-cache-directory data/binned_bams/' if snakemake.config['strain_analysis'] is True else ''), shell=True).wait()
    elif snakemake.config["long_read_type"][0] in ["rs", "sq", "ccs"]:
        subprocess.Popen("coverm genome -t %d -d data/galah_bins/ --single %s -p minimap2-pb --min-covered-fraction 0.0 -x fa %s --discard-unmapped > data/long_abundances.tsv" %
                         (snakemake.threads, " ".join(snakemake.config["long_reads"]),
                          '--bam-file-cache-directory data/binned_bams/' if snakemake.config['strain_analysis'] is True else ''), shell=True).wait()
    else:
        subprocess.Popen("coverm genome -t %d -d data/galah_bins/ --single %s -p minimap2-ont --min-covered-fraction 0.0 -x fa %s --discard_unmapped > data/long_abundances.tsv" %
                         (snakemake.threads, " ".join(snakemake.config["long_reads"]),
                          '--bam-file-cache-directory data/binned_bams/' if snakemake.config['strain_analysis'] is True else ''), shell=True).wait()

if snakemake.config['short_reads_2'] != 'none':
    subprocess.Popen("coverm genome -t %d -d data/galah_bins/ -1 %s -2 %s --min-covered-fraction 0.0 -x fa %s --discard-unmapped > data/short_abundances.tsv" %
                     (snakemake.threads, " ".join(snakemake.config["short_reads_1"]), " ".join(snakemake.config["short_reads_2"]),
                      '--bam-file-cache-directory data/binned_bams/' if snakemake.config['strain_analysis'] is True else ''), shell=True).wait()

elif snakemake.config['short_reads_1'] != 'none':
    subprocess.Popen("coverm genome -t %d -d data/galah_bins/ --interleaved %s --min-covered-fraction 0.0 -x fa %s --discard-unmapped > data/short_abundances.tsv" %
                     (snakemake.threads, " ".join(snakemake.config["short_reads_1"]),
                      '--bam-file-cache-directory data/binned_bams/' if snakemake.config['strain_analysis'] is True else ''), shell=True).wait()

# Concatenate the two coverage files if both long and short exist
if snakemake.config["long_reads"] != "none" and (snakemake.config["short_reads_1"] != "none"):
    with open('data/coverm_abundances.tsv', 'w') as file3:
        with open('data/short_abundances.tsv', 'r') as file1:
            with open('data/long_abundances.tsv', 'r') as file2:
                for line1, line2 in zip(file1, file2):
                    print(line1.strip(), "\t".join(line2.strip().split()[1::]), file=file3)
elif snakemake.config["long_reads"] != "none":  # rename long reads cov if only it exists
    os.rename("data/long_abundances.tsv", "data/coverm_abundances.tsv")
elif snakemake.config["short_reads_1"] != "none":  # rename shrot reads cov if only they exist
    os.rename("data/short_abundances.tsv", "data/coverm_abundances.tsv")
