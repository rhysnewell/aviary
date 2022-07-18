import subprocess
import os
import sys


if snakemake.config["long_reads"] != "none":
    if snakemake.config["long_read_type"][0] in ["ont", "ont_hq"]:
        subprocess.Popen("coverm genome -t %d -d bins/final_bins/ -m relative_abundance covered_fraction --single %s -p minimap2-ont --min-covered-fraction 0.0 -x fna %s > data/long_abundances.tsv" %
                         (snakemake.threads, " ".join(snakemake.config["long_reads"]),
                          '--bam-file-cache-directory data/reads_mapped_to_mags/long/ --discard-unmapped' if snakemake.config['strain_analysis'] is True else ''), shell=True).wait()
    elif snakemake.config["long_read_type"][0] in ["rs", "sq", "ccs"]:
        subprocess.Popen("coverm genome -t %d -d bins/final_bins/ -m relative_abundance covered_fraction --single %s -p minimap2-pb --min-covered-fraction 0.0 -x fna %s > data/long_abundances.tsv" %
                         (snakemake.threads, " ".join(snakemake.config["long_reads"]),
                          '--bam-file-cache-directory data/reads_mapped_to_mags/long/ --discard-unmapped' if snakemake.config['strain_analysis'] is True else ''), shell=True).wait()
    else:
        subprocess.Popen("coverm genome -t %d -d bins/final_bins/ -m relative_abundance covered_fraction --single %s -p minimap2-ont --min-covered-fraction 0.0 -x fna %s > data/long_abundances.tsv" %
                         (snakemake.threads, " ".join(snakemake.config["long_reads"]),
                          '--bam-file-cache-directory data/reads_mapped_to_mags/long/ --discard-unmapped' if snakemake.config['strain_analysis'] is True else ''), shell=True).wait()

if snakemake.config['short_reads_2'] != 'none':
    subprocess.Popen("coverm genome -t %d -d bins/final_bins/ -m relative_abundance covered_fraction -1 %s -2 %s --min-covered-fraction 0.0 -x fna %s > data/short_abundances.tsv" %
                     (snakemake.threads, " ".join(snakemake.config["short_reads_1"]), " ".join(snakemake.config["short_reads_2"]),
                      '--bam-file-cache-directory data/reads_mapped_to_mags/short/ --discard-unmapped' if snakemake.config['strain_analysis'] is True else ''), shell=True).wait()

elif snakemake.config['short_reads_1'] != 'none':
    subprocess.Popen("coverm genome -t %d -d bins/final_bins/ -m relative_abundance covered_fraction --interleaved %s --min-covered-fraction 0.0 -x fna %s > data/short_abundances.tsv" %
                     (snakemake.threads, " ".join(snakemake.config["short_reads_1"]),
                      '--bam-file-cache-directory data/reads_mapped_to_mags/short/ --discard-unmapped' if snakemake.config['strain_analysis'] is True else ''), shell=True).wait()

# Concatenate the two coverage files if both long and short exist
if snakemake.config["long_reads"] != "none" and (snakemake.config["short_reads_1"] != "none"):
    with open('data/coverm_abundances.tsv', 'w') as file3:
        with open('data/short_abundances.tsv', 'r') as file1:
            with open('data/long_abundances.tsv', 'r') as file2:
                for line1, line2 in zip(file1, file2):
                    long_cov_line = "\t".join([l.strip() for l in line2.strip().split('\t')[1::]])
                    print(line1.strip(), "\t", long_cov_line, file=file3)
elif snakemake.config["long_reads"] != "none":  # rename long reads cov if only it exists
    os.rename("data/long_abundances.tsv", "data/coverm_abundances.tsv")
elif snakemake.config["short_reads_1"] != "none":  # rename shrot reads cov if only they exist
    os.rename("data/short_abundances.tsv", "data/coverm_abundances.tsv")
