import sys
import os

os.makedirs('www/fraction_recovered', exist_ok=True)

if snakemake.config["long_reads"] != "none" and not os.path.exists("data/long_cov.tsv"):
    if snakemake.config["long_read_type"][0] in ["ont", "ont_hq"]:
        subprocess.Popen("coverm genome -t %d -r %s --single %s -p minimap2-ont --min-read-percent-identity 0.85 > www/fraction_recovered/long_fraction_recovered.tsv" %
                         (snakemake.threads, snakemake.input.fasta, " ".join(snakemake.config["long_reads"])), shell=True).wait()
    elif snakemake.config["long_read_type"][0] in ["rs", "sq", "ccs"]:
        subprocess.Popen("coverm genome -t %d -r %s --single %s -p minimap2-pb --min-read-percent-identity 0.85 > www/fraction_recovered/long_fraction_recovered.tsv" %
                         (snakemake.threads, snakemake.input.fasta, " ".join(snakemake.config["long_reads"])), shell=True).wait()
    else:
        subprocess.Popen(
            "coverm genome -t %d -r %s --single %s -p minimap2-ont > www/fraction_recovered/long_fraction_recovered.tsv" %
            (snakemake.threads, snakemake.input.fasta, " ".join(snakemake.config["long_reads"])), shell=True).wait()

if snakemake.config['short_reads_2'] != 'none' and not os.path.exists("data/short_cov.tsv"):
    subprocess.Popen("coverm genome -t %d -r %s -1 %s -2 %s > www/fraction_recovered/short_fraction_recovered.tsv" %
                     (snakemake.threads, snakemake.input.fasta, " ".join(snakemake.config["short_reads_1"]), " ".join(snakemake.config["short_reads_2"])), shell=True).wait()

elif snakemake.config['short_reads_1']  != 'none' and not os.path.exists("data/short_cov.tsv"):
    subprocess.Popen("coverm genome -t %d -r %s --interleaved %s > www/fraction_recovered/short_fraction_recovered.tsv" %
                     (snakemake.threads, snakemake.input.fasta, " ".join(snakemake.config["short_reads_1"])), shell=True).wait()