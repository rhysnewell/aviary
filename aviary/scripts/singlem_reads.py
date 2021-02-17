import subprocess
import os
import sys

try:
    os.mkdir("data/singlem_out/")
except OSError:
    print("Using prexisting directory: data/singlem_out/")

if snakemake.config["long_reads"] != "none":
        subprocess.Popen(
            "singlem pipe --threads %d --sequences %s --otu_table data/singlem_out/metagenome.longread_otu_table.csv" %
            (snakemake.config["pplacer_threads"], " ".join(snakemake.config["long_reads"])), shell=True).wait()


if snakemake.config["short_reads_2"] != "none":
        subprocess.Popen(
            "singlem pipe --threads %d --forward %s --reverse %s --otu_table data/singlem_out/metagenome.shortread_otu_table.csv" %
            (snakemake.config["pplacer_threads"],
             " ".join(snakemake.config["short_reads_1"]),
             " ".join(snakemake.config["short_reads_2"])), shell=True).wait()
elif snakemake.config["short_reads_1"] != "none":
        subprocess.Popen(
            "singlem pipe --threads %d --sequences %s --otu_table data/singlem_out/metagenome.shortread_otu_table.csv" %
            (snakemake.config["pplacer_threads"],
             " ".join(snakemake.config["short_reads_1"])), shell=True).wait()

subprocess.Popen(
            "singlem summarise --input_otu_tables data/singlem_out/*.csv --output_otu_table data/singlem_out/metagenome.combined_otu_table.csv", shell=True).wait()