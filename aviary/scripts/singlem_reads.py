from subprocess import CalledProcessError, run, STDOUT
import os
from pathlib import Path


def run_singlem(
    long_reads,
    short_reads_1,
    short_reads_2,
    pplacer_threads: int,
    log: str,
):
    try:
        os.mkdir("data/singlem_out/")
    except OSError:
        with open(log, "a") as logf:
            logf.write("Using prexisting directory: data/singlem_out/\n")

    singlem_output_list = []
    if long_reads != "none":
        singlem_pipe_cmd = f"singlem pipe --threads {pplacer_threads} --sequences {' '.join(long_reads)} --otu-table data/singlem_out/metagenome.longread_otu_table.csv".split()
        with open(log, "a") as logf:
            run(singlem_pipe_cmd, stdout=logf, stderr=STDOUT)

    if short_reads_2 != "none":
        singlem_pipe_cmd = f"singlem pipe --threads {pplacer_threads} --forward {' '.join(short_reads_1)} --reverse {' '.join(short_reads_2)} --otu-table data/singlem_out/metagenome.shortread_otu_table.csv".split()
        with open(log, "a") as logf:
            run(singlem_pipe_cmd, stdout=logf, stderr=STDOUT)

    elif short_reads_1 != "none":
        singlem_pipe_cmd = f"singlem pipe --threads {pplacer_threads} --sequences {' '.join(short_reads_1)} --otu-table data/singlem_out/metagenome.shortread_otu_table.csv".split()
        with open(log, "a") as logf:
            run(singlem_pipe_cmd, stdout=logf, stderr=STDOUT)


    # if file exists then add it to otu output list
    if os.path.exists("data/singlem_out/metagenome.longread_otu_table.csv"):
        singlem_output_list.append("data/singlem_out/metagenome.longread_otu_table.csv")
    
    if os.path.exists("data/singlem_out/metagenome.shortread_otu_table.csv"):
        singlem_output_list.append("data/singlem_out/metagenome.shortread_otu_table.csv")

    summarise_cmd = f"singlem summarise --input-otu-tables {' '.join(singlem_output_list)} --output-otu-table data/singlem_out/metagenome.combined_otu_table.csv".split()
    
    try:
        with open(log, "a") as logf:
            run(summarise_cmd, stdout=logf, stderr=STDOUT)
        Path("data/singlem_out/metagenome.combined_otu_table.csv").touch()
    except CalledProcessError as e:
        with open(log, "a") as logf:
            logf.write(e)
            logf.write("\nSingleM summarise failed. Exiting.\n")
        Path("data/singlem_out/metagenome.combined_otu_table.csv").touch()


if __name__ == '__main__':
    long_reads = snakemake.config['long_reads']
    short_reads_1 = snakemake.config['short_reads_1']
    short_reads_2 = snakemake.config['short_reads_2']
    pplacer_threads = snakemake.config["pplacer_threads"]
    log = snakemake.log[0]

    with open(log, "w") as logf: pass

    run_singlem(
        long_reads,
        short_reads_1,
        short_reads_2,
        pplacer_threads,
        log,
    )
