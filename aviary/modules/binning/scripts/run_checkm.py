import subprocess
import shutil
import os
from pathlib import Path

def checkm(checkm2_db, bin_folder, bin_ext, refinery_max_iterations, output_folder, output_file, threads, log):
    if len([f for f in os.listdir(bin_folder) if f.endswith(bin_ext)]) == 0:
        with open(log, "a") as logf:
            logf.write(f"No bins found in {bin_folder}\n")
        os.makedirs(output_folder)
        Path(output_file).touch()
    elif refinery_max_iterations == 0:
        with open(log, "a") as logf:
            logf.write("Skipping pre-refinery CheckM2 rules\n")
        os.makedirs(output_folder)
        Path(output_file).touch()
    else:
        os.environ["CHECKM2DB"] = f"{checkm2_db}/uniref100.KO.1.dmnd"
        with open(log, "a") as logf:
            logf.write(f"Using CheckM2 database {checkm2_db}/uniref100.KO.1.dmnd\n")
            subprocess.run(
                f"checkm2 predict -i {bin_folder}/ -x {bin_ext} -o {output_folder} -t {threads} --force".split(),
                env=os.environ,
                stdout=logf,
                stderr=subprocess.STDOUT
                )
        shutil.copy(f"{output_folder}/quality_report.tsv", output_file)


if __name__ == '__main__':
    checkm2_db = snakemake.params.checkm2_db_path
    bin_folder = snakemake.params.bin_folder
    bin_ext = snakemake.params.extension
    refinery_max_iterations = snakemake.params.refinery_max_iterations
    output_folder = snakemake.output.output_folder
    output_file = snakemake.output.output_file
    threads = snakemake.threads
    log = snakemake.log[0]

    with open(log, "w") as logf: pass

    checkm(
        checkm2_db,
        bin_folder,
        bin_ext,
        refinery_max_iterations,
        output_folder,
        output_file,
        threads,
        log,
    )
