import subprocess
import shutil
import os
from pathlib import Path

def checkm(checkm2_db, bin_folder, bin_ext, refinery_max_iterations, output_folder, output_file, threads):
    if len([f for f in os.listdir(bin_folder) if f.endswith(bin_ext)]) == 0:
        print(f"No bins found in {bin_folder}")
        os.makedirs(output_folder)
        Path(output_file).touch()
    elif refinery_max_iterations == 0:
        print("Skipping pre-refinery CheckM2 rules")
        os.makedirs(output_folder)
        Path(output_file).touch()
    else:
        print(f"Using CheckM2 database {checkm2_db}/uniref100.KO.1.dmnd")
        os.environ["CHECKM2DB"] = f"{checkm2_db}/uniref100.KO.1.dmnd"
        subprocess.run(f"checkm2 predict -i {bin_folder}/ -x {bin_ext} -o {output_folder} -t {threads} --force".split(), env=os.environ)
        shutil.copy(f"{output_folder}/quality_report.tsv", output_file)


if __name__ == '__main__':
    checkm2_db = snakemake.params.checkm2_db_path
    bin_folder = snakemake.params.bin_folder
    bin_ext = snakemake.params.extension
    refinery_max_iterations = snakemake.params.refinery_max_iterations
    output_folder = snakemake.output.output_folder
    output_file = snakemake.output.output_file
    threads = snakemake.threads

    checkm(checkm2_db, bin_folder, bin_ext, refinery_max_iterations, output_folder, output_file, threads)
