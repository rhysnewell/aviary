#!/usr/bin/env python3
 
import subprocess
import shutil
import os
from pathlib import Path
import argparse

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
    parser = argparse.ArgumentParser(description="Run CheckM2 with specified parameters.")
    parser.add_argument("--checkm2-db", required=True, help="Path to the CheckM2 database.")
    parser.add_argument("--bin-folder", required=True, help="Path to the folder containing bins.")
    parser.add_argument("--bin-ext", required=True, help="Extension of bin files.")
    parser.add_argument("--refinery-max-iterations", type=int, required=True, help="Maximum iterations for refinery.")
    parser.add_argument("--output-folder", required=True, help="Path to the output folder.")
    parser.add_argument("--output-file", required=True, help="Path to the output file.")
    parser.add_argument("--threads", type=int, required=True, help="Number of threads to use.")
    parser.add_argument("--log", required=True, help="Path to the log file.")

    args = parser.parse_args()

    with open(args.log, "w") as logf:
        pass

    checkm(
        args.checkm2_db,
        args.bin_folder,
        args.bin_ext,
        args.refinery_max_iterations,
        args.output_folder,
        args.output_file,
        args.threads,
        args.log,
    )
