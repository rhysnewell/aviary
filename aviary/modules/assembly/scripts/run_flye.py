#!/usr/bin/env python3

from subprocess import run, STDOUT
import os
from pathlib import Path
import argparse

def dummy_fasta(output_fasta: str):
    with open(output_fasta, "w") as out:
        out.write(">dummy\n")
        out.write("ACGT\n")

def create_dummy_output(output_dir: str):
    # need to create dummy output files
    os.makedirs(output_dir, exist_ok=True)
    output_fasta = f"{output_dir}/assembly.fasta"
    dummy_fasta(output_fasta)
    Path.joinpath(Path(output_dir), "assembly_info.txt").touch()
    Path.joinpath(Path(output_dir), "assembly_graph.gfa").touch()
    os.makedirs(f"{output_dir}/00-assembly", exist_ok=True)
    os.makedirs(f"{output_dir}/10-consensus", exist_ok=True)
    os.makedirs(f"{output_dir}/20-repeat", exist_ok=True)
    os.makedirs(f"{output_dir}/30-contigger", exist_ok=True)
    os.makedirs(f"{output_dir}/40-polishing", exist_ok=True)

def run_flye(
    long_read_type: str,
    input_fastq: str,
    output_dir: str,
    meta_flag: bool,
    threads: int,
    log: str,
):
    # check if input_fastq has any reads
    if os.path.getsize(input_fastq) == 0:
        with open(log, "a") as logf:
            logf.write(f"Input fastq file {input_fastq} is empty\n")
            logf.write(f"Skipping flye assembly\n")
        create_dummy_output(output_dir)
        return

    meta = ""
    flye_type = "--nano-raw"
    if meta_flag:
        meta = "--meta"
    if long_read_type == 'ont':
        flye_type = "--nano-raw"
    elif long_read_type == 'ont_hq':
        flye_type = "--nano-hq"
    elif long_read_type == 'ccs' or long_read_type == 'hifi':
        flye_type = "--pacbio-hifi"
    else:
        flye_type = "--pacbio-raw"
    
    flye_cmd = f"flye {flye_type} {input_fastq} {meta} -o {output_dir} -t {threads}".split()
    with open(log, "a") as logf:
        run(flye_cmd, stdout=logf, stderr=STDOUT)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run Flye assembly.")
    parser.add_argument("--long-read-type", required=True, help="Type of long reads (e.g., ont, ont_hq, ccs, hifi).")
    parser.add_argument("--input-fastq", required=True, help="Path to the input FASTQ file.")
    parser.add_argument("--output-dir", default="data/flye", help="Directory for Flye output.")
    parser.add_argument("--meta-flag", action="store_true", help="Enable meta assembly mode.")
    parser.add_argument("--threads", type=int, required=True, help="Number of threads to use.")
    parser.add_argument("--log", required=True, help="Path to the log file.")

    args = parser.parse_args()

    # Clear the log file
    with open(args.log, "w") as logf:
        pass

    run_flye(
        long_read_type=args.long_read_type,
        input_fastq=args.input_fastq,
        output_dir=args.output_dir,
        meta_flag=args.meta_flag,
        threads=args.threads,
        log=args.log,
    )