from subprocess import run, STDOUT
import os
from pathlib import Path

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
    long_read_type = snakemake.params.long_read_type
    input_fastq = snakemake.input.fastq
    output_dir = "data/flye"
    meta_flag = True
    threads = snakemake.threads
    log = snakemake.log[0]

    with open(log, "w") as logf: pass

    run_flye(
        long_read_type=long_read_type,
        input_fastq=input_fastq,
        output_dir=output_dir,
        meta_flag=meta_flag,
        threads=threads,
        log=log,
    )