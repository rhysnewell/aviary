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
    # Path.joinpath(Path(output_dir), "assembly_info.txt").touch()
    # Path.joinpath(Path(output_dir), "assembly_graph.gfa").touch()
    # os.makedirs(f"{output_dir}/00-assembly", exist_ok=True)
    # os.makedirs(f"{output_dir}/10-consensus", exist_ok=True)
    # os.makedirs(f"{output_dir}/20-repeat", exist_ok=True)
    # os.makedirs(f"{output_dir}/30-contigger", exist_ok=True)
    # os.makedirs(f"{output_dir}/40-polishing", exist_ok=True)


def run_metamdbg(
    long_read_type: str,
    input_fastq: str,
    output_dir: str,
    threads: int,
    log: str,
):
    # check if input_fastq has any reads
    if os.path.getsize(input_fastq) == 0:
        with open(log, "a") as logf:
            logf.write(f"Input fastq file {input_fastq} is empty\n")
            logf.write(f"Skipping metamdbg assembly\n")
        create_dummy_output(output_dir)
        return

    if long_read_type == 'ont':
        read_type = "--in-ont"
    elif long_read_type == 'ont_hq':
        read_type = "--in-ont"
    elif long_read_type == 'ccs' or long_read_type == 'hifi':
        read_type = "--in-hifi"
    else:
        raise ValueError(f"Invalid long read type: {long_read_type}")
    
    cmd = f"metaMDBG asm {read_type} {input_fastq} --out-dir {output_dir} --threads {threads}".split()
    with open(log, "a") as logf:
        run(cmd, stdout=logf, stderr=STDOUT)


if __name__ == '__main__':
    long_read_type = snakemake.params.long_read_type
    input_fastq = snakemake.input.fastq
    output_dir = "data/metamdbg"
    threads = snakemake.threads
    log = snakemake.log[0]

    with open(log, "w") as logf: pass

    run_metamdbg(
        long_read_type=long_read_type,
        input_fastq=input_fastq,
        output_dir=output_dir,
        threads=threads,
        log=log,
    )
