import argparse
import os
import subprocess
from pathlib import Path

from Bio import SeqIO


# Match Flye's assembly_info.txt header so downstream parsing remains consistent
HEADER = "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"


def write_assembly_info(fasta_path: str, info_path: str) -> None:
    Path(info_path).parent.mkdir(parents=True, exist_ok=True)
    with open(info_path, "w") as info_handle:
        info_handle.write(HEADER)
        for record in SeqIO.parse(fasta_path, "fasta"):
            info_handle.write(
                f"{record.id}\t{len(record.seq)}\t0\tN\tN\t1\tNA\tNA\n"
            )


def ensure_placeholder(path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).touch()


def run_myloasm(input_fastq: str, output_dir: str, long_read_type: str, threads: int, log: str) -> None:
    os.makedirs(output_dir, exist_ok=True)

    with open(log, "w") as logf:
        cmd = [
            "myloasm",
            "--reads",
            input_fastq,
            "--out-dir",
            output_dir,
            "--threads",
            str(threads),
        ]
        if long_read_type:
            cmd.extend(["--long-read-type", long_read_type])

        subprocess.run(cmd, check=True, stdout=logf, stderr=subprocess.STDOUT)

    assembly_fasta = os.path.join(output_dir, "assembly.fasta")
    assembly_graph = os.path.join(output_dir, "assembly_graph.gfa")
    assembly_info = os.path.join(output_dir, "assembly_info.txt")

    if not os.path.exists(assembly_fasta):
        raise FileNotFoundError(f"Myloasm did not create an assembly at {assembly_fasta}")

    if not os.path.exists(assembly_graph):
        ensure_placeholder(assembly_graph)

    if not os.path.exists(assembly_info):
        write_assembly_info(assembly_fasta, assembly_info)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run myloasm assembler")
    parser.add_argument("--input-fastq", required=True, help="Input long read FASTQ")
    parser.add_argument("--output-dir", default="data/myloasm", help="Output directory")
    parser.add_argument("--long-read-type", default="ont", help="Long read technology")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    parser.add_argument("--log", required=True, help="Log file")

    args = parser.parse_args()

    run_myloasm(
        input_fastq=args.input_fastq,
        output_dir=args.output_dir,
        long_read_type=args.long_read_type,
        threads=args.threads,
        log=args.log,
    )
