#!/usr/bin/env python3

import argparse
import os
import subprocess
import tempfile

# myloasm's own default kmer size (see `myloasm <file> --markdown-help`); kept in
# sync here since aviary doesn't expose --kmer-size for myloasm itself.
KMER_SIZE = 21


def run_myloasm_kmc(
    input_fastq: str, db_base: str, threads: int, log: str,
    tmp_dir: str | None = None, ram_gb: int | None = None,
) -> None:
    if tmp_dir is None:
        tmp_dir = os.path.join(tempfile.gettempdir(), "aviary_myloasm_kmc")
    os.makedirs(tmp_dir, exist_ok=True)
    os.makedirs(os.path.dirname(db_base) or ".", exist_ok=True)

    cmd = [
        "myloasm-kmc-v1",
        "-o", db_base,
        "--tmp-dir", tmp_dir,
        "-k", str(KMER_SIZE),
        "-t", str(threads),
    ]
    if ram_gb is not None:
        cmd.extend(["--ram-gb", str(ram_gb)])
    cmd.append(input_fastq)

    with open(log, "w") as logf:
        subprocess.run(cmd, check=True, stdout=logf, stderr=subprocess.STDOUT)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count k-mers on disk with myloasm-kmc-v1")
    parser.add_argument("--input-fastq", required=True, help="Input long read FASTQ")
    parser.add_argument("--db-base", required=True, help="Output database base path (produces <path>.kmc_pre/.kmc_suf)")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    parser.add_argument("--tmp-dir", default=None, help="Existing (or creatable) directory for counting scratch space")
    parser.add_argument("--ram-gb", type=int, default=None, help="RAM budget in GB")
    parser.add_argument("--log", required=True, help="Log file")

    args = parser.parse_args()

    run_myloasm_kmc(
        input_fastq=args.input_fastq,
        db_base=args.db_base,
        threads=args.threads,
        log=args.log,
        tmp_dir=args.tmp_dir,
        ram_gb=args.ram_gb,
    )
