#!/usr/bin/env python3
from subprocess import run, Popen, STDOUT
import os
import argparse
from pathlib import Path


def spawn_rastqc(files, log, threads=1):
    rastqc_cmd = ['rastqc', '-o', 'www/fastqc/', '-t', str(threads)] + files
    with open(log, "a") as logf:
        run(rastqc_cmd, stdout=logf, stderr=STDOUT)


def run_rastqc(
    short_reads_1,
    short_reads_2,
    threads: int,
    log: str,
):
    with open(log, "w") as logf: pass

    if os.path.exists('data/short_reads.fastq.gz'):
        rastqc_cmd = f"rastqc -o www/fastqc/ -t {threads} data/short_reads.fastq.gz".split()
        with open(log, "a") as logf:
            run(rastqc_cmd, stdout=logf, stderr=STDOUT)

    elif short_reads_2 != 'none': # paired end
        threads = max(len(short_reads_1 + short_reads_2) // threads, 1)
        reads = short_reads_1 + short_reads_2
        spawn_rastqc(reads, log, threads)
    elif short_reads_1 != 'none': # interleaved
        if isinstance(short_reads_1, str):
            reads = [short_reads_1]
        else:
            threads = max(len(short_reads_1) // threads, 1)
            reads = short_reads_1
        spawn_rastqc(reads, log, threads)
    else:
        echo_cmd = 'echo "no short reads"'.split()

        with open(log, 'a') as f:
            Popen(echo_cmd, stdout=f).wait()

        with open('www/fastqc/short_reads_fastqc.html', 'w') as f:
            Popen(echo_cmd, stdout=f).wait()

    Path('www/fastqc/done').touch()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run RastQC on reads')
    parser.add_argument('--short-reads-1', required=True, nargs='+', help='Path to short reads file 1 or "none"')
    parser.add_argument('--short-reads-2', default='none', nargs='+', help='Path to short reads file 2 or "none"')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    parser.add_argument('--log', required=True, help='Path to log file')

    args = parser.parse_args()

    run_rastqc(
        args.short_reads_1,
        args.short_reads_2,
        args.threads,
        args.log,
    )
