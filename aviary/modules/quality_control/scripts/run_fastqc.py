from subprocess import run, Popen, STDOUT
import multiprocessing as mp
import os
from pathlib import Path


def spawn_fastqc(file, log, threads=1):
    fastq_cmd = f"fastqc -o www/fastqc/ -t {threads} {file}".split()
    with open(log, "a") as logf:
        run(fastq_cmd, stdout=logf, stderr=STDOUT)


def run_fastqc(
    short_reads_1,
    short_reads_2,
    threads: int,
    log: str,
):
    with open(log, "w") as logf: pass

    if os.path.exists('data/short_reads.fastq.gz'):
        fastq_cmd = f"fastqc -o www/fastqc/ -t {threads} data/short_reads.fastq.gz".split()
        with open(log, "a") as logf:
            run(fastq_cmd, stdout=logf, stderr=STDOUT)

    elif short_reads_2 != 'none': # paired end
        pool = mp.Pool(threads)
        if isinstance(short_reads_1, str):
            reads = [short_reads_1, short_reads_2]
        else:
            threads = max(len(short_reads_1 + short_reads_2) // threads, 1)
            reads = short_reads_1 + short_reads_2

        mp_results = [pool.apply_async(spawn_fastqc, args=(read, log, threads))
                    for read in
                    reads]

        for result in mp_results:
            result.get()

        pool.close()
        pool.join()
    elif short_reads_1 != 'none': # interleaved
        pool = mp.Pool(threads)
        if isinstance(short_reads_1, str):
            reads = [short_reads_1]
        else:
            threads = max(len(short_reads_1) // threads,
                        1)
            reads = short_reads_1

        mp_results = [pool.apply_async(spawn_fastqc, args=(read, log, threads))
                    for read in
                    reads]

        for result in mp_results:
            result.get()

        pool.close()
        pool.join()
    else:
        echo_cmd = 'echo "no short reads"'.split()

        with open(log, 'a') as f:
            Popen(echo_cmd, stdout=f).wait()

        with open('www/fastqc/short_reads_fastqc.html', 'w') as f:
            Popen(echo_cmd, stdout=f).wait()

    Path('www/fastqc/done').touch()


if __name__ == '__main__':
    run_fastqc(
        snakemake.config['short_reads_1'],
        snakemake.config['short_reads_2'],
        snakemake.threads,
        snakemake.log[0],
    )
