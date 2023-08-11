from subprocess import run, Popen
import multiprocessing as mp
import os
from pathlib import Path


def spawn_fastqc(file, threads=1):
    fastq_cmd = f"fastqc -o www/fastqc/ -t {threads} {file}".split()
    run(fastq_cmd)


def run_fastqc(
    short_reads_1,
    short_reads_2,
    threads: int,
):
    if os.path.exists('data/short_reads.fastq.gz'):
        fastq_cmd = f"fastqc -o www/fastqc/ -t {threads} data/short_reads.fastq.gz".split()
        run(fastq_cmd)

    elif short_reads_2 != 'none': # paired end
        pool = mp.Pool(threads)
        if isinstance(short_reads_1, str):
            reads = [short_reads_2, short_reads_2]
        else:
            threads = max(len(short_reads_1 + short_reads_2) // threads, 1)
            reads = short_reads_1 + short_reads_2
        mp_results = [pool.apply_async(spawn_fastqc, args=(read, threads))
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

        mp_results = [pool.apply_async(spawn_fastqc, args=(read, threads))
                    for read in
                    reads]

        for result in mp_results:
            result.get()

        pool.close()
        pool.join()
    else:
        echo_cmd = 'echo "no short reads"'.split()
        
        with open('www/fastqc/short_reads_fastqc.html', 'w') as f:
            Popen(echo_cmd, stdout=f).wait()

    Path('www/fastqc/done').touch()


if __name__ == '__main__':
    run_fastqc(
        snakemake.config['short_reads_1'],
        snakemake.config['short_reads_2'],
        snakemake.threads,
    )
