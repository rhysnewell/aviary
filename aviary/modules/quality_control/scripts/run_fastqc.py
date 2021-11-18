import subprocess
import multiprocessing as mp
import os
import logging


def spawn_fastqc(file, threads=1):
    subprocess.Popen('fastqc -o www/fastqc/ -t %d %s' % (threads, file), shell=True).wait()


if os.path.exists('data/short_reads.fastq.gz'):
    subprocess.Popen("fastqc -o www/fastqc/ -t %d data/short_reads.fastq.gz" % (snakemake.threads))
elif snakemake.config['short_reads_2'] != 'none': # paired end
    pool = mp.Pool(snakemake.threads)
    threads = max(len(snakemake.config['short_reads_1'] + snakemake.config['short_reads_2']) // snakemake.threads, 1)
    mp_results = [pool.apply_async(spawn_fastqc, args=(reads, threads))
                  for reads in
                  snakemake.config['short_reads_1'] + snakemake.config['short_reads_2']]

    for result in mp_results:
        result.get()

    pool.close()
    pool.join()
elif snakemake.config['short_reads_1'] != 'none': # interleaved
    pool = mp.Pool(snakemake.threads)
    threads = max(len(snakemake.config['short_reads_1'] + snakemake.config['short_reads_2']) // snakemake.threads, 1)

    mp_results = [pool.apply_async(spawn_fastqc, args=(reads, threads))
                  for reads in
                  snakemake.config['short_reads_1']]

    for result in mp_results:
        result.get()

    pool.close()
    pool.join()
else:
    subprocess.Popen('echo "no short reads" > www/fastqc/short_reads_fastqc.html', shell=True).wait()

subprocess.Popen('touch www/fastqc/done', shell=True).wait()
