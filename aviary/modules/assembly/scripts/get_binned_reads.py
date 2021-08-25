import subprocess
import os
import multiprocessing as mp
import glob
import logging

def process_short_reads_interleaved_filtered(reads, file):
    subprocess.Popen('seqtk seq -1 %s | mfqe --fastq-read-name-lists %s --output-fastq-files %s.1.fastq.gz' %
                     (reads, file, file[:-5]), shell=True).wait()
    subprocess.Popen('seqtk seq -2 %s | mfqe --fastq-read-name-lists %s --output-fastq-files %s.2.fastq.gz' %
                     (reads, file, file[:-5]), shell=True).wait()

def process_short_reads_interleaved(reads, file, threads):
    subprocess.Popen('seqtk seq -1 %s | seqkit -j %d grep -f %s | pigz -p %d > %s.1.fastq.gz' %
                     (reads, threads, file, threads, file[:-5]), shell=True).wait()
    subprocess.Popen('seqtk seq -2 %s | seqkit -j %d grep -f %s | pigz -p %d > %s.2.fastq.gz' %
                     (reads, threads, file, threads, file[:-5]), shell=True).wait()

def process_short_reads_paired(reads_1, reads_2, file, threads):
    subprocess.Popen('seqkit -j %d grep -f %s %s | pigz -p %d > %s.1.fastq.gz' %
                     (threads, file, reads_1, threads, file[:-5]), shell=True).wait()
    subprocess.Popen('seqkit -j %d grep -f %s %s | pigz -p %d > %s.2.fastq.gz' %
                     (threads, file, reads_2, threads, file[:-5]), shell=True).wait()

def process_long_reads(reads, file, threads):
    subprocess.Popen('seqtk subseq %s %s | pigz -p %d > %s.fastq.gz' %
                     (reads, file, threads, file[:-5]), shell=True).wait()


def get_index(n_files, current):
    if current >= n_files / 2:
        return 2
    else:
        return 1


if snakemake.config['long_reads'] != 'none':
    logging.info("Extracting binned long reads...")

    pool = mp.Pool(snakemake.threads)

    mp_results = [pool.apply_async(process_long_reads,
                                    args=(snakemake.config["long_reads"][0], file, 2))
                                    for file in glob.glob('data/binned_reads/*.long.list')]
    for result in mp_results:
        result.get()

    pool.close()
    pool.join()

if os.path.exists('data/short_reads.fastq.gz'):
    logging.info("Extracting binned interleaved filtered short reads...")
    pool = mp.Pool(snakemake.threads)
    # extract reads from filtered short reads
    # take end 1 and 2 separately
    mp_results = [pool.apply_async(process_short_reads_interleaved_filtered,
                                    args=('data/short_reads.fastq.gz', file))
                                    for file in glob.glob('data/binned_reads/*.short.list')]
    for result in mp_results:
        result.get()

    pool.close()
    pool.join()

elif snakemake.config['short_reads_2'] != 'none':
    logging.info("Extracting binned paired end short reads...")
    pool = mp.Pool(snakemake.threads // 2)
    # extract reads from all short read files
    mp_results = [pool.apply_async(process_short_reads_paired,
                                    args=(' '.join(snakemake.config['short_reads_1']),
                                          ' '.join(snakemake.config['short_reads_2']),
                                          file, 2))
                                    for file in glob.glob('data/binned_reads/*.short.list')]

    for result in mp_results:
        result.get()

    pool.close()
    pool.join()

elif snakemake.config['short_reads_1'] != 'none':
    logging.info("Extracting binned interleaved short reads...")
    pool = mp.Pool(snakemake.threads // 2)
    # extract reads from all short read files
    mp_results = [pool.apply_async(process_short_reads_interleaved,
                                   args=(' '.join(snakemake.config['short_reads_1']),
                                         file, 2))
                  for file in glob.glob('data/binned_reads/*.short.list')]

    for result in mp_results:
        result.get()

    pool.close()
    pool.join()

subprocess.Popen('touch %s' % (snakemake.output), shell=True).wait()

