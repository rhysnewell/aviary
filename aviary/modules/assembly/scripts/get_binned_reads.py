from subprocess import Popen, PIPE, run
import os
import multiprocessing as mp
import glob
import logging
from pathlib import Path

def run_mfqe(
    read_pair: str, # either '1' or '2'
    input_reads: str,
    output_file: str,
):
    output_prefix = output_file[:-5]
    seqtk_cmd = f"seqtk seq -{read_pair} {input_reads}".split()
    mfqe_cmd = f"mfqe --fastq-read-name-lists {output_file} --output-fastq-files {output_prefix}.{read_pair}.fastq.gz".split()

    seqtk = Popen(seqtk_cmd, stdout=PIPE)
    mfqe = Popen(mfqe_cmd, stdin=seqtk.stdout)
    mfqe.wait()

    # thoretically p1 and p2 may still be running, this ensures we are collecting their return codes
    seqtk.wait()
    print("seqtk return: ", seqtk.returncode)
    print("mfqe return: ", mfqe.returncode)

def process_short_reads_interleaved_filtered(reads, file):
    run_mfqe(read_pair='1', input_reads=reads, output_file=file)
    run_mfqe(read_pair='2', input_reads=reads, output_file=file)


def run_seqkit(
    read_pair: str, # either '1' or '2'
    reads: str,
    output_file: str, # also the pattern file for seqkit
    threads: int
):
    output_prefix = output_file[:-5]
    with open(f"{output_prefix}.{read_pair}.fastq.gz", 'w') as out:
        seqtk_cmd = f"seqtk seq -{read_pair} {reads}".split()
        seqkit_cmd = f"seqkit -j {threads} grep -f {output_file}".split()
        pigz_cmd = f"pigz -p {threads}".split()

        seqtk = Popen(seqtk_cmd, stdout=PIPE)
        seqkit = Popen(seqkit_cmd, stdin=seqtk.stdout, stdout=PIPE)
        pigz = Popen(pigz_cmd, stdin=seqkit.stdout, stdout=out)
        pigz.wait()

        seqtk.wait()
        seqkit.wait()
        print("seqtk return: ", seqtk.returncode)
        print("seqkit return: ", seqkit.returncode)
        print("pigz return: ", pigz.returncode)


def process_short_reads_interleaved(reads, file, threads):
    run_seqkit(read_pair='1', reads=reads, output_file=file, threads=threads)
    run_seqkit(read_pair='2', reads=reads, output_file=file, threads=threads)


def run_seqkit_without_seqtk(
    read_pair: str, # either '1' or '2'
    reads: str,
    output_file: str, # also the pattern file for seqkit
    threads: int
):
    output_prefix = output_file[:-5]
    with open(f"{output_prefix}.{read_pair}.fastq.gz", 'w') as out:
        seqkit_cmd = f"seqkit -j {threads} grep -f {output_file} {reads}".split()
        pigz_cmd = f"pigz -p {threads}".split()

        seqkit = Popen(seqkit_cmd, stdout=PIPE)
        pigz = Popen(pigz_cmd, stdin=seqkit.stdout, stdout=out)
        pigz.wait()

        seqkit.wait()
        print("seqkit return: ", seqkit.returncode)
        print("pigz return: ", pigz.returncode)

def process_short_reads_paired(reads_1, reads_2, file, threads):
    run_seqkit_without_seqtk(read_pair='1', reads=reads_1, output_file=file, threads=threads)
    run_seqkit_without_seqtk(read_pair='2', reads=reads_2, output_file=file, threads=threads)



def process_long_reads(reads, file, threads):

    output_prefix = file[:-5]
    with open(f"{output_prefix}.fastq.gz", 'w') as out:
        seqtk_cmd = f"seqtk subseq {reads} {file}".split()
        pigz_cmd = f"pigz -p {threads}".split()

        seqtk = Popen(seqtk_cmd, stdout=PIPE)
        pigz = Popen(pigz_cmd, stdin=seqtk.stdout, stdout=out)
        pigz.wait()

        seqtk.wait()
        print("seqtk return: ", seqtk.returncode)
        print("pigz return: ", pigz.returncode)


def get_index(n_files, current):
    if current >= n_files / 2:
        return 2
    else:
        return 1

def get_binned_reads(
    long_reads,
    short_reads_1,
    short_reads_2,
    threads: int,
    output_file: str
):
    if long_reads != 'none':
        logging.info("Extracting binned long reads...")

        pool = mp.Pool(threads)

        mp_results = [pool.apply_async(process_long_reads,
                                        args=(long_reads[0], file, 2))
                                        for file in glob.glob('data/binned_reads/*.long.list')]
        for result in mp_results:
            result.get()

        pool.close()
        pool.join()

    if os.path.exists('data/short_reads.fastq.gz'):
        logging.info("Extracting binned interleaved filtered short reads...")
        pool = mp.Pool(threads)
        # extract reads from filtered short reads
        # take end 1 and 2 separately
        mp_results = [pool.apply_async(process_short_reads_interleaved_filtered,
                                        args=('data/short_reads.fastq.gz', file))
                                        for file in glob.glob('data/binned_reads/*.short.list')]
        for result in mp_results:
            result.get()

        pool.close()
        pool.join()

    elif short_reads_2 != 'none':
        logging.info("Extracting binned paired end short reads...")
        pool = mp.Pool(threads // 2)
        # extract reads from all short read files
        mp_results = [pool.apply_async(process_short_reads_paired,
                                        args=(' '.join(short_reads_1),
                                            ' '.join(short_reads_2),
                                            file, 2))
                                        for file in glob.glob('data/binned_reads/*.short.list')]

        for result in mp_results:
            result.get()

        pool.close()
        pool.join()

    elif short_reads_1 != 'none':
        logging.info("Extracting binned interleaved short reads...")
        pool = mp.Pool(threads // 2)
        # extract reads from all short read files
        mp_results = [pool.apply_async(process_short_reads_interleaved,
                                    args=(' '.join(short_reads_1),
                                            file, 2))
                    for file in glob.glob('data/binned_reads/*.short.list')]

        for result in mp_results:
            result.get()

        pool.close()
        pool.join()

    Path(output_file).touch()


if __name__ == '__main__':
    get_binned_reads(
        snakemake.config['long_reads'],
        snakemake.config['short_reads_1'],
        snakemake.config['short_reads_2'],
        snakemake.threads,
        snakemake.output
    )
