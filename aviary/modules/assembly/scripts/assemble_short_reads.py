import subprocess
import os
import shutil
import sys
from typing import List



def assemble_short_reads(
        read_set1, read_set2,
        reference_filter: str,
        max_memory: int,
        use_megahit: bool,
        coassemble: bool,
        threads: int,
        tmp_dir: str,
        kmer_sizes: List[int],
        log: str):
    '''
    Assemble short reads using either megahit or spades
    :param read_set1: list of short reads 1, or 'none'
    :param read_set2: list of short reads 2, or 'none'
    :param max_memory: maximum memory to use
    :param use_megahit: use megahit or not
    :param coassemble: coassemble or not
    :param threads: number of threads
    :param tmp_dir: temporary directory
    :param kmer_sizes: list of kmer sizes
    :return:
    '''

    if not tmp_dir:
        try:
            tmp_dir = os.environ["TMPDIR"]
        except KeyError:
            tmp_dir_arg = ""

    if tmp_dir:
        tmp_dir_arg = f"--tmp-dir {tmp_dir}"

    # deal with read sets i.e. are we coassembling? Which assembler are we using?
    # Non co-assembled reads are handled the same for each assembler
    if read_set1 != 'none':
        if not coassemble or len(read_set1) == 1:
            read_set1 = read_set1[0]
            if read_set2 != 'none':
                read_set2 = read_set2[0]

        elif not use_megahit and reference_filter == 'none':
            if read_set2 != 'none':
                for reads1, reads2 in zip(read_set1, read_set2):
                    with open(log, 'a') as logf:
                        with open('data/short_reads.1.fastq.gz', 'a') as out1:
                            subprocess.run(f"cat {reads1}".split(), stdout=out1, stderr=logf)

                        with open('data/short_reads.2.fastq.gz', 'a') as out2:
                            subprocess.run(f"cat {reads2}".split(), stdout=out2, stderr=logf)

                read_set1 = 'data/short_reads.1.fastq.gz'
                read_set2 = 'data/short_reads.2.fastq.gz'
            else:
                for reads1 in read_set1:
                    with open(log, 'a') as logf:
                        with open('data/short_reads.1.fastq.gz', 'a') as out1:
                            subprocess.run(f"cat {reads1}".split(), stdout=out1, stderr=logf)

                read_set1 = 'data/short_reads.1.fastq.gz'

        else:
            read_set1 = ",".join(read_set1)
            if read_set2 != 'none':
                read_set2 = ",".join(read_set2)
    else:
        # forward reads must always be present as either
        # forward or interleaved or single ended
        with open(log, 'a') as logf:
            logf.write(f"============= ERROR =============\n\n")
            logf.write(f"Invalid read sets provided \n\n "
                f"for short read assembly: \n\n")
            logf.write(f"    Set 1: {read_set1}\n\n")
            logf.write(f"    Set 2: {read_set2}\n\n")
            logf.write("Validate read sets and resubmit\n")
        sys.exit(1)


    # designate input read string
    read_string = f"--12 {read_set1}"
    if read_set2 != 'none':
        read_string = f"-1 {read_set1} -2 {read_set2}"
    
    if reference_filter != 'none':
        read_string = f"--12 data/short_reads.fastq.gz"


    # Run chosen assembler
    if use_megahit:
        max_memory_in_bytes = max_memory * 1024*1024*1024
        command = f"megahit {read_string} -t {threads} -m {max_memory_in_bytes} -o data/megahit_assembly {tmp_dir_arg}"

        with open(log, 'a') as logf:
            logf.write(f"Queueing command {command}\n")
            subprocess.run(command.split(), stdout=logf, stderr=subprocess.STDOUT)
        os.makedirs("data/short_read_assembly", exist_ok=True)
        shutil.copyfile("data/megahit_assembly/final.contigs.fa", "data/short_read_assembly/scaffolds.fasta")

    else:
        kmers = " ".join(kmer_sizes)
        command = f"spades.py --memory {max_memory} --meta -t {threads} " \
                f"-o data/short_read_assembly {read_string} -k {kmers} {tmp_dir_arg}"
        with open(log, 'a') as logf:
            logf.write(f"Queueing command {command}\n")
            subprocess.run(command.split(), stdout=logf, stderr=subprocess.STDOUT)


if __name__ == '__main__':
    read_set1 = snakemake.config['short_reads_1']
    read_set2 = snakemake.config['short_reads_2']
    log = snakemake.log[0]

    with open(log, 'w') as logf: pass

    assemble_short_reads(
        read_set1,
        read_set2,
        snakemake.config['reference_filter'],
        snakemake.config['max_memory'],
        snakemake.params.use_megahit,
        snakemake.params.coassemble,
        snakemake.threads,
        snakemake.params.tmpdir,
        snakemake.params.kmer_sizes,
        log,
    )