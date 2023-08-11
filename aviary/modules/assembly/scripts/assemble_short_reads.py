import subprocess
import os
import shutil
import sys
from typing import List



def assemble_short_reads(
        read_set1, read_set2,
        max_memory: int,
        use_megahit: bool,
        coassemble: bool,
        threads: int,
        tmp_dir: str,
        kmer_sizes: List[int]):
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

    # deal with read sets i.e. are we coassembling? Which assembler are we using?
    # Non co-assembled reads are handled the same for each assembler
    if read_set1 != 'none':
        if not coassemble or len(read_set1) == 1:
            read_set1 = read_set1[0]
            if read_set2 != 'none':
                read_set2 = read_set2[0]

        elif not use_megahit:
            if read_set2 != 'none':
                for reads1, reads2 in zip(read_set1, read_set2):
                    with open('data/short_reads.1.fastq.gz', 'a') as out1:
                        subprocess.run(f"cat {reads1}".split(), stdout=out1)

                    with open('data/short_reads.2.fastq.gz', 'a') as out2:
                        subprocess.run(f"cat {reads2}".split(), stdout=out2)

                read_set1 = 'data/short_reads.1.fastq.gz'
                read_set2 = 'data/short_reads.2.fastq.gz'
            else:
                for reads1 in read_set1:
                    with open('data/short_reads.1.fastq.gz', 'a') as out1:
                        subprocess.run(f"cat {reads1}".split(), stdout=out1)

                read_set1 = 'data/short_reads.1.fastq.gz'

        else:
            read_set1 = ",".join(read_set1)
            if read_set2 != 'none':
                read_set2 = ",".join(read_set2)
    else:
        # forward reads must always be present as either
        # forward or interleaved or single ended
        print(f"============= ERROR =============\n")
        print(f"Invalid read sets provided \n "
            f"for short read assembly: \n")
        print(f"    Set 1: {read_set1}\n")
        print(f"    Set 2: {read_set2}\n")
        print("Validate read sets and resubmit")
        sys.exit(1)


    # designate input read string
    read_string = f"--12 {read_set1}"
    if read_set2 != 'none':
        read_string = f"-1 {read_set1} -2 {read_set2}"


    # Run chosen assembler
    if use_megahit:

        command = f"megahit {read_string} -t {threads} -m {max_memory} -o data/megahit_assembly --tmp-dir {tmp_dir}"
        print(f"Queueing command {command}")

        subprocess.run(command.split())
        os.makedirs("data/short_read_assembly", exist_ok=True)
        shutil.copyfile("data/megahit_assembly/final.contigs.fa", "data/short_read_assembly/scaffolds.fasta")

    else:
        kmers = " ".join(kmer_sizes)
        command = f"spades.py --memory {max_memory} --meta -t {threads} " \
                f"-o data/short_read_assembly {read_string} -k {kmers} --tmp-dir {tmp_dir}"
        print(f"Queueing command {command}")
        subprocess.run(command.split())


if __name__ == '__main__':
    read_set1 = snakemake.config['short_reads_1']
    read_set2 = snakemake.config['short_reads_2']

    assemble_short_reads(
        read_set1,
        read_set2,
        snakemake.config['max_memory'],
        snakemake.params.use_megahit,
        snakemake.params.coassemble,
        snakemake.threads,
        snakemake.params.tmpdir,
        snakemake.params.kmer_sizes,
    )