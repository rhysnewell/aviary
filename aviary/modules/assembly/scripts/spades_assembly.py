#!/usr/bin/env python3
import subprocess
import os
import sys
import argparse
from typing import List

def spades_assembly(
    input_fastq: str,
    input_long_reads: str,
    output_fasta: str,
    output_spades_folder: str,
    max_memory: int,
    threads: int,
    kmer_sizes: List[int],
    tmp_dir: str,
    long_read_type: str,
    log: str):

    '''
    Assemble short reads and long reads (if any) using spades
    :param input_fastq: short reads fastq file
    :param input_long_reads: long reads fastq file
    :param output_fasta: output fasta file
    :param output_spades_folder: output spades folder
    :param max_memory: maximum memory to use
    :param threads: number of threads
    :param kmer_sizes: list of kmer sizes
    :param tmp_dir: temporary directory
    :param long_read_type: type of long reads
    :param log: log file
    :return:
    '''
    if tmp_dir:
        tmp_dir_arg = f"--tmp-dir {tmp_dir}"
    else:
        try:
            tmp_dir = os.environ["TMPDIR"]
            tmp_dir_arg = f"--tmp-dir {tmp_dir}"
        except KeyError:
            tmp_dir_arg = ""
    
    kmer_str = ",".join(str(k) for k in kmer_sizes)
    
    if os.path.exists("data/spades_assembly/tmp"):
        with open(log, 'a') as logf:
            subprocess.run("rm -rf data/spades_assembly/tmp".split(), stdout=logf, stderr=subprocess.STDOUT)   
    # remove existing temporary directory
    minimumsize=500000
    actualsize = int(subprocess.check_output('stat -c%s data/short_reads.filt.fastq.gz', shell=True))
    # check if directory exists
    if os.path.exists("data/spades_assembly"):
        # resume previous assembly
        command = f"spades.py --restart-from last --memory {max_memory} -t {threads} " \
                f"-o data/spades_assembly -k {kmer_str} {tmp_dir_arg}" 
        # run cmd
        with open(log, 'a') as logf:
            logf.write(f"Queueing command {command}\n")
            subprocess.run(command.split(), stdout=logf, stderr=subprocess.STDOUT)
            subprocess.run(f"cp data/spades_assembly/scaffolds.fasta {output_fasta}".split(), stdout=logf, stderr=subprocess.STDOUT)
    elif actualsize >= minimumsize:
        if long_read_type in ["ont", "ont_hq"]:
            command = f"spades.py --checkpoints all --memory {max_memory} --meta --nanopore {input_long_reads} --12 {input_fastq} "\
                f"-o {output_spades_folder} -t {threads} -k {kmer_str} {tmp_dir_arg} "       
        else:
            command = f"spades.py --checkpoints all --memory {max_memory} --meta --pacbio {input_long_reads} --12 {input_fastq} "\
                f"-o {output_spades_folder} -t {threads} -k {kmer_str} {tmp_dir_arg} "
        # run cmd
        with open(log, 'a') as logf:
            logf.write(f"Queueing command {command}\n")
            subprocess.run(command.split(), stdout=logf, stderr=subprocess.STDOUT)        
            subprocess.run(f"cp {output_spades_folder}/scaffolds.fasta {output_fasta}".split(), stdout=logf, stderr=subprocess.STDOUT)
    else:
        with open(log, 'a') as logf:
            subprocess.run(f"mkdir -p {output_spades_folder} && touch {output_fasta}", stdout=logf, stderr=subprocess.STDOUT, shell=True)


def main():
    parser = argparse.ArgumentParser(description='Assemble reads using SPAdes')
    parser.add_argument('--input-fastq', required=True, help='Path to short reads fastq file')
    parser.add_argument('--input-long-reads', required=True, help='Path to long reads fastq file')
    parser.add_argument('--output-fasta', required=True, help='Path to output fasta file')
    parser.add_argument('--output-spades-folder', required=True, help='Path to output SPAdes folder')
    parser.add_argument('--max-memory', type=int, required=True, help='Maximum memory to use (GB)')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--kmer-sizes', nargs='+', required=True, help='Kmer sizes for assembly')
    parser.add_argument('--tmp-dir', default='', help='Temporary directory')
    parser.add_argument('--long-read-type', choices=['ont', 'ont_hq', 'pacbio', 'pacbio_hifi'], 
                      default='pacbio', help='Type of long reads')
    parser.add_argument('--log', required=True, help='Path to log file')
    
    args = parser.parse_args()
    
    with open(args.log, "w") as logf:
        pass  # Initialize empty log file

    spades_assembly(
        args.input_fastq,
        args.input_long_reads,
        args.output_fasta,
        args.output_spades_folder,
        args.max_memory,
        args.threads,
        args.kmer_sizes,
        args.tmp_dir,
        args.long_read_type,
        args.log
    )


if __name__ == "__main__":
    main()