#!/usr/bin/env python3

import os
import argparse
from subprocess import Popen, PIPE
from pathlib import Path
from typing import List

def run_skip_qc(
    long_reads: List[str],
    output_long_reads: str,
    coassemble: bool,
    log_file: str,
    threads: int
):
    """
    Skip quality control
    :param long_reads: list of long reads
    :param output_long_reads: output long reads
    :param coassemble: coassemble or not, if true we will filter all reads into the same file
    :param log_file: file to write logs to
    :param threads: number of threads for compression
    :return:
    """
    with open(log_file, 'a') as logf:
        logf.write(f"Skipping quality control\n")
        logf.write(f"Coassemble: {coassemble}\n")
        if coassemble:
            # cat files into the same file
            with open(output_long_reads, 'w') as out:
                for reads in long_reads:
                    if not os.path.exists(reads):
                        logf.write(f"Long read file {reads} does not exist\n")
                        exit(1)

                    cat_or_zcat = 'zcat' if reads.endswith('.gz') else 'cat'
                    cat_cmd = f'{cat_or_zcat} {reads}'.split()
                    pigz_cmd = f'pigz -p {threads}'.split()

                    logf.write(f"Shell style : {' '.join(cat_cmd)} | {' '.join(pigz_cmd)} > {output_long_reads}\n")

                    cat_p1 = Popen(cat_cmd, stdout=PIPE, stderr=logf)
                    pigz_p2 = Popen(pigz_cmd, stdin=cat_p1.stdout, stdout=out, stderr=logf)

                    pigz_p2.wait()
                    cat_p1.wait()
                    logf.write(f"cat return: {cat_p1.returncode}\n")
                    logf.write(f"pigz return: {pigz_p2.returncode}\n")


        elif len(long_reads) > 0:
            # sym link the first file as the output file
            logf.write(f"Symlinking {long_reads[0]} to {output_long_reads}\n")
            if not os.path.exists(long_reads[0]):
                logf.write(f"Long read file {long_reads[0]} does not exist\n")
                exit(1)

            os.symlink(long_reads[0], output_long_reads)
        else:
            logf.write(f"No long reads to quality control\n")
            Path(output_long_reads).touch()

def qc_long_reads(
    long_reads: List[str], # long reads to quality control
    reference_filter: List[str], # one or more references to filter against
    coassemble: bool,
    skip_qc: bool,
    min_length: int,
    min_quality: int,
    keep_percent: int,
    threads: int,
    output_long_reads: str,
    log_file: str,
):
    """
    Quality control long reads using chopper
    :param long_reads: list of long reads
    :param reference_filter: list of reference genomes to filter against
    :param coassemble: coassemble or not, if true we will filter all reads into the same file
    :param skip_qc: skip quality control if true
    :param min_length: minimum length of reads to keep
    :param min_quality: minimum mean quality of reads to keep
    :param keep_percent: percent of reads to keep
    :param threads: number of threads
    :param output_long_reads: output long reads
    :param log_file: file to write logs to
    :return:
    """

    if skip_qc or len(long_reads) == 0:
        run_skip_qc(long_reads, output_long_reads, coassemble, log_file, threads)
        return

    # if we have more than one reference_filter file, we need to concatenate them into a single temp file
    with open(log_file, 'a') as logf:
        reference_filter_file_string = ''
        if len(reference_filter) > 1:
            with open(f'{output_long_reads}.reference_filter.fasta', 'w') as out:
                for reference in reference_filter:
                    # check if file exists
                    if not os.path.exists(reference):
                        logf.write(f"Reference filter file {reference} does not exist\n")
                        exit(1)

                    # concatenate accoutning for gzipped files
                    cat_or_zcat = 'zcat' if reference.endswith('.gz') else 'cat'
                    cat_cmd = f'{cat_or_zcat} {reference}'.split()

                    logf.write(f"Shell style : {' '.join(cat_cmd)} > {output_long_reads}\n")

                    cat_p1 = Popen(cat_cmd, stdout=out, stderr=logf)
                    cat_p1.wait()
                    logf.write(f"cat return: {cat_p1.returncode}\n")

            # gzip the concatenated file
            pigz_cmd = f'pigz -p {threads} {output_long_reads}.reference_filter.fasta'.split()
            logf.write(f"Shell style : {' '.join(pigz_cmd)}\n")

            pigz_p1 = Popen(pigz_cmd, stderr=logf)
            pigz_p1.wait()
            logf.write(f"pigz return: {pigz_p1.returncode}\n")

            reference_filter_file_string = f'-c {output_long_reads}.reference_filter.fasta.gz'
        elif len(reference_filter) == 1:
            # make sure file exists
            if not os.path.exists(reference_filter[0]):
                logf.write(f"Reference filter file {reference_filter[0]} does not exist\n")
                exit(1)

            reference_filter_file_string = f'-c {reference_filter[0]}'


        # run chopper    
        # chopper reads on stdin and write to stdout
        logf.write(f"Running chopper on {len(long_reads)} files\n")
        logf.write(f"Reference filter: {reference_filter_file_string}\n")
        with open(output_long_reads, 'w') as out:
            for long_read in long_reads:
                # make sure file exists
                if not os.path.exists(long_read):
                    logf.write(f"Long read file {long_read} does not exist\n")
                    exit(1)

                cat_or_zcat = 'zcat' if long_read.endswith('.gz') else 'cat'
                cat_cmd = f'{cat_or_zcat} {long_read}'.split()
                chopper_cmd = f'chopper -l {min_length} -q {min_quality} -t {threads} {reference_filter_file_string}'.split()
                pigz_cmd = f'pigz -p {threads}'.split()

                logf.write(f"Shell style : {' '.join(cat_cmd)} | {' '.join(chopper_cmd)} | {' '.join(pigz_cmd)} > {output_long_reads}\n")
                
                cat_p1 = Popen(cat_cmd, stdout=PIPE, stderr=logf)
                chopper_p2 = Popen(chopper_cmd, stdin=cat_p1.stdout, stdout=PIPE, stderr=logf)
                pigz_p3 = Popen(pigz_cmd, stdin=chopper_p2.stdout, stdout=out, stderr=logf)

                pigz_p3.wait()
                chopper_p2.wait()
                cat_p1.wait()
                logf.write(f"cat return: {cat_p1.returncode}\n")
                logf.write(f"chopper return: {chopper_p2.returncode}\n")
                logf.write(f"pigz return: {pigz_p3.returncode}\n")

                if not coassemble:
                    break
        
        # clean up reference filter if we concatenated
        if os.path.exists(f'{output_long_reads}.reference_filter.fasta.gz'):
            os.remove(f'{output_long_reads}.reference_filter.fasta.gz')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Quality control long reads using chopper')
    parser.add_argument('--long-reads', required=True, nargs='+', help='Long read files')
    parser.add_argument('--reference-filter', nargs='*', default=[], help='Reference genome files to filter against')
    parser.add_argument('--coassemble', type=str, choices=['True', 'False'], help='Coassemble reads into one file')
    parser.add_argument('--skip-qc', type=str, choices=['True', 'False'], help='Skip quality control')
    parser.add_argument('--min-length', type=int, default=1000, help='Minimum read length')
    parser.add_argument('--min-quality', type=int, default=10, help='Minimum mean quality')
    parser.add_argument('--keep-percent', type=int, default=90, help='Percentage of reads to keep')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--output-long', required=True, help='Output file for quality controlled long reads')
    parser.add_argument('--log', required=True, help='Log file')
    
    args = parser.parse_args()
    
    # Convert string boolean arguments to actual booleans
    coassemble = args.coassemble == 'True'
    skip_qc = args.skip_qc == 'True'
    
    # Initialize log file
    with open(args.log, "w") as logf:
        pass
    
    qc_long_reads(
        long_reads=args.long_reads,
        reference_filter=args.reference_filter,
        coassemble=coassemble,
        skip_qc=skip_qc,
        min_length=args.min_length,
        min_quality=args.min_quality,
        keep_percent=args.keep_percent,
        threads=args.threads,
        output_long_reads=args.output_long,
        log_file=args.log,
    )