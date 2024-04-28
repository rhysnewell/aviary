import subprocess
import os
import sys
from typing import List

def spades_asssembly(
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
    :param tmpdir: temporary directory
    :param kmer_sizes: list of kmer sizes
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
    
    if os.path.exists("data/spades_assembly/tmp"):
        with open(log, 'a') as logf:
            subprocess.run("rm -rf data/spades_assembly/tmp".split(), stdout=logf, stderr=subprocess.STDOUT)   
    # remove existing temporary directory
    minimumsize=500000
    actualsize = int(subprocess.check_output('stat -c%s data/short_reads.filt.fastq.gz', shell=True))
    # check if directory exists
    if  os.path.exists("data/spades_assembly"):
        # resume previous assembly
        command = f"spades.py --restart-from last --memory {max_memory} -t {threads} " \
                f"-o data/spades_assembly -k {kmer_sizes}  {tmp_dir_arg}" 
        # run cmd
        with open(log, 'a') as logf:
            logf.write(f"Queueing command {command}\n")
            subprocess.run(command.split(), stdout=logf, stderr=subprocess.STDOUT)
            subprocess.run("cp data/spades_assembly/scaffolds.fasta data/spades_assembly.fasta".split(), stdout=logf, stderr=subprocess.STDOUT)
    elif actualsize >= minimumsize:
        if long_read_type in ["ont","ont_hq"]:
            command = f"spades.py --checkpoints all --memory {max_memory} --meta --nanopore {input_long_reads} --12 {input_fastq} "\
                f"-o data/spades_assembly -t {threads}  -k {kmer_sizes} {tmp_dir_arg} "       
        else:
            command = f"spades.py --checkpoints all --memory {max_memory} --meta --pacbio {input_long_reads} --12 {input_fastq} "\
                f"-o data/spades_assembly -t {threads}  -k {kmer_sizes} {tmp_dir_arg} "
        # run cmd
        with open(log, 'a') as logf:
            logf.write(f"Queueing command {command}\n")
            subprocess.run(command.split(), stdout=logf, stderr=subprocess.STDOUT)        
            subprocess.run("cp data/spades_assembly/scaffolds.fasta data/spades_assembly.fasta".split(), stdout=logf, stderr=subprocess.STDOUT)
    else:
        with open(log, 'a') as logf:
            subprocess.run(f"mkdir -p {output.spades_folder} && touch {output.fasta}".split(), stdout=logf, stderr=subprocess.STDOUT)
    

if __name__ == '__main__':
    log = snakemake.log[0]
    with open(log, 'w') as logf: pass

    spades_asssembly(
        snakemake.input.fastq,
        snakemake.input.long_reads,
        snakemake.output.fasta,
        snakemake.output.spades_folder,
        snakemake.params.max_memory,
        snakemake.threads,
        snakemake.params.kmer_sizes,
        snakemake.params.tmpdir,
        snakemake.params.long_read_type,
        log
    )