import subprocess
import sys

read_set1 = snakemake.config['short_reads_1']
read_set2 = snakemake.config['short_reads_2']

# deal with read sets i.e. are we coassembling? Which assembler are we using?
# Non co-assembled reads are handled the same for each assembler
if read_set1 != 'none':
    if not snakemake.params.coassemble or len(read_set1) == 1:
        read_set1 = read_set1[0]
        if read_set2 != 'none':
            read_set2 = read_set2[0]

    elif not snakemake.params.use_megahit:
        if read_set2 != 'none':
            for reads1, reads2 in zip(read_set1, read_set2):
                subprocess.Popen(f"cat {reads1} >> data/short_reads.1.fastq.gz", shell=True).wait()
                subprocess.Popen(f"cat {reads2} >> data/short_reads.2.fastq.gz", shell=True).wait()
            read_set1 = 'data/short_reads.1.fastq.gz'
            read_set2 = 'data/short_reads.2.fastq.gz'
        else:
            for reads1 in read_set1:
                subprocess.Popen(f"cat {reads1} >> data/short_reads.1.fastq.gz", shell=True).wait()
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
if snakemake.params.use_megahit:
    command = f"megahit {read_string} -t {snakemake.threads} -m {snakemake.config['max_memory']} " \
              f"-o data/megahit_assembly --tmp-dir {snakemake.params.tmpdir}"
    print(f"Queueing command {command}")

    subprocess.Popen(
        command,
        shell=True
    ).wait()
    subprocess.Popen(
        "mkdir -p data/short_read_assembly; cp data/megahit_assembly/final.contigs.fa data/short_read_assembly/scaffolds.fasta",
        shell=True
    ).wait()
else:
    kmers = " ".join(snakemake.params.kmer_sizes)
    command = f"spades.py --memory {snakemake.config['max_memory']} --meta -t {snakemake.threads} " \
              f"-o data/short_read_assembly {read_string} -k {kmers} --tmp-dir {snakemake.params.tmpdir}"
    print(f"Queueing command {command}")

    subprocess.Popen(
        command,
        shell=True).wait()