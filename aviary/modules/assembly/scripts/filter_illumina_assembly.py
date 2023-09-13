from subprocess import Popen, run, PIPE, STDOUT
import os


def run_mapping_process(
    reads_string: str, # combination of reads1 and reads2 or just reads1
    input_fasta: str,
    output_bam: str,
    output_fastq: str,
    threads: int,
    log: str,
):
    """
    :param reads_string: combination of reads1 and reads2 or just reads1
    :param input_fasta: input fasta file
    :param output_bam: output bam file
    :param output_fastq: output fastq file
    :param threads: number of threads
    :return:
    """
    with open(log, 'a') as logf:
        minimap_cmd = f"minimap2 -ax sr -t {threads} {input_fasta} {reads_string}".split()
        samtools_view_cmd = f"samtools view -b -@ {threads}".split()
        samtools_sort_cmd = f"samtools sort -@ {threads} -o {output_bam}".split()
        logf.write(f"Shell style : {' '.join(minimap_cmd)} | {' '.join(samtools_view_cmd)} | {' '.join(samtools_sort_cmd)}\n")

        minimap_p1 = Popen(minimap_cmd, stdout=PIPE, stderr=logf) # stderr=PIPE optional, minimap2 is chatty
        samtools_view_p2 = Popen(samtools_view_cmd, stdin=minimap_p1.stdout, stdout=PIPE, stderr=logf)
        samtools_sort_p3 = Popen(samtools_sort_cmd, stdin=samtools_view_p2.stdout, stderr=logf)
        samtools_sort_p3.wait()

        # thoretically p1 and p2 may still be running, this ensures we are collecting their return codes
        minimap_p1.wait()
        samtools_view_p2.wait()
        logf.write(f"minimap return: {minimap_p1.returncode}\n")
        logf.write(f"samtools view return: {samtools_view_p2.returncode}\n")
        logf.write(f"samtools sort return: {samtools_sort_p3.returncode}\n")

        # samtools index
        samtools_index_cmd = f"samtools index -@ {threads} {output_bam}".split()
        run(samtools_index_cmd, stderr=logf)
        
        # samtools bam2fq
        samtools_bam2fq_cmd = f"samtools bam2fq -@ {threads} -f 12 {output_bam}".split()
        pigz_cmd = f"pigz -p {threads}".split()

        logf.write(f"Shell style : {' '.join(samtools_bam2fq_cmd)} | {' '.join(pigz_cmd)} > {output_fastq}\n")
        with open(output_fastq, 'w') as output_fq:
            samtools_bam2fq_p1 = Popen(samtools_bam2fq_cmd, stdout=PIPE, stderr=logf)
            pigz_p2 = Popen(pigz_cmd, stdin=samtools_bam2fq_p1.stdout, stdout=output_fq, stderr=logf)

            # thoretically p1 and p2 may still be running, this ensures we are collecting their return codes
            samtools_bam2fq_p1.wait()
            pigz_p2.wait()
            logf.write(f"samtools bam2fq return: {samtools_bam2fq_p1.returncode}\n")
            logf.write(f"pigz return: {pigz_p2.returncode}\n")


def filter_illumina_assembly(
        short_reads_1,
        short_reads_2,
        input_fasta: str,
        output_bam: str,
        output_fastq: str,
        threads: int,
        coassemble: bool,
        log: str,
):
    """
    :param short_reads_1: list of short reads 1, or 'none'
    :param short_reads_2: list of short reads 2, or 'none'
    :param input_fasta: input fasta file
    :param output_bam: output bam file
    :param output_fastq: output fastq file
    :param threads: number of threads
    :param coassemble: coassemble or not
    :return:
    """
    if os.path.exists('data/short_reads.fastq.gz'):
        run_mapping_process(
            reads_string='data/short_reads.fastq.gz',
            input_fasta=input_fasta,
            output_bam=output_bam,
            output_fastq=output_fastq,
            threads=threads,
            log=log,
        )
  
    elif short_reads_2 != 'none':
        if len(short_reads_2) == 1 or not coassemble:
            pe1 = short_reads_1[0]
            pe2 = short_reads_2[0]
        else:
            if not os.path.exists("data/short_reads.1.fastq.gz"):
                for reads1, reads2 in zip(short_reads_1, short_reads_2):
                    with open(log, "a") as logf:
                        with open("data/short_reads.1.fastq.gz", "a") as f:
                            run(f"cat {reads1}", stdout=f, stderr=logf)
                        
                        with open("data/short_reads.2.fastq.gz", "a") as f:
                            run(f"cat {reads2}", stdout=f, stderr=logf)
            pe1 = "data/short_reads.1.fastq.gz"
            pe2 = "data/short_reads.2.fastq.gz"

        reads_string = f"{pe1} {pe2}"
        run_mapping_process(
            reads_string=reads_string,
            input_fasta=input_fasta,
            output_bam=output_bam,
            output_fastq=output_fastq,
            threads=threads,
            log=log,
        )

        if os.path.exists("data/short_reads.1.fastq.gz"):
            os.remove("data/short_reads.1.fastq.gz")
            os.remove("data/short_reads.2.fastq.gz")

    elif short_reads_1 != 'none':
        if len(short_reads_1) == 1:
            pe1 = short_reads_1[0]
        else:
            if not os.path.exists("data/short_reads.1.fastq.gz") or not coassemble:
                for reads1 in short_reads_1:
                    with open(log, "a") as logf:
                        with open("data/short_reads.1.fastq.gz", "a") as f:
                            run(f"cat {reads1}", stdout=f, stderr=logf)
            pe1 = "data/short_reads.1.fastq.gz"


        run_mapping_process(
            reads_string=pe1,
            input_fasta=input_fasta,
            output_bam=output_bam,
            output_fastq=output_fastq,
            threads=threads,
            log=log,
        )

        if os.path.exists("data/short_reads.1.fastq.gz"):
            os.remove("data/short_reads.1.fastq.gz")

if __name__ == "__main__":
    short_reads_1 = snakemake.config['short_reads_1']
    short_reads_2 = snakemake.config['short_reads_2']

    input_fasta = snakemake.input.fasta
    output_bam = snakemake.output.bam
    output_fastq = snakemake.output.fastq

    threads = snakemake.threads
    coassemble = snakemake.params.coassemble
    log = snakemake.log[0]

    with open(log, "a") as logf: pass

    filter_illumina_assembly(
        short_reads_1=short_reads_1,
        short_reads_2=short_reads_2,
        input_fasta=input_fasta,
        output_bam=output_bam,
        output_fastq=output_fastq,
        threads=threads,
        coassemble=coassemble,
        log=log,
    )
