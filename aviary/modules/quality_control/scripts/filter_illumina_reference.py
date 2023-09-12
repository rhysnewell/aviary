from subprocess import Popen, PIPE, run, STDOUT
from pathlib import Path
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
    with open(log, "a") as logf:
        minimap_cmd = f"minimap2 -ax sr -t {threads} {input_fasta} {reads_string}".split()
        samtools_view_cmd = f"samtools view -b -f 12 -@ {threads} -o {output_bam}".split()
        logf.write(f"Shell style : {' '.join(minimap_cmd)} | {' '.join(samtools_view_cmd)}\n")

        minimap_p1 = Popen(minimap_cmd, stdout=PIPE, stderr=logf) # stderr=PIPE optional, dd is chatty
        samtools_view_p2 = Popen(samtools_view_cmd, stdin=minimap_p1.stdout, stderr=logf)
        samtools_view_p2.wait()

        # thoretically p1 and p2 may still be running, this ensures we are collecting their return codes
        minimap_p1.wait()
        logf.write(f"minimap return: {minimap_p1.returncode}\n")
        logf.write(f"samtools view return: {samtools_view_p2.returncode}\n")

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

def filter_illumina_reference(
    short_reads_1,
    short_reads_2,
    reference_filter: str,
    output_bam: str,
    output_fastq: str,
    threads: int,
    coassemble: bool,
    filtered: str,
    log: str,
):
    if os.path.exists('data/short_reads.fastq.gz'):
        run_mapping_process(
            reads_string='data/short_reads.fastq.gz',
            input_fasta=reference_filter,
            output_bam=output_bam,
            output_fastq=output_fastq,
            threads=threads,
            log=log,
        )

    elif short_reads_1 != 'none':
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
            input_fasta=reference_filter,
            output_bam=output_bam,
            output_fastq=output_fastq,
            threads=threads,
            log=log,
        )
        if os.path.exists("data/short_reads.1.fastq.gz"):
            os.remove("data/short_reads.1.fastq.gz")
            os.remove("data/short_reads.2.fastq.gz")

    elif short_reads_1 != 'none':
        if len(short_reads_1) == 1 or not coassemble:
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
            input_fasta=reference_filter,
            output_bam=output_bam,
            output_fastq=output_fastq,
            threads=threads,
            log=log,
        )
        if os.path.exists("data/short_reads.1.fastq.gz"):
            os.remove("data/short_reads.1.fastq.gz")

    Path(filtered).touch()

if __name__ == '__main__':
    reference_filter = snakemake.input.reference_filter
    short_reads_1 = snakemake.config['short_reads_1']
    short_reads_2 = snakemake.config['short_reads_2']
    output_bam = snakemake.output.bam
    output_fastq = snakemake.output.fastq
    output_filtered = snakemake.output.filtered
    coassemble = snakemake.params.coassemble

    threads = snakemake.threads
    log = snakemake.log[0]

    with open(log, "w") as logf: pass

    filter_illumina_reference(
        short_reads_1=short_reads_1,
        short_reads_2=short_reads_2,
        reference_filter=reference_filter,
        output_bam=output_bam,
        output_fastq=output_fastq,
        threads=threads,
        coassemble=coassemble,
        filtered=output_filtered,
        log=log,
    )
