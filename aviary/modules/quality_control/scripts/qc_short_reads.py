from subprocess import Popen, PIPE, run, STDOUT
from pathlib import Path
import os
from typing import List

def cat_reads(read_path: str, output_path: str, threads: int, log: str):
    """
    :param read_path: path to reads
    :param output_path: path to output
    :param threads: number of threads
    :param log: path to log file
    :return:
    """
    with open(log, "a") as logf:
        with open(output_path, 'a') as out:
            cat_or_zcat = 'zcat' if read_path.endswith('.gz') else 'cat'
            cat_cmd = f'{cat_or_zcat} {read_path}'.split()
            pigz_cmd = f'pigz -p {threads}'.split()

            logf.write(f"Shell style : {' '.join(cat_cmd)} | {' '.join(pigz_cmd)} > {output_path}\n")

            cat_p1 = Popen(cat_cmd, stdout=PIPE, stderr=logf)
            pigz_p2 = Popen(pigz_cmd, stdin=cat_p1.stdout, stdout=out, stderr=logf)

            pigz_p2.wait()
            cat_p1.wait()
            logf.write(f"cat return: {cat_p1.returncode}\n")
            logf.write(f"pigz return: {pigz_p2.returncode}\n")

def combine_reads(
    short_reads_1,
    short_reads_2,
    output_fastq: str,
    coassemble: bool,
    log_file: str
):
    """
    Combine reads before QC
    :param long_reads: list of long reads
    :param output_long_reads: output long reads
    :param coassemble: coassemble or not, if true we will filter all reads into the same file
    :return:
    """
    with open(log_file, 'a') as logf:
        logf.write(f"Combining reads before quality control\n")
        logf.write(f"Coassemble: {coassemble}\n")

        if "none" in short_reads_1 and "none" in short_reads_2:
            logf.write(f"Creating dummy files")
            Path(output_fastq).touch()


        # we've got to concatenate the files together
        if coassemble:
            # reads 1 first
            for reads in short_reads_1:
                if reads == "none":
                    continue
                
                if not os.path.exists(reads):
                    logf.write(f"Short read file {reads} does not exist\n")
                    exit(1)

                cat_reads(reads, output_fastq, threads, log)

            # reads 2 second
            for reads in short_reads_2:
                if reads == "none":
                    continue

                if not os.path.exists(reads):
                    logf.write(f"Short read file {reads} does not exist\n")
                    exit(1)

                cat_reads(reads, output_fastq, threads, log)
        
        else:
            # if we have paired reads, we need to concatenate them together
            if "none" not in short_reads_2 and "none" not in short_reads_1:
                for reads1, reads2 in zip(short_reads_1, short_reads_2):
                    if not os.path.exists(reads1):
                        logf.write(f"Short read file {reads1} does not exist\n")
                        exit(1)

                    if not os.path.exists(reads2):
                        logf.write(f"Short read file {reads2} does not exist\n")
                        exit(1)

                    cat_reads(reads1, output_fastq, threads, log)
                    cat_reads(reads2, output_fastq, threads, log)
                    break
            elif "none" not in short_reads_1:
                # otherwise we just need to symlink the first file
                logf.write(f"Symlinking {short_reads_1[0]} to {output_fastq}\n")
                if not os.path.exists(short_reads_1[0]):
                    logf.write(f"Short read file {short_reads_1[0]} does not exist\n")
                    exit(1)

                os.symlink(short_reads_1[0], output_fastq)
            elif "none" not in short_reads_2:
                # otherwise we just need to symlink the first file
                logf.write(f"Symlinking {short_reads_2[0]} to {output_fastq}\n")
                if not os.path.exists(short_reads_2[0]):
                    logf.write(f"Short read file {short_reads_2[0]} does not exist\n")
                    exit(1)

                os.symlink(short_reads_2[0], output_fastq)
            else:
                logf.write(f"Both reads_1 and reads_2 are None, Error has occured in read concatenation.\n")
                exit(1)



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

        minimap_p1 = Popen(minimap_cmd, stdout=PIPE, stderr=logf)
        samtools_view_p2 = Popen(samtools_view_cmd, stdin=minimap_p1.stdout, stderr=logf)
        samtools_view_p2.wait()

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

            samtools_bam2fq_p1.wait()
            pigz_p2.wait()
            logf.write(f"samtools bam2fq return: {samtools_bam2fq_p1.returncode}\n")
            logf.write(f"pigz return: {pigz_p2.returncode}\n")

def run_fastp(
    reads_1,
    reads_2,
    output_fastq: str,
    threads: int,
    disable_adapter_trimming: bool,
    quality_cutoff: int,
    unqualified_percent_limit: int,
    min_length: int,
    max_length: int,
    extra_fastp_params: str,
    log: str,
) -> str:
    """
    :param input_fastq: input fastq file
    :param output_fastq: output fastq file
    :param threads: number of threads
    :param disable_adapter_trimming: disable adapter trimming
    :param quality_cutoff: quality cutoff
    :param unqualified_percent_limit: unqualified percent limit
    :param min_length: minimum length
    :param max_length: maximum length
    :return:
    """
    with open(log, "a") as logf:
        if reads_1 is None and reads_2 is None:
            logf.write(f"Both reads_1 and reads_2 are None, Error has occured in read concatenation.\n")
            exit(1)

        fastp_cmd_list = ["fastp", "--stdout", "-w", str(threads), "-q", str(quality_cutoff), "-u", str(unqualified_percent_limit), "-l", str(min_length), "--length_limit", str(max_length)]
        if disable_adapter_trimming:
            fastp_cmd_list.append("--disable_adapter_trimming")

        if reads_2 is not None:
            fastp_cmd_list.extend(["-i", reads_1, "-I", reads_2])
        
        if reads_1 is not None and reads_2 is None:
            fastp_cmd_list.extend(["-i", reads_1])
        
        fastp_cmd_list.extend(extra_fastp_params.split())

        pigz_cmd = f"pigz -p {threads}".split()

        logf.write(f"Shell style : {' '.join(fastp_cmd_list)} | {' '.join(pigz_cmd)} > {output_fastq}\n")
        
        with open(output_fastq, 'w') as output_fq:
            fastp_p1 = Popen(fastp_cmd_list, stdout=PIPE, stderr=logf)
            pigz_p2 = Popen(pigz_cmd, stdin=fastp_p1.stdout, stdout=output_fq, stderr=logf)

            pigz_p2.wait()
            fastp_p1.wait()
            logf.write(f"fastp return: {fastp_p1.returncode}\n")
            logf.write(f"pigz return: {pigz_p2.returncode}\n")

    return output_fastq

def filter_illumina_reference(
    short_reads_1,
    short_reads_2,
    disable_adapter_trimming: bool,
    quality_cutoff: int,
    unqualified_percent_limit: int,
    min_length: int,
    max_length: int,
    extra_fastp_params: str,
    reference_filter: List[str],
    output_bam: str,
    output_fastq: str,
    threads: int,
    coassemble: bool,
    filtered: str,
    log: str,
    skip_qc: bool,
):
    
    if skip_qc or ("none" in short_reads_1 and "none" in short_reads_2):
        combine_reads(short_reads_1, short_reads_2, output_fastq, coassemble, log)
        with open(log, 'a') as logf:
            logf.write(f"Skipping quality control\n")
        Path(filtered).touch()
        Path(output_bam).touch()
        return

    # run fastp for quality control of reads
    se1_string = None
    if "none" not in short_reads_1:
        combine_reads(short_reads_1, ["none"], "data/short_reads.pre_qc.1.fastq.gz", coassemble, log)
        se1_string = "data/short_reads.pre_qc.1.fastq.gz"
    
    se2_string = None
    if "none" not in short_reads_2:
        combine_reads(["none"], short_reads_2, "data/short_reads.pre_qc.2.fastq.gz", coassemble, log)
        se2_string = "data/short_reads.pre_qc.2.fastq.gz"
    
    run_fastp(
        reads_1=se1_string,
        reads_2=se2_string,
        output_fastq="data/short_reads.fastq.gz",
        threads=threads,
        disable_adapter_trimming=disable_adapter_trimming,
        quality_cutoff=quality_cutoff,
        unqualified_percent_limit=unqualified_percent_limit,
        min_length=min_length,
        max_length=max_length,
        extra_fastp_params=extra_fastp_params,
        log=log,
    )

    # Move fastp.json and fastp.html to www/
    if os.path.exists("fastp.json") and os.path.exists("fastp.html"):
        # make sure www/ exists
        if not os.path.exists("www/"):
            os.mkdir("www/")
        os.rename("fastp.html", "www/fastp.html")
        os.rename("fastp.json", "www/fastp.json")

    # remove the pre qc files
    if os.path.exists("data/short_reads.pre_qc.1.fastq.gz"):
        os.remove("data/short_reads.pre_qc.1.fastq.gz")
    if os.path.exists("data/short_reads.pre_qc.2.fastq.gz"):
        os.remove("data/short_reads.pre_qc.2.fastq.gz")

    

    reference_filter_file_string = ''
    with open(log, "a") as logf:
        if len(reference_filter) == 0:
            logf.write(f"Not performing reference filtering: {reference_filter}\n")
            Path(filtered).touch()
            Path(output_bam).touch()
            return
        
        if len(reference_filter) > 1:
            with open(f'data/reference_filter.fasta', 'w') as out:
                for reference in reference_filter:
                    # check if file exists
                    if not os.path.exists(reference):
                        logf.write(f"Reference filter file {reference} does not exist\n")
                        exit(1)

                    # concatenate accoutning for gzipped files
                    cat_or_zcat = 'zcat' if reference.endswith('.gz') else 'cat'
                    cat_cmd = f'{cat_or_zcat} {reference}'.split()

                    logf.write(f"Shell style : {' '.join(cat_cmd)} > data/reference_filter.fasta\n")

                    cat_p1 = Popen(cat_cmd, stdout=out, stderr=logf)
                    cat_p1.wait()
                    logf.write(f"cat return: {cat_p1.returncode}\n")

            # gzip the concatenated file
            pigz_cmd = f'pigz -p {threads} data/reference_filter.fasta'.split()
            logf.write(f"Shell style : {' '.join(pigz_cmd)}\n")

            pigz_p1 = Popen(pigz_cmd, stderr=logf)
            pigz_p1.wait()
            logf.write(f"pigz return: {pigz_p1.returncode}\n")

            reference_filter_file_string = f'data/reference_filter.fasta.gz'
        else:
            # make sure file exists
            if not os.path.exists(reference_filter[0]):
                logf.write(f"Reference filter file {reference_filter[0]} does not exist\n")
                exit(1)

            reference_filter_file_string = f'{reference_filter[0]}'


    run_mapping_process(
        reads_string='data/short_reads.fastq.gz',
        input_fasta=reference_filter_file_string,
        output_bam=output_bam,
        output_fastq=output_fastq,
        threads=threads,
        log=log,
    )

    # remove the reference filter file
    if os.path.exists("data/reference_filter.fasta.gz"):
        os.remove("data/reference_filter.fasta.gz")

    Path(filtered).touch()

if __name__ == '__main__':
    short_reads_1 = snakemake.input.short_reads_1
    short_reads_2 = snakemake.input.short_reads_2
    output_bam = snakemake.output.bam
    output_fastq = snakemake.output.fastq
    output_filtered = snakemake.output.filtered

    # fastp parameters
    disable_adapter_trimming = snakemake.params.disable_adapter_trimming
    quality_cutoff = snakemake.params.quality_cutoff
    unqualified_percent_limit = snakemake.params.unqualified_percent_limit
    min_length = snakemake.params.min_length
    max_length = snakemake.params.max_length
    extra_fastp_params = snakemake.params.extra_fastp_params

    coassemble = snakemake.params.coassemble
    reference_filter = snakemake.params.reference_filter
    skip_qc = snakemake.params.skip_qc

    threads = snakemake.threads
    log = snakemake.log[0]

    with open(log, "w") as logf: pass

    filter_illumina_reference(
        short_reads_1=short_reads_1,
        short_reads_2=short_reads_2,
        disable_adapter_trimming=disable_adapter_trimming,
        quality_cutoff=quality_cutoff,
        unqualified_percent_limit=unqualified_percent_limit,
        min_length=min_length,
        max_length=max_length,
        extra_fastp_params=extra_fastp_params,
        reference_filter=reference_filter,
        output_bam=output_bam,
        output_fastq=output_fastq,
        threads=threads,
        coassemble=coassemble,
        filtered=output_filtered,
        log=log,
        skip_qc=skip_qc,
    )
