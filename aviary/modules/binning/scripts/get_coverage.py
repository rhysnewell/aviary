#!/usr/bin/env python3

import argparse
from subprocess import run, STDOUT
import os

def path_exists(file_path: str) -> bool:
    return os.path.exists(file_path) and os.path.getsize(file_path) > 0

def get_coverage(
    long_reads,
    short_reads_1,
    short_reads_2,
    long_read_type: str,
    input_fasta: str,
    bam_cache: str,
    working_dir: str,
    coverm_output: str,
    maxbin_output: str,
    tmpdir: str,
    threads: int,
    log: str,
):
    if tmpdir: os.environ["TMPDIR"] = tmpdir
    os.makedirs(working_dir, exist_ok=True)

    if long_reads != "none" and not path_exists(f"{working_dir}/long_cov.tsv"):
        if long_read_type in ["ont", "ont_hq"]:
            coverm_cmd = f"coverm contig -t {threads} -r {input_fasta} --single {' '.join(long_reads)} -p minimap2-ont -m length trimmed_mean variance --bam-file-cache-directory {bam_cache} --discard-unmapped --min-read-percent-identity 0.85 --output-file {working_dir}/long_cov.tsv".split()
            with open(log, "a") as logf:
                run(coverm_cmd, stdout=logf, stderr=STDOUT, check=True)

        elif long_read_type in ["rs", "sq", "ccs", "hifi"]:
            coverm_cmd = f"coverm contig -t {threads} -r {input_fasta} --single {' '.join(long_reads)} -p minimap2-pb -m length trimmed_mean variance --bam-file-cache-directory {bam_cache} --discard-unmapped --min-read-percent-identity 0.9 --output-file {working_dir}/long_cov.tsv".split()

            with open(log, "a") as logf:
                run(coverm_cmd, stdout=logf, stderr=STDOUT, check=True)

        else:
            coverm_cmd = f"coverm contig -t {threads} -r {input_fasta} --single {' '.join(long_reads)} -p minimap2-ont -m length trimmed_mean variance --bam-file-cache-directory {bam_cache} --discard-unmapped --min-read-percent-identity 0.85 --output-file {working_dir}/long_cov.tsv".split()

            with open(log, "a") as logf:
                run(coverm_cmd, stdout=logf, stderr=STDOUT, check=True)

    if short_reads_2 != 'none' and not path_exists(f"{working_dir}/short_cov.tsv"):
        coverm_cmd = f"coverm contig -t {threads} -r {input_fasta} -1 {' '.join(short_reads_1)} -2 {' '.join(short_reads_2)} -m metabat --bam-file-cache-directory {bam_cache} --discard-unmapped --output-file {working_dir}/short_cov.tsv".split()
        with open(log, "a") as logf:
                run(coverm_cmd, stdout=logf, stderr=STDOUT, check=True)

    elif short_reads_1  != 'none' and not path_exists(f"{working_dir}/short_cov.tsv"):
        coverm_cmd = f"coverm contig -t {threads} -r {input_fasta} --interleaved {' '.join(short_reads_1)} -m metabat --bam-file-cache-directory {bam_cache} --discard-unmapped --output-file {working_dir}/short_cov.tsv".split()

        with open(log, "a") as logf:
                run(coverm_cmd, stdout=logf, stderr=STDOUT, check=True)

    # Concatenate the two coverage files if both long and short exist
    if long_reads != "none" and (short_reads_1 != "none"):
        short_count = len(short_reads_1)
        long_count = len(long_reads)
        with open(coverm_output, 'w') as file3:
            with open(f'{working_dir}/short_cov.tsv', 'r') as file1:
                with open(f'{working_dir}/long_cov.tsv', 'r') as file2:
                    for idx, (line1, line2) in enumerate(zip(file1, file2)):
                        long_values = line2.strip().split("\t")[1::]
                        del long_values[0::3] # delete length value
                        short_values = line1.strip().split("\t")
                        if idx != 0:
                            long_cov = sum([float(x) for x in long_values[0::2]])
                            short_cov = sum([float(x) for x in short_values[3::2]])
                            tot_depth = (long_cov + short_cov) / (short_count + long_count)
                            line = (
                                "{contig_name}\t"
                                "{contig_length}\t"
                                "{total_avg_depth}\t"
                                "{short}\t"
                                "{long}"
                            ).format(
                                contig_name = short_values[0],
                                contig_length = short_values[1],
                                total_avg_depth = tot_depth,
                                short = "\t".join(short_values[3::]),
                                long = "\t".join(long_values)
                            )
                            print(line, file=file3)
                        else:
                            line = (
                                "{short}\t{long}"
                            ).format(
                                short = "\t".join(short_values),
                                long = "\t".join(long_values)
                            )
                            print(line, file=file3)
                            
    elif long_reads != "none":
        with open(coverm_output, 'w') as file3:
            with open(f'{working_dir}/long_cov.tsv', 'r') as file1:
                for idx, line1 in enumerate(file1):
                    long_values = line1.strip().split("\t")
                    del long_values[4::3] # delete extra length values
                    long_count = len(long_values[2::2])
                    if idx != 0:
                        long_cov = sum([float(x) for x in long_values[2::2]])
                        tot_depth = long_cov / long_count
                        line = (
                            "{contig_name}\t"
                            "{contig_length}\t"
                            "{total_avg_depth}\t"
                            "{long}"
                        ).format(
                            contig_name = long_values[0],
                            contig_length = long_values[1],
                            total_avg_depth = tot_depth,
                            long = "\t".join(long_values[2::])
                        )
                        print(line, file=file3)
                    else:
                        line = (
                            "contigName\tcontigLen\ttotalAvgDepth\t{samples}"
                        ).format(
                            samples = "\t".join(long_values[2::])
                        )
                        print(line, file=file3)
                            
    elif short_reads_1 != "none":  # rename short reads cov if only they exist
        os.rename(f"{working_dir}/short_cov.tsv", coverm_output)

    if long_reads != "none":
        with open(f'{working_dir}/long.cov', 'w') as file3:
            with open(f'{working_dir}/long_cov.tsv', 'r') as file:
                for idx, line1 in enumerate(file):
                    long_values = line1.strip().split("\t")
                    del long_values[4::3] # delete extra length values
                    long_count = len(long_values[2::2])
                    if idx != 0:
                        long_cov = sum([float(x) for x in long_values[2::2]])
                        tot_depth = long_cov / long_count
                        line = (
                            "{contig_name}\t"
                            "{contig_length}\t"
                            "{total_avg_depth}\t"
                            "{long}"
                        ).format(
                            contig_name = long_values[0],
                            contig_length = long_values[1],
                            total_avg_depth = tot_depth,
                            long = "\t".join(long_values[2::])
                        )
                        print(line, file=file3)
                    else:
                        line = (
                            "contigName\tcontigLen\ttotalAvgDepth\t{samples}"
                        ).format(
                            samples = "\t".join(long_values[2::])
                        )
                        print(line, file=file3)
            # os.remove(f"{working_dir}/long_cov.tsv")

    try:
        os.makedirs(f"{working_dir}/maxbin_cov/")
    except OSError:
        pass
    with open(coverm_output) as f, open(maxbin_output, "w") as o:
        cov_list = []
        contig_list = []
        f.readline()
        for line in f:
            contig = line.split()[0]
            contig_list.append(contig)
            depths = line.split()[3:]
            cov_list.append([])
            for i in range(0, len(depths), 2):
                cov_list[-1].append(depths[i])
        for i in range(len(cov_list[0])):
            with open(f"{working_dir}/maxbin_cov/p%d.cov" % i, 'w') as oo:
                for j in range(len(contig_list)):
                    oo.write(contig_list[j] + '\t' + cov_list[j][i] + '\n')
            o.write(f"{working_dir}/maxbin_cov/p%d.cov\n" % i)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate coverage for contigs.")
    parser.add_argument("--long-reads", type=str, nargs='*', required=True, help="Path to long reads.")
    parser.add_argument("--short-reads-1", type=str, nargs='*', required=True, help="Path to first set of short reads.")
    parser.add_argument("--short-reads-2", type=str, nargs='*', required=True, help="Path to second set of short reads.")
    parser.add_argument("--long-read-type", type=str, required=True, help="Type of long reads (e.g., ont, rs, etc.).")
    parser.add_argument("--input-fasta", type=str, required=True, help="Path to input FASTA file.")
    parser.add_argument("--tmpdir", type=str, help="Temporary directory.")
    parser.add_argument("--threads", type=int, required=True, help="Number of threads to use.")
    parser.add_argument("--log", type=str, required=True, help="Path to log file.")
    parser.add_argument("--bam-cache", type=str, required=True, help="Path to BAM cache directory.")
    parser.add_argument("--working-dir", type=str, required=True, help="Path to working directory.")
    parser.add_argument("--coverm-output", type=str, required=True, help="Path to CoverM output file for metabat coverage.")
    parser.add_argument("--maxbin-output", type=str, required=True, help="Path to MaxBin coverage output file.")

    args = parser.parse_args()

    get_coverage(
        long_reads="none" if args.long_reads == ["none"] or args.long_reads == [] else args.long_reads,
        short_reads_1="none" if args.short_reads_1 == ["none"] or args.short_reads_1 == [] else args.short_reads_1,
        short_reads_2="none" if args.short_reads_2 == ["none"] or args.short_reads_2 == [] else args.short_reads_2,
        long_read_type=args.long_read_type,
        input_fasta=args.input_fasta,
        bam_cache=args.bam_cache,
        working_dir=args.working_dir,
        coverm_output=args.coverm_output,
        maxbin_output=args.maxbin_output,
        tmpdir=args.tmpdir,
        threads=args.threads,
        log=args.log,
    )
