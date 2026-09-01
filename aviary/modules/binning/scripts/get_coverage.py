#!/usr/bin/env python3

import argparse
from subprocess import run, STDOUT
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
from aviary import resolve_mapper_model, short_read_mapper_extra_params

def path_exists(file_path: str) -> bool:
    return os.path.exists(file_path) and os.path.getsize(file_path) > 0


def _fasta_lengths(fasta_path: str) -> dict:
    lengths = {}
    name = None
    length = 0
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line)
        if name is not None:
            lengths[name] = length
    return lengths


def _reshape_aemb_to_metabat(aemb_tsv: str, fasta_path: str, output_path: str) -> None:
    """`coverm contig -m strobealign-aemb` emits a two-column TSV (contig,
    mean coverage) -- it shells out to `strobealign --aemb` directly and
    skips CoverM's normal BAM/pileup pipeline, so unlike every other -m value
    it has no contig length or variance column. CoverM 0.8.0 additionally
    prepends a header row ("Contig\t<sample> Strobealign aemb"); rows whose
    first column isn't a real contig (i.e. that header) are skipped rather
    than looked up in `lengths`, so this doesn't depend on the header's exact
    wording and stays correct if a future CoverM version drops it again.
    Everything downstream of short_cov.tsv in this file (the long+short
    concat logic and the maxbin writer) assumes the metabat depth-file shape
    (contigName, contigLen, totalAvgDepth, <sample>.bam, <sample>.bam-var,
    ...), and so do the binners that read data/coverm.cov directly. Reshape
    aemb's output into that shape so strobealign-aemb is a drop-in short-read
    coverage source rather than needing its own binner wiring. Variance is
    not computed by aemb, so it is filled with 0.
    """
    lengths = _fasta_lengths(fasta_path)
    with open(aemb_tsv) as f, open(output_path, "w") as out:
        print("contigName\tcontigLen\ttotalAvgDepth\tshort_reads.bam\tshort_reads.bam-var", file=out)
        for line in f:
            contig, depth = line.rstrip("\n").split("\t")
            if contig not in lengths:
                continue
            print(f"{contig}\t{lengths[contig]}\t{depth}\t{depth}\t0", file=out)

def get_coverage(
    long_reads,
    short_reads_1,
    short_reads_2,
    long_read_type: str,
    long_read_mapper: str,
    short_read_mapper: str,
    input_fasta: str,
    bam_cache: str,
    working_dir: str,
    coverm_output: str,
    maxbin_output: str,
    tmpdir: str,
    threads: int,
    log: str,
    long_read_mapper_model: str = None,
    minibwa_params: str = None,
    bwa_params: str = None,
    strobealign_params: str = None,
    minimap2_params: str = None,
    rammap_params: str = None,
):
    if tmpdir: os.environ["TMPDIR"] = tmpdir
    os.makedirs(working_dir, exist_ok=True)
    os.makedirs(bam_cache, exist_ok=True)

    if long_reads != "none" and not path_exists(f"{working_dir}/long_cov.tsv"):
        # minibwa identity cutoff mirrors get_abundances.py/fraction_recovered.py:
        # 0.85 for ONT-family reads, 0.9 for PacBio-family reads.
        min_identity = 0.85 if long_read_type in ["ont", "ont_hq"] else 0.9
        if long_read_mapper == "minibwa":
            mapper_p = "minibwa"
            extra_params = ["--minibwa-params", minibwa_params] if minibwa_params else []
        else:
            mapper_p = resolve_mapper_model(long_read_mapper, long_read_mapper_model, long_read_type)
            if long_read_mapper == "minimap2" and minimap2_params:
                extra_params = ["--minimap2-params", minimap2_params]
            elif long_read_mapper == "rammap" and rammap_params:
                extra_params = ["--rammap-params", rammap_params]
            else:
                extra_params = []
        # extra_params is appended as separate argv tokens rather than folded
        # into the f-string below, since run() is called without a shell and a
        # naive .split() on an embedded quoted value would mis-tokenize it.
        coverm_cmd = (
            f"coverm contig -t {threads} -r {input_fasta} --single {' '.join(long_reads)} "
            f"-p {mapper_p} -m length trimmed_mean variance "
            f"--bam-file-cache-directory {bam_cache} --discard-unmapped "
            f"--min-read-percent-identity {min_identity} --output-file {working_dir}/long_cov.tsv"
        ).split() + extra_params
        with open(log, "a") as logf:
            # log the coverm command
            print("Running command:", " ".join(coverm_cmd), file=logf)
            run(coverm_cmd, stdout=logf, stderr=STDOUT, check=True)

    is_aemb = short_read_mapper == "strobealign-aemb"
    aemb_raw_output = f"{working_dir}/short_cov.aemb_raw.tsv"
    short_extra_params = short_read_mapper_extra_params(
        short_read_mapper, bwa_params, strobealign_params, minimap2_params, rammap_params)

    if short_reads_2 != 'none' and not path_exists(f"{working_dir}/short_cov.tsv"):
        if is_aemb:
            # No -p (strobealign-aemb ignores it and always runs strobealign
            # itself -- see aviary/__init__.py's SHORT_READ_MAPPERS comment),
            # no --bam-file-cache-directory/--discard-unmapped (no BAM is
            # produced), and it must be the only -m value given.
            coverm_cmd = f"coverm contig -t {threads} -r {input_fasta} -1 {' '.join(short_reads_1)} -2 {' '.join(short_reads_2)} -m strobealign-aemb --output-file {aemb_raw_output}".split()
        else:
            coverm_cmd = f"coverm contig -t {threads} -r {input_fasta} -1 {' '.join(short_reads_1)} -2 {' '.join(short_reads_2)} -p {short_read_mapper} -m metabat --bam-file-cache-directory {bam_cache} --discard-unmapped --output-file {working_dir}/short_cov.tsv".split() + short_extra_params
        with open(log, "a") as logf:
            print("Running command:", " ".join(coverm_cmd), file=logf)
            run(coverm_cmd, stdout=logf, stderr=STDOUT, check=True)
        if is_aemb:
            _reshape_aemb_to_metabat(aemb_raw_output, input_fasta, f"{working_dir}/short_cov.tsv")

    elif short_reads_1  != 'none' and not path_exists(f"{working_dir}/short_cov.tsv"):
        if is_aemb:
            coverm_cmd = f"coverm contig -t {threads} -r {input_fasta} --interleaved {' '.join(short_reads_1)} -m strobealign-aemb --output-file {aemb_raw_output}".split()
        else:
            coverm_cmd = f"coverm contig -t {threads} -r {input_fasta} --interleaved {' '.join(short_reads_1)} -p {short_read_mapper} -m metabat --bam-file-cache-directory {bam_cache} --discard-unmapped --output-file {working_dir}/short_cov.tsv".split() + short_extra_params

        with open(log, "a") as logf:
            print("Running command:", " ".join(coverm_cmd), file=logf)
            run(coverm_cmd, stdout=logf, stderr=STDOUT, check=True)
        if is_aemb:
            _reshape_aemb_to_metabat(aemb_raw_output, input_fasta, f"{working_dir}/short_cov.tsv")

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
    parser.add_argument("--long-read-mapper", type=str, default="rammap",
                        help="Aligner family for long reads; the -x preset is chosen from --long-read-type.")
    parser.add_argument("--long-read-mapper-model", type=str, default=None,
                        help="Explicit CoverM preset for --long-read-mapper, overriding the "
                             "--long-read-type default. Ignored for minibwa.")
    parser.add_argument("--minibwa-params", type=str, default=None,
                        help="Raw passthrough params for minibwa (e.g. '-x lr'); only used "
                             "when --long-read-mapper is minibwa.")
    parser.add_argument("--bwa-params", type=str, default=None,
                        help="Raw CoverM --bwa-params passthrough; only used when "
                             "--short-read-mapper is bwa-mem or bwa-mem2.")
    parser.add_argument("--strobealign-params", type=str, default=None,
                        help="Raw CoverM --strobealign-params passthrough; only used when "
                             "--short-read-mapper is strobealign (not strobealign-aemb, which "
                             "does not accept it).")
    parser.add_argument("--minimap2-params", type=str, default=None,
                        help="Raw CoverM --minimap2-params passthrough; used when either "
                             "--short-read-mapper or --long-read-mapper is minimap2.")
    parser.add_argument("--rammap-params", type=str, default=None,
                        help="Raw CoverM --rammap-params passthrough; used when either "
                             "--short-read-mapper or --long-read-mapper is rammap.")
    parser.add_argument("--short-read-mapper", type=str, default="strobealign",
                        help="CoverM -p value for short reads, or 'strobealign-aemb' to use "
                             "CoverM's -m strobealign-aemb direct abundance estimator instead "
                             "of a normal alignment (coverm-contig only; reshaped into the "
                             "metabat depth-file shape here).")
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
        long_read_mapper=args.long_read_mapper,
        long_read_mapper_model=None if args.long_read_mapper_model in (None, "none") else args.long_read_mapper_model,
        minibwa_params=None if args.minibwa_params in (None, "none") else args.minibwa_params,
        bwa_params=None if args.bwa_params in (None, "none") else args.bwa_params,
        strobealign_params=None if args.strobealign_params in (None, "none") else args.strobealign_params,
        minimap2_params=None if args.minimap2_params in (None, "none") else args.minimap2_params,
        rammap_params=None if args.rammap_params in (None, "none") else args.rammap_params,
        short_read_mapper=args.short_read_mapper,
        input_fasta=args.input_fasta,
        bam_cache=args.bam_cache,
        working_dir=args.working_dir,
        coverm_output=args.coverm_output,
        maxbin_output=args.maxbin_output,
        tmpdir=args.tmpdir,
        threads=args.threads,
        log=args.log,
    )
