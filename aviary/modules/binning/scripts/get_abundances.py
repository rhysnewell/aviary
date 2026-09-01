#!/usr/bin/env python

import argparse
from subprocess import run, STDOUT
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
from aviary import resolve_mapper_model, LONG_READ_MODELS, short_read_mapper_extra_params

# Long-read aligners coverm accepts. rammap is a minimap2-compatible Rust
# implementation and takes the same -x presets; --long-read-mapper selects
# which family is used. minibwa has no preset suffix of its own here -- CoverM
# just runs `minibwa map` with default parameters (its own adaptive "adap"
# mode, which auto-detects short vs long reads via a length threshold rather
# than taking an explicit preset), so it is passed through as a bare
# "minibwa" -p value rather than a "<family>-<model>" one. Verified against
# the installed minibwa 0.6 binary: adap and an explicit -x lr produce
# equivalent alignments on long reads (same reads mapped, same mean MAPQ), so
# this is a real, working long-read mapper choice, not a fallback that
# silently produces bad numbers -- unlike an unlisted mapper name, which
# CoverM would default to strobealign (short-read only) for.
#
# Built from LONG_READ_MODELS (aviary/__init__.py), the single source of
# truth --long-read-mapper-model is validated against in processor.py, rather
# than a second hand-maintained tuple: a previous hardcoded version of this
# list omitted the "no-preset" model, so --long-read-mapper-model no-preset
# passed CLI validation but was then rejected here as "not an explicit
# long-read mapper", even though resolve_mapper_model() had correctly
# resolved it to e.g. "minimap2-no-preset".
LONG_READ_MAPPERS = ("minibwa",) + tuple(
    f"{family}-{model}"
    for family in ("minimap2", "rammap")
    for model in LONG_READ_MODELS
)


def run_coverm(
    reads: str,
    minimap2_type: str,
    output_file: str,
    read_type: str,
    threads: int,
    strain_analysis: bool,
    output_dir: str,
    log: str,
    bins_dir: str,
    long_reads: bool = False,
    extra_params: list = None,
):
    # CoverM >=0.7.0 defaults to strobealign, which is SHORT READ ONLY. Long
    # reads must therefore always name a minimap2 mapper explicitly -- omitting
    # -p would map them with strobealign, which returns a well-formed coverage
    # table of near-zero depths rather than failing, so the mistake would show
    # up as wrong abundances rather than an error.
    if long_reads and minimap2_type not in LONG_READ_MAPPERS:
        raise ValueError(
            "Long reads require an explicit long-read mapper "
            f"({' or '.join(LONG_READ_MAPPERS)}), got {minimap2_type!r}."
        )

    strain_analysis_flag = f"--bam-file-cache-directory {output_dir} --discard-unmapped" if strain_analysis else ""

    # extra_params (e.g. --minibwa-params '-x lr') is appended as separate argv
    # tokens rather than folded into the f-string below, since run() is called
    # without a shell and a naive .split() on an embedded quoted value would
    # mis-tokenize it.
    coverm_cmd = f"coverm genome -t {threads} {strain_analysis_flag} -d {bins_dir} -m relative_abundance covered_fraction {read_type} {reads} -p {minimap2_type} --output-file {output_file} --min-covered-fraction 0.0 -x fna".split()
    if extra_params:
        coverm_cmd += extra_params

    with open(log, "a") as logf:
        # Logged so the mapper actually used is recoverable from the run, the
        # same as get_coverage.py does.
        print("Running command:", " ".join(coverm_cmd), file=logf)
        # check=True because CoverM creates its --output-file before it fails:
        # a mid-run crash leaves a 0-byte abundance table behind, which
        # satisfies snakemake's missing-output check, so the rule is marked
        # successful and finalise_stats folds an empty table into bin_info.tsv.
        # Verified against coverm 0.8: exit 1, output file present, 0 bytes.
        run(coverm_cmd, stdout=logf, stderr=STDOUT, check=True)

def get_abundances(
    long_reads,
    short_reads_1,
    short_reads_2,
    long_read_type: str,
    long_read_mapper: str,
    short_read_mapper: str,
    threads: int,
    strain_analysis: bool,
    log: str,
    bins_dir: str,
    long_read_mapper_model: str = None,
    minibwa_params: str = None,
    bwa_params: str = None,
    strobealign_params: str = None,
    minimap2_params: str = None,
    rammap_params: str = None,
):
    if long_reads != "none":
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
        run_coverm(
            reads=" ".join(long_reads),
            minimap2_type=mapper_p,
            long_reads=True,
            output_file="data/long_abundances.tsv",
            read_type="--single",
            threads=threads,
            strain_analysis=strain_analysis,
            output_dir="data/reads_mapped_to_mags/long/",
            log=log,
            bins_dir=bins_dir,
            extra_params=extra_params,
        )

    short_extra_params = short_read_mapper_extra_params(
        short_read_mapper, bwa_params, strobealign_params, minimap2_params, rammap_params)

    if short_reads_2 != 'none':
        reads_str = []
        for r1, r2 in zip(short_reads_1, short_reads_2):
            reads_str.append(f"{r1} {r2}")
        reads_str = " ".join(reads_str)

        run_coverm(
            reads=reads_str,
            minimap2_type=short_read_mapper,
            output_file="data/short_abundances.tsv",
            read_type="--coupled",
            threads=threads,
            strain_analysis=strain_analysis,
            output_dir="data/reads_mapped_to_mags/short/",
            log=log,
            bins_dir=bins_dir,
            extra_params=short_extra_params,
        )

    elif short_reads_1 != 'none':
        run_coverm(
            reads=" ".join(short_reads_1),
            minimap2_type=short_read_mapper,
            output_file="data/short_abundances.tsv",
            read_type="--interleaved",
            threads=threads,
            strain_analysis=strain_analysis,
            output_dir="data/reads_mapped_to_mags/short/",
            log=log,
            bins_dir=bins_dir,
            extra_params=short_extra_params,
        )

    # Concatenate the two coverage files if both long and short exist
    if long_reads != "none" and short_reads_1 != "none":
        with open('data/coverm_abundances.tsv', 'w') as file3:
            with open('data/short_abundances.tsv', 'r') as file1:
                with open('data/long_abundances.tsv', 'r') as file2:
                    for line1, line2 in zip(file1, file2):
                        long_cov_line = "\t".join([l.strip() for l in line2.strip().split('\t')[1::]])
                        print(line1.strip(), "\t", long_cov_line, file=file3)
    elif long_reads != "none":  # rename long reads cov if only it exists
        os.rename("data/long_abundances.tsv", "data/coverm_abundances.tsv")
    elif short_reads_1 != "none":  # rename short reads cov if only they exist
        os.rename("data/short_abundances.tsv", "data/coverm_abundances.tsv")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get abundances using CoverM.")
    parser.add_argument("--long-reads", nargs="+", default="none", help="Paths to long reads files.")
    parser.add_argument("--short-reads-1", nargs="+", default="none", help="Paths to first set of short reads files.")
    parser.add_argument("--short-reads-2", nargs="+", default="none", help="Paths to second set of short reads files.")
    parser.add_argument("--short-read-mapper", type=str, default="strobealign",
                        help="CoverM -p value for short reads.")
    parser.add_argument("--long-read-mapper", type=str, default="rammap",
                        help="Aligner family for long reads (rammap, minimap2 or minibwa).")
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
                             "--short-read-mapper is strobealign.")
    parser.add_argument("--minimap2-params", type=str, default=None,
                        help="Raw CoverM --minimap2-params passthrough; used when either "
                             "--short-read-mapper or --long-read-mapper is minimap2.")
    parser.add_argument("--rammap-params", type=str, default=None,
                        help="Raw CoverM --rammap-params passthrough; used when either "
                             "--short-read-mapper or --long-read-mapper is rammap.")
    parser.add_argument("--long-read-type", type=str, required=True, help="Type of long reads (e.g., ont, rs, etc.).")
    parser.add_argument("--threads", type=int, required=True, help="Number of threads to use.")
    parser.add_argument("--strain-analysis", type=lambda x: x.lower() == 'true', nargs='?', const=True, default=False, help="Enable strain analysis (True/False).")
    parser.add_argument("--log", type=str, required=True, help="Path to log file.")
    parser.add_argument("--bins-dir", type=str, required=True, help="Directory containing bins to quantify.")

    args = parser.parse_args()

    with open(args.log, "w") as logf:
        pass

    get_abundances(
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
        threads=args.threads,
        strain_analysis=args.strain_analysis,
        log=args.log,
        bins_dir=args.bins_dir,
    )
