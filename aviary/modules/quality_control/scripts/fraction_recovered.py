#!/usr/bin/env python3
from subprocess import run, STDOUT
import os
import sys
import argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
from aviary import resolve_mapper_model, short_read_mapper_extra_params


def get_fraction_recovered(
    long_reads,
    short_reads_1,
    short_reads_2,
    input_fasta: str,
    long_read_type: str,
    long_read_mapper: str,
    short_read_mapper: str,
    threads: int,
    log: str,
    long_read_mapper_model: str = None,
    minibwa_params: str = None,
    bwa_params: str = None,
    strobealign_params: str = None,
    minimap2_params: str = None,
    rammap_params: str = None,
):
    if long_reads != "none" and not os.path.exists("data/long_cov.tsv"):
        # Identity cutoff mirrors get_coverage.py/get_abundances.py: 0.85 for
        # ONT-family reads, 0.9 for PacBio-family reads.
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
            f"coverm genome -t {threads} --single-genome -r {input_fasta} --single {' '.join(long_reads)} "
            f"-p {mapper_p} --min-read-percent-identity {min_identity} "
            f"-o www/fraction_recovered/long_fraction_recovered.tsv"
        ).split() + extra_params

        with open(log, "a") as logf:
            print("Running command:", " ".join(coverm_cmd), file=logf)
            run(coverm_cmd, stdout=logf, stderr=STDOUT, check=True)

    short_extra_params = short_read_mapper_extra_params(
        short_read_mapper, bwa_params, strobealign_params, minimap2_params, rammap_params)

    if short_reads_2 != 'none' and not os.path.exists("data/short_cov.tsv"):
        coverm_cmd = (f"coverm genome -t {threads} --single-genome -r {input_fasta} -1 {' '.join(short_reads_1)} -2 {' '.join(short_reads_2)} -p {short_read_mapper} -o www/fraction_recovered/short_fraction_recovered.tsv".split()
                       + short_extra_params)

        with open(log, "a") as logf:
            print("Running command:", " ".join(coverm_cmd), file=logf)
            run(coverm_cmd, stdout=logf, stderr=STDOUT, check=True)

    elif short_reads_1 != 'none' and not os.path.exists("data/short_cov.tsv"):
        coverm_cmd = (f"coverm genome -t {threads} --single-genome -r {input_fasta} --interleaved {' '.join(short_reads_1)} -p {short_read_mapper} -o www/fraction_recovered/short_fraction_recovered.tsv".split()
                       + short_extra_params)

        with open(log, "a") as logf:
            print("Running command:", " ".join(coverm_cmd), file=logf)
            run(coverm_cmd, stdout=logf, stderr=STDOUT, check=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate fraction recovered from sequencing reads')
    parser.add_argument('--input-fasta', required=True, help='Input FASTA file')
    parser.add_argument('--long-reads', nargs='+', default=['none'], help='Long read files')
    parser.add_argument('--short-reads-1', nargs='+', default=['none'], help='Short reads (first pair or interleaved)')
    parser.add_argument('--short-reads-2', nargs='+', default=['none'], help='Short reads (second pair)')
    parser.add_argument('--short-read-mapper', default='strobealign',
                        help='CoverM -p value for short reads')
    parser.add_argument('--long-read-mapper', default='rammap',
                        help='Aligner family for long reads (rammap, minimap2 or minibwa)')
    parser.add_argument('--long-read-mapper-model', default=None,
                        help="Explicit CoverM preset for --long-read-mapper, overriding the "
                             "--long-read-type default. Ignored for minibwa.")
    parser.add_argument('--minibwa-params', default=None,
                        help="Raw passthrough params for minibwa (e.g. '-x lr'); only used "
                             "when --long-read-mapper is minibwa.")
    parser.add_argument('--bwa-params', default=None,
                        help="Raw CoverM --bwa-params passthrough; only used when "
                             "--short-read-mapper is bwa-mem or bwa-mem2.")
    parser.add_argument('--strobealign-params', default=None,
                        help="Raw CoverM --strobealign-params passthrough; only used when "
                             "--short-read-mapper is strobealign.")
    parser.add_argument('--minimap2-params', default=None,
                        help="Raw CoverM --minimap2-params passthrough; used when either "
                             "--short-read-mapper or --long-read-mapper is minimap2.")
    parser.add_argument('--rammap-params', default=None,
                        help="Raw CoverM --rammap-params passthrough; used when either "
                             "--short-read-mapper or --long-read-mapper is rammap.")
    parser.add_argument('--long-read-type', default='ont', help='Long read type (ont, ont_hq, rs, sq, ccs, hifi)')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    parser.add_argument('--log', default='get_fraction_recovered.log', help='Log file')
    
    args = parser.parse_args()
    
    os.makedirs('www/fraction_recovered', exist_ok=True)

    with open(args.log, "w") as logf: 
        pass

    get_fraction_recovered(
        long_reads="none" if args.long_reads in (["none"], [], "none") else args.long_reads,
        short_reads_1="none" if args.short_reads_1 in (["none"], [], "none") else args.short_reads_1,
        short_reads_2="none" if args.short_reads_2 in (["none"], [], "none") else args.short_reads_2,
        input_fasta=args.input_fasta,
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
        log=args.log,
    )