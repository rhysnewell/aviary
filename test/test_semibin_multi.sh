#!/usr/bin/env bash
# Quick test of --semibin-mode multi using SRR13153254 (smallest sample ~800M)
# with 3 cobinning samples for speed. Run from the aviary repo root.
# Usage: bash test/test_semibin_multi.sh

set -euo pipefail

DATA="/work/microbiome/aviary_module_benchmarking/data/test_data"
OUTPUT="test_semibin_multi/aviary_out"


pixi run --frozen --manifest-path aviary/pixi.toml -e dev \
    aviary recover \
    -o "${OUTPUT}" \
    -1 "${DATA}/SRR13153254_1.fastq.gz" \
       "${DATA}/ERR12120022_1.fastq.gz" \
       "${DATA}/ERR12120023_1.fastq.gz" \
    -2 "${DATA}/SRR13153254_2.fastq.gz" \
       "${DATA}/ERR12120022_2.fastq.gz" \
       "${DATA}/ERR12120023_2.fastq.gz" \
    --semibin-mode multi \
    --use-megahit \
    --request-gpu --strict \
    -n 32 -t 32 --local-cores 1 \
    -m 280 --coassemble no \
    --snakemake-profile aqua --cluster-retries 3

echo "Done. Check ${OUTPUT}/bins/bin_info.tsv for results."
