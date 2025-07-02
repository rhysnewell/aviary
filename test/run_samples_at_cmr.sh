#!/bin/bash

set -euo pipefail

# Run from base dir of repo
# $ bash test/run_samples_at_cmr.sh
#
# For now, need to make sure manually that all the mqsub'd jobs finish without error at the end.
# Randomly selected via Sandpiper
# Gut: SRR29134037, SRR6028613, SRR17283995
# Soil: SRR28059992, SRR30718229, SRR26545921
# Ocean: ERR599108, SRR13153254, ERR12716859

MQSUB_OUT="mqsub_sample.out"
DATA_DIR="/work/microbiome/aviary_module_benchmarking/data/test_data"
if [[ -f "$MQSUB_OUT" ]]; then
    echo "Removing old mqsub output file: $MQSUB_OUT"
    rm "$MQSUB_OUT"
fi

DATE=$(date +%Y%m%d)
OUTPUT_DIR="/work/microbiome/aviary_module_benchmarking/results/test_data/${DATE}"
mkdir -p "$OUTPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

for sample in \
    "SRR29134037" \
    "SRR6028613" \
    "SRR17283995" \
    "SRR28059992" \
    "SRR30718229" \
    "SRR26545921" \
    "ERR599108" \
    "SRR13153254" \
    "ERR12716859"
do
    echo "Processing sample $sample"
    mqsub -t 32 -m 256 --segregated-log-files --bg -- pixi run --frozen --manifest-path aviary/pixi.toml \
        aviary recover \
        -o "$OUTPUT_DIR/$sample" \
        -1 "$DATA_DIR/${sample}_1.fastq.gz" \
        -2 "$DATA_DIR/${sample}_2.fastq.gz" \
        --binning-only \
        -t 32 \
        -m 256 \
        '&>' "$OUTPUT_DIR/$sample.log" &>> $MQSUB_OUT
done

echo
echo "Need to mqwait for these jobs to finish and check they succeeded. TODO: Make this automatic."
cat $MQSUB_OUT
