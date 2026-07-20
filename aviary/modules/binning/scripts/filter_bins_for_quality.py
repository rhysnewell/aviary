#!/usr/bin/env python3

import argparse
import os
import shutil

import polars as pl


def filter_bins(checkm_path, bins_dir, output_dir, min_completeness, max_contamination):
    checkm_output = pl.read_csv(checkm_path, separator="\t")
    bin_column = checkm_output.columns[0]

    filtered = (
        checkm_output.filter(
            (pl.col("Completeness") >= min_completeness)
            & (pl.col("Contamination") <= max_contamination)
        )
        .get_column(bin_column)
        .to_list()
    )

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    for bin_name in filtered:
        source = os.path.join(bins_dir, f"{bin_name}.fna")
        destination = os.path.join(output_dir, f"{bin_name}.fna")
        if os.path.exists(source):
            os.symlink(os.path.abspath(source), destination)

    with open(os.path.join(output_dir, "filtered_bins.txt"), "w") as handle:
        for bin_name in filtered:
            handle.write(f"{bin_name}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter bins for quality using CheckM2 completeness/contamination thresholds."
    )
    parser.add_argument("--checkm", required=True, help="Path to CheckM2 quality report TSV")
    parser.add_argument("--bins-dir", required=True, help="Directory containing final bin .fna files")
    parser.add_argument("--output-dir", required=True, help="Directory to write filtered bins")
    parser.add_argument(
        "--min-completeness", type=float, default=50.0, help="Minimum completeness threshold"
    )
    parser.add_argument(
        "--max-contamination", type=float, default=5.0, help="Maximum contamination threshold"
    )

    args = parser.parse_args()
    filter_bins(
        args.checkm,
        args.bins_dir,
        args.output_dir,
        args.min_completeness,
        args.max_contamination,
    )
