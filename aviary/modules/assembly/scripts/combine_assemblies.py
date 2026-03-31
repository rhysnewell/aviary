#!/usr/bin/env python3

import subprocess
import os
import logging
import shutil
import argparse

def combine_assemblies(flye_fasta, short_fasta, output_fasta: str):
    """
    Combines Flye and metaSPAdes/Unicycler assemblies if both are present.
    Otherwise, treats Flye assembly as final assembly.
    """
    print(flye_fasta, short_fasta, output_fasta)

    if os.path.exists(output_fasta):
        # remove output_fasta if it already exists
        os.remove(output_fasta)

    if flye_fasta is None:
        shutil.copyfile(short_fasta, output_fasta)
        logging.info("Treating Spades assembly as final assembly...")
    elif short_fasta is None:
        shutil.copyfile(flye_fasta, output_fasta)
        logging.info("Treating Flye assembly as final assembly...")
    else:
        with open(output_fasta, 'w') as output:
            subprocess.run(f"cat {flye_fasta} {short_fasta}", shell=True, stdout=output)

        logging.info("Flye and metaSPAdes/Unicycler assemblies combined...")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine Flye and metaSPAdes/Unicycler assemblies.")
    parser.add_argument(
        "--flye-fasta",
        type=str,
        required=False,
        help="Path to the Flye assembly FASTA file."
    )
    parser.add_argument(
        "--short-fasta",
        type=str,
        required=False,
        help="Path to the metaSPAdes/Unicycler assembly FASTA file."
    )
    parser.add_argument(
        "--output-fasta",
        type=str,
        required=True,
        help="Path to the output combined FASTA file."
    )

    args = parser.parse_args()

    combine_assemblies(
        flye_fasta=args.flye_fasta,
        short_fasta=args.short_fasta,
        output_fasta=args.output_fasta
    )