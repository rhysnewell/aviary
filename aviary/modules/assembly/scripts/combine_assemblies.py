import subprocess
import os
import logging
import shutil


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
            subprocess.run(f"cat {flye_fasta} {short_fasta}", stdout=output_fasta)

        logging.info("Flye and metaSPAdes/Unicycler assemblies combined...")


if __name__ == "__main__":
    try:
        flye_fasta = snakemake.input.flye_fasta
    except AttributeError:
        flye_fasta = None

    try:
        short_fasta = snakemake.input.short_fasta
    except AttributeError:
        short_fasta = None

    output_fasta = snakemake.output.output_fasta

    combine_assemblies(
        flye_fasta=flye_fasta,
        short_fasta=short_fasta,
        output_fasta=output_fasta
    )