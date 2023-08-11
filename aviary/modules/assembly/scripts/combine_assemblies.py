import subprocess
import os
import logging


def combine_assemblies(flye_fasta: str, short_fasta: str, input_fasta: str, output_fasta: str):
    """
    Combines Flye and metaSPAdes/Unicycler assemblies if both are present.
    Otherwise, treats Flye assembly as final assembly.
    """
    if not os.path.exists(output_fasta):
        try:
            with open(output_fasta, 'w') as output:
                subprocess.run(f"cat {flye_fasta} {short_fasta}", stdout=output_fasta)

            logging.info("Flye and metaSPAdes/Unicycler assemblies combined...")
        except AttributeError:
            os.symlink(input_fasta, output_fasta)
            logging.info("Treating Flye assembly as final assembly...")


if __name__ == "__main__":
    flye_fasta = snakemake.input.flye_fasta
    short_fasta = snakemake.input.short_fasta
    input_fasta = snakemake.input.input_fasta
    output_fasta = snakemake.output.output_fasta

    combine_assemblies(
        flye_fasta=flye_fasta,
        short_fasta=short_fasta,
        input_fasta=input_fasta,
        output_fasta=output_fasta
    )