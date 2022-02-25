import subprocess
import os
import logging


if not os.path.exists(snakemake.output.fasta):
    try:
        subprocess.Popen("cat %s %s > %s" % (snakemake.input.flye_fasta,
                                             snakemake.input.short_fasta,
                                             snakemake.output.fasta), shell=True).wait()
        logging.info("Flye and metaSPAdes/Unicycler assemblies combined...")
    except AttributeError:
        subprocess.Popen("ln %s %s" % (snakemake.input.fasta,
                                       snakemake.output.fasta), shell=True).wait()
        logging.info("Treating Flye assembly as final assembly...")