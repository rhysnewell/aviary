import subprocess
import os

long_reads = None
short_reads = None

if os.path.isdir("data/reads_mapped_to_mags/long/"):
    long_reads = "-l data/reads_mapped_to_mags/long/*.bam"
elif snakemake.config["long_reads"] != "none":
    long_reads = f"--longreads {' '.join(snakemake.config['long_reads'])}"
else:
    long_reads = ""

if os.path.isdir("data/reads_mapped_to_mags/short/"):
    short_reads = "-b data/reads_mapped_to_mags/short/*.bam"
elif snakemake.config["short_reads_2"] != "none":
    short_reads = f"-1 {' '.join(snakemake.config['short_reads_1'])} -2 {' '.join(snakemake.config['short_reads_2'])}"
elif snakemake.config["short_reads_1"] != "none":
    short_reads = f"--interleaved {' '.join(snakemake.config['short_reads_1'])}"
else:
    short_reads = ""

subprocess.Popen(f"lorikeet call -t {snakemake.threads} -d {snakemake.config['mag_directory']} "
                 f"-x {snakemake.params.mag_extension} -p {snakemake.params.parallel_genomes} {short_reads} {long_reads} "
                 f"-o {snakemake.output.output_directory} --calculate-dnds --calculate-fst", shell=True).wait()
