from subprocess import run
import os

long_reads = None
short_reads = None

def run_lorikeet(
    long_reads,
    short_reads_1,
    short_reads_2,
    mag_directory: str,
    output_directory: str,
    mag_extension: str,
    parallel_genomes: int,
    threads: int,
):
    if os.path.isdir("data/reads_mapped_to_mags/long/"):
        long_reads = "-l data/reads_mapped_to_mags/long/*.bam"
    elif long_reads != "none":
        long_reads = f"--longreads {' '.join(long_reads)}"
    else:
        long_reads = ""

    if os.path.isdir("data/reads_mapped_to_mags/short/"):
        short_reads = "-b data/reads_mapped_to_mags/short/*.bam"
    elif short_reads_2 != "none":
        short_reads = f"-1 {' '.join(short_reads_1)} -2 {' '.join(short_reads_2)}"
    elif short_reads_1 != "none":
        short_reads = f"--interleaved {' '.join(short_reads_1)}"
    else:
        short_reads = ""

    lorikeet_cmd = f"lorikeet call -t {threads} -d {mag_directory} -x {mag_extension} -P {parallel_genomes} {short_reads} {long_reads} -o {output_directory} --do-not-call-svs --calculate-dnds --calculate-fst".split()

    run(lorikeet_cmd)


if __name__ == '__main__':
    long_reads = snakemake.config['long_reads']
    short_reads_1 = snakemake.config['short_reads_1']
    short_reads_2 = snakemake.config['short_reads_2']
    mag_directory = snakemake.config['mag_directory']
    output_directory = snakemake.output.output_directory
    mag_extension = snakemake.params.mag_extension
    parallel_genomes = snakemake.params.parallel_genomes
    threads = snakemake.threads
    run_lorikeet(
        long_reads,
        short_reads_1,
        short_reads_2,
        mag_directory,
        output_directory,
        mag_extension,
        parallel_genomes,
        threads,
    )
