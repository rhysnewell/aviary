import subprocess
import os
import sys

with open(snakemake.config["batch_file"]) as batch_file:
    for line in batch_file:

        line = line.strip().split()
        assembly = line[0]
        identifier = line[1]
        # Assuming that batch file has reads written as interleaved files e.g. 1 2 1 2...
        # File paths must also be absolute
        forward = line[2::2]
        reverse = line[3::2]
        os.mkdir("data/" + identifier)

        # Symbolically link to main BinSnek folder
        subprocess.Popen("ln -s Snakefile data/%s/" % identifier, shell=True).wait()
        subprocess.Popen("cp config.yaml data/%s/" % identifier, shell=True).wait()
        subprocess.Popen("ln -s envs data/%s/" % identifier, shell=True).wait()
        subprocess.Popen("ln -s scripts data/%s/" % identifier, shell=True).wait()
        subprocess.Popen("ln -s .snakemake/ data/%s/" % identifier, shell=True).wait()

        # Flags to specify when to change the next line
        changing_fasta = False
        changing_reads_1 = False
        changing_reads_2 = False
        changing_batch = False
        with open("data/%s/config.yaml" % identifier, "w") as config:
            for config_line in config:
                if config_line.startswith("fasta"):
                    changing_fasta = True
                elif changing_fasta:
                    config_line = "\t" + assembly
                    changing_fasta = False
                elif config_line.startswith("short_reads_1"):
                    changing_reads_1 = True
                elif changing_reads_1:
                    config_line = "\t" + " ".join(forward)
                    changing_reads_1 = False
                elif config_line.startswith("short_reads_2"):
                    changing_reads_2 = True
                elif changing_reads_1:
                    config_line = "\t" + " ".join(reverse)
                    changing_reads_2 = False
                elif config_line.startswith("batch_file"):
                    changing_batch = True
                elif changing_batch:
                    config_line = "\t" + "none"
                    changing_batch = False

        os.chdir("data/%s" % identifier)
        # Run a new snakemake process using the updated config.yaml
        subprocess.Popen("snakemake --use-conda --cores %d recover_mags" % snakemake.threads, shell=True).wait()
        os.chdir("../../")

