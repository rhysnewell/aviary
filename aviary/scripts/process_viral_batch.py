from subprocess import run
import os
from pathlib import Path


def process_batch(
    batch_file_path: str,
    threads: int,
    ):
    main_directory = os.getcwd()
    with open(batch_file_path, "r") as batch_file:
        for line in batch_file:
            line = line.strip().split()
            assembly = line[0]
            identifier = line[1]
            # Assuming that batch file has reads written as interleaved files e.g. 1 2 1 2...
            # File paths must also be absolute
            forward = line[2::2]
            reverse = line[3::2]
            try:
                os.mkdir("data/" + identifier)
            except FileExistsError:
                print("Directory already exists for sample %s" % identifier)
            # Symbolically link to main aviary folder
            os.symlink(f"{main_directory}/envs", f"{main_directory}/data/{identifier}/envs")
            os.symlink(f"{main_directory}/scripts", f"{main_directory}/data/{identifier}/scripts")
            os.symlink(f"{main_directory}/.snakemake", f"{main_directory}/data/{identifier}/.snakemake")

            # Flags to specify when to change the next line
            changing_fasta = False
            changing_reads_1 = False
            changing_reads_2 = False
            changing_batch = False
            with open("%s/data/%s/template_config.yaml" % (main_directory, identifier), "w+") as new_config_file:
                with open("%s/template_config.yaml" % main_directory, 'r') as config_file:
                    for config_line in config_file:
                        if config_line.startswith("fasta"):
                            new_config_file.write(config_line)
                            changing_fasta = True
                        elif changing_fasta:
                            new_config_file.write(" " + assembly + "\n")
                            changing_fasta = False
                        elif config_line.startswith("short_reads_1"):
                            new_config_file.write(config_line)
                            changing_reads_1 = True
                        elif changing_reads_1:
                            config_line = " " + " ".join(forward) + "\n"
                            new_config_file.write(config_line)
                            changing_reads_1 = False
                        elif config_line.startswith("short_reads_2"):
                            new_config_file.write(config_line)
                            changing_reads_2 = True
                        elif changing_reads_2:
                            config_line = " " + " ".join(reverse) + "\n"
                            new_config_file.write(config_line)
                            changing_reads_2 = False
                        elif config_line.startswith("batch_file"):
                            new_config_file.write(config_line)
                            changing_batch = True
                        elif changing_batch:
                            config_line = " " + "none" + "\n"
                            new_config_file.write(config_line)
                            changing_batch = False
                        else:
                            new_config_file.write(config_line)


            os.chdir("%s/data/%s" % (main_directory, identifier))
            # Run a new snakemake process using the updated template_config.yam
            unlock_cmd = f"snakemake --unlock --use-conda --conda-prefix {main_directory}/.snakemake/ -s {main_directory}/assembly.smk recover_mags".split()
            run(unlock_cmd)

            run_cmd = f"snakemake --use-conda --conda-prefix {main_directory}/.snakemake/ -s {main_directory}/assembly.smk --cores {threads} recover_mags".split()
            run(run_cmd)
            os.chdir(main_directory)


if __name__ == "__main__":
    process_batch(snakemake.config["batch_file"], snakemake.threads)
    Path("data/done").touch()