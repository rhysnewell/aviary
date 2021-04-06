import subprocess
import os
import sys


def process_batch(batch_file_path):
    main_directory = os.getcwd()
    install_directory = "~/git/aviary"
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
            # Symbolically link to main aviary folder === BAD IDEA ===
            # subprocess.Popen("ln -s %s/annotation.smk %s/data/%s/"
            #                  % (main_directory, main_directory, identifier), shell=True).wait()
            # subprocess.Popen("cp %s/template_config.yaml %s/data/%s/"
                             # % (main_directory, main_directory, identifier), shell=True).wait()
            subprocess.Popen("ln -s %s/envs %s/data/%s/"
                             % (main_directory, main_directory, identifier), shell=True).wait()
            subprocess.Popen("ln -s %s/scripts %s/data/%s/"
                             % (main_directory, main_directory, identifier), shell=True).wait()
            subprocess.Popen("ln -s %s/.snakemake/ %s/data/%s/"
                             % (main_directory, main_directory, identifier), shell=True).wait()

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
            # Run a new snakemake process using the updated template_config.yaml
            subprocess.Popen("snakemake --unlock --use-conda --conda-prefix %s/.snakemake/ -s %s/annotation.smk recover_mags"
                             % (install_directory, install_directory), shell=True).wait()
            subprocess.Popen("snakemake --use-conda --conda-prefix %s/.snakemake/ -s %s/annotation.smk --cores %d recover_mags"
                             % (install_directory, install_directory, snakemake.threads), shell=True).wait()
            os.chdir(main_directory)


if __name__ == "__main__":
    process_batch(snakemake.config["batch_file"])
    subprocess.Popen("touch data/done", shell=True).wait()