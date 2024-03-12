from subprocess import CalledProcessError, run, STDOUT
import os
from pathlib import Path
import logging
import gzip
import subprocess
import tempfile
import glob


class SingleMContainer:
    def __init__(self, threads: int, output_dir: str, genomes: str, assembly: str, pipe_results: str, logf: str):
        self.commands = []
        self.threads = threads
        self.pipe_results = pipe_results
        self.genomes = []
        # for each file ending in .fna in genomes folder, add to self.genomes
        for file in os.listdir(genomes):
            if file.endswith(".fna"):
                self.genomes.append(f"{genomes}/{file}")
        self.assembly = assembly
        self.output_dir = output_dir
        self.intermediate_dir = os.path.join(self.output_dir, "intermediate_genomes")
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.intermediate_dir, exist_ok=True)
        self.logf = logf
        self.process_queue = []

    
    def run(self):
        with open(self.logf, "a") as logf:
            logf.write("generating SingleM commands\n")
            self.create_commands()
            for command in self.commands:
                logf.write(" ".join(command) + "\n")
            logf.write("running SingleM commands\n")
            self.run_commands(logf)
            self.appraise_otu_tables(logf)
    
    def appraise_otu_tables(self, logf):
        logf.write("Appraising SingleM otu tables\n")
        genome_otu_tables = glob.glob(os.path.join(self.intermediate_dir, "*genome_single*.csv"))
        assembly_otu_table = os.path.join(self.intermediate_dir, "metagenome.assembly_0_otu_table.csv")
        appraise_cmd = f"singlem appraise --genome-otu-tables {' '.join(genome_otu_tables)} \
            --assembly-otu-table {assembly_otu_table} \
            --metagenome-otu-tables {self.pipe_results} \
            --plot {os.path.join(self.output_dir, 'singlem_appraise.svg')} \
            --output-binned-otu-table {os.path.join(self.output_dir, 'binned.otu_table.csv')} \
            --output-unbinned-otu-table {os.path.join(self.output_dir, 'unbinned.otu_table.csv')} \
            --output-assembled-otu-table {os.path.join(self.output_dir, 'assembled.otu_table.csv')}".split()
        try:
            # create output file: data/singlem_out/singlem_appraisal.tsv
            output_file = os.path.join(self.output_dir, "singlem_appraisal.tsv")
            with open(output_file, "w") as outf:
                with open(self.logf, "a") as logf:
                    run(appraise_cmd, stdout=outf, stderr=logf)
            Path("data/singlem_out/metagenome.combined_otu_table.csv").touch()
        except CalledProcessError as e:
            with open(log, "a") as logf:
                logf.write(e)
                logf.write("\nSingleM appraise failed. Exiting.\n")
            Path("data/singlem_out/metagenome.combined_otu_table.csv").touch()
   
    def create_commands(self):
        self._create_genome_commands()
        self._create_assembly_commands()


    def run_commands(self, logf):
        process_index = 0
        for command in self.commands:
            f = tempfile.TemporaryFile()
            logf.write(f"running process {process_index} {' '.join(command)}\n")
            p = subprocess.Popen(command, stdout=f, stderr=logf)
            self.process_queue.append((p, f))
            process_index += 1
            if len(self.process_queue) >= self.threads:
                self._check_processes(self.threads + 1, logf)
        
        # write how many processes are left
        logf.write(f"waiting for {len(self.process_queue)} processes to finish\n")
        while len(self.process_queue) > 0:
            self._check_processes(0, logf)

    def _check_processes(self, max_processes: int, logf):
        while len(self.process_queue) > max_processes:
            for i, (p, f) in enumerate(self.process_queue):
                if p.poll() is not None:
                    for line in f:
                        logf.write(line)
                    f.close()
                    self.process_queue.pop(i)
                    break


    def _create_assembly_commands(self):
        threads = max(self.threads // (len(self.genomes) + 1), 1)
        command = f"singlem pipe --threads {threads} --genome-fasta-files {self.assembly} --otu-table {self.intermediate_dir}/metagenome.assembly_0_otu_table.csv".split()
        self.commands.append(command)
    
    def _create_genome_commands(self):
        threads = max(self.threads // (len(self.genomes) + 1), 1)
        for i, genome in enumerate(self.genomes):
            command = f"singlem pipe --threads {threads} --genome-fasta-files {genome} --otu-table {self.intermediate_dir}/metagenome.genome_single_{i}_otu_table.csv".split()
            self.commands.append(command)


def run_singlem(
    genomes_folder: str,
    assembly: str,
    pipe_results: str,
    threads: int,
    log: str,
):
    output_dir = "data/singlem_out"
    singlem_container = SingleMContainer(threads, output_dir, genomes_folder, assembly, pipe_results, log)
    singlem_container.run()

def valid_path(path: str) -> bool:
    return os.path.exists(path)

if __name__ == '__main__':
    # check if SINGLEM_METAPACKAGE_PATH environment variable is set and path is valid
    # if not then, error and exit
    os.environ["SINGLEM_METAPACKAGE_PATH"] = snakemake.params.package_path
    if "SINGLEM_METAPACKAGE_PATH" not in os.environ or not valid_path(os.environ["SINGLEM_METAPACKAGE_PATH"]):
        raise ValueError("SINGLEM_METAPACKAGE_PATH environment variable not set. Please set using 'aviary configure' or manually. Exiting.")

    assembly = snakemake.input.assembly
    genomes = snakemake.params.genomes_folder
    pipe_results = snakemake.input.pipe_results
    threads = snakemake.threads
    log = snakemake.log[0]

    with open(log, "w") as logf: pass

    run_singlem(
        genomes,
        assembly,
        pipe_results,
        threads,
        log,
    )
