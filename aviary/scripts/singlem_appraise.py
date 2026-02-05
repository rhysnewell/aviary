#!/usr/bin/env python3

import argparse
from subprocess import run
import os
import sys
import extern

class SingleMContainer:
    def __init__(self, threads: int, output_dir: str, genomes: str, assembly: str, pipe_results: str):
        self.commands = []
        self.threads = threads
        self.pipe_results = pipe_results
        self.genomes_dir = genomes
        self.assembly = assembly
        self.output_dir = output_dir
        self.intermediate_dir = os.path.join(self.output_dir, "intermediate_genomes")
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.intermediate_dir, exist_ok=True)

    def run(self):
        print("generating SingleM commands", file=sys.stderr)
        self.create_commands()
        for command in self.commands:
            print(command, file=sys.stderr)
        print("running SingleM commands", file=sys.stderr)
        self.run_commands()
        self.appraise_otu_tables()

    def appraise_otu_tables(self):
        print("Appraising SingleM otu tables", file=sys.stderr)
        genome_otu_tables = os.path.join(self.intermediate_dir, "*genomes_otu_table.csv")
        assembly_otu_table = os.path.join(self.intermediate_dir, "metagenome.assembly_0_otu_table.csv")
        output_file = os.path.join(self.output_dir, "singlem_appraisal.tsv")
        appraise_cmd = f"singlem appraise --genome-otu-tables {genome_otu_tables} \
            --assembly-otu-table {assembly_otu_table} \
            --metagenome-otu-tables {self.pipe_results} \
            --plot {os.path.join(self.output_dir, 'singlem_appraise.svg')} \
            --output-binned-otu-table {os.path.join(self.output_dir, 'binned.otu_table.csv')} \
            --output-unbinned-otu-table {os.path.join(self.output_dir, 'unbinned.otu_table.csv')} \
            --output-unaccounted-for-otu-table {os.path.join(self.output_dir, 'unaccounted_for.otu_table.csv')} \
            > {output_file}"
        extern.run(appraise_cmd)
        print(f"SingleM appraisal complete. Results in {output_file}", file=sys.stderr)
   
    def create_commands(self):
        self._create_genome_commands()
        self._create_assembly_commands()


    def run_commands(self):
        process_index = 0
        for command in self.commands:
            print(f"running process {process_index} {command}", file=sys.stderr)
            extern.run(command)

    def _create_assembly_commands(self):
        threads = self.threads
        command = f"singlem pipe --threads {threads} --genome-fasta-files {self.assembly} --otu-table {self.intermediate_dir}/metagenome.assembly_0_otu_table.csv"
        self.commands.append(command)
    
    def _create_genome_commands(self):
        threads = self.threads
        command = f"singlem pipe --threads {threads} --genome-fasta-directory {self.genomes_dir} --otu-table {self.intermediate_dir}/metagenome.genomes_otu_table.csv"
        self.commands.append(command)


def run_singlem(
    genomes_folder: str,
    assembly: str,
    pipe_results: str,
    threads: int,
    package_path: str,
):
    # Set SINGLEM_METAPACKAGE_PATH environment variable
    os.environ["SINGLEM_METAPACKAGE_PATH"] = package_path
    if not valid_path(os.environ["SINGLEM_METAPACKAGE_PATH"]):
        raise ValueError("SINGLEM_METAPACKAGE_PATH environment variable not valid. Please set using 'aviary configure' or manually. Exiting.")
        
    output_dir = "data/singlem_out"
    singlem_container = SingleMContainer(threads, output_dir, genomes_folder, assembly, pipe_results)
    singlem_container.run()

def valid_path(path: str) -> bool:
    return os.path.exists(path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run SingleM appraise on genomes and assembly')
    parser.add_argument('--genomes-folder', required=True, help='Folder containing genome fasta files')
    parser.add_argument('--assembly', required=True, help='Assembly file path')
    parser.add_argument('--pipe-results', required=True, help='Path to SingleM pipe results')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    parser.add_argument('--package-path', required=True, help='Path to SingleM metapackage')
    
    args = parser.parse_args()

    run_singlem(
        args.genomes_folder,
        args.assembly,
        args.pipe_results,
        args.threads,
        args.package_path,
    )
