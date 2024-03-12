from subprocess import CalledProcessError, run, STDOUT
import os
from pathlib import Path
import logging
import gzip
import subprocess
import tempfile
import glob

class ReadContainer:
    def __init__(self, long_reads, short_reads_1, short_reads_2):
        if long_reads == "none":
            long_reads = []

        if short_reads_1 == "none":
            short_reads_1 = []
        
        if short_reads_2 == "none":
            short_reads_2 = []

        self.long_reads = long_reads
        self.short_reads_1 = short_reads_1
        self.short_reads_2 = short_reads_2
    
    def _check_paired_reads(self):
        if len(self.short_reads_2) == 0:
            return False
        if len(self.short_reads_1) == 0:
            return False
        if len(self.short_reads_1) != len(self.short_reads_2):
            logging.warning("Short read 1 and short read 2 are not the same length but paired reads assumed, please check input")
            return False
        return True

    def get_paired_reads(self):
        if not self._check_paired_reads():
            return []

        paired_reads = []
        for i in range(len(self.short_reads_1)):
            paired_reads.append((self.short_reads_1[i], self.short_reads_2[i]))
        return paired_reads
    
    def get_paired_read_count(self):
        if not self._check_paired_reads() or self._check_interleaved():
            return 0
        return len(self.short_reads_1)

    def get_single_reads(self):
        if self._check_paired_reads() or self._check_interleaved():
            return []
        return self.short_reads_1 if len(self.short_reads_1) > 0 else self.short_reads_2

    def get_single_read_count(self):
        if self._check_paired_reads() or self._check_interleaved():
            return 0
        return len(self.short_reads_1) if len(self.short_reads_1) > 0 else len(self.short_reads_2)

    def get_long_read_count(self):
        return len(self.long_reads)

    def get_long_reads(self):
        return self.long_reads
    
    def get_interleaved_reads(self):
        if not self._check_interleaved() or self._check_paired_reads():
            return []
        return self.short_reads_1 if len(self.short_reads_1) > 0 else self.short_reads_2
    
    def get_interleaved_read_count(self):
        if not self._check_interleaved() or self._check_paired_reads():
            return 0
        return len(self.short_reads_1) if len(self.short_reads_1) > 0 else len(self.short_reads_2)
    
    def _check_interleaved(self):
        if self._check_paired_reads():
            return False
        if len(self.short_reads_1) == 0 and len(self.short_reads_2) == 0:
            return False

        short_reads = self.short_reads_1 if len(self.short_reads_1) > 0 else self.short_reads_2
        # open the reads and check if they are interleaved
        # if they are then return True
        # if they are not then return False
        
        # check if first file is gzipped
        interleaved = False
        if short_reads[0].endswith(".gz"):
            with gzip.open(short_reads[0], "rt") as read_file:
                interleaved = self._forward_and_reverse_present(read_file)
        else:
            with open(short_reads[0], "r") as read_file:
                interleaved = self._forward_and_reverse_present(read_file)

        return interleaved
    
    def _forward_and_reverse_present(self, file) -> bool:
        read_pattern = []
        for i, line in enumerate(file):
            if i > 32:
                break
            if i % 4 == 0:
                if line.endswith("1"):
                    read_pattern.append(1)
                elif line.endswith("2"):
                    read_pattern.append(2)
                else:
                    read_pattern.append(0)
        
        
        for r1, r2 in zip(read_pattern[::2], read_pattern[1::2]):
            if r1 == 0 or r2 == 0:
                return False
            if r1 == r2:
                return False
            if r1 == 1 and r2 == 2:
                return True
        
        return False
    
    def get_total_read_count(self):
        return self.get_long_read_count() + self.get_single_read_count() + self.get_paired_read_count() + self.get_interleaved_read_count()

class SingleMContainer:
    def __init__(self, threads: int, output_dir: str, read_container: ReadContainer, logf: str):
        self.commands = []
        self.threads = threads
        self.read_container = read_container
        self.output_dir = output_dir
        self.intermediate_dir = os.path.join(self.output_dir, "intermediate")
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.intermediate_dir, exist_ok=True)
        self.logf = logf
        self.process_queue = []
        self.total_reads = self.read_container.get_total_read_count()

    
    def run(self):
        with open(self.logf, "a") as logf:
            logf.write("generating SingleM commands\n")
            self.create_commands()
            for command in self.commands:
                logf.write(" ".join(command) + "\n")
            logf.write("running SingleM commands\n")
            self.run_commands(logf)
            self.combine_otu_tables(logf)
    
    def combine_otu_tables(self, logf):
        logf.write("combining SingleM otu tables\n")
        intermediate_otu_tables = glob.glob(os.path.join(self.intermediate_dir, "*.csv"))
        summarise_cmd = f"singlem summarise --input-otu-tables {' '.join(intermediate_otu_tables)} --output-otu-table {os.path.join(self.output_dir, 'metagenome.combined_otu_table.csv')}".split()
        try:
            with open(self.logf, "a") as logf:
                run(summarise_cmd, stdout=logf, stderr=STDOUT)
            Path("data/singlem_out/metagenome.combined_otu_table.csv").touch()
        except CalledProcessError as e:
            with open(log, "a") as logf:
                logf.write(e)
                logf.write("\nSingleM summarise failed. Exiting.\n")
            Path("data/singlem_out/metagenome.combined_otu_table.csv").touch()
   
    def create_commands(self):
        self._create_longread_commands()
        self._create_shortread_commands()


    def run_commands(self, logf):
        process_index = 0
        for command in self.commands:
            f = tempfile.TemporaryFile()
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


    def _create_longread_commands(self):
        threads = max(self.threads // self.total_reads, 1)
        for i, long_read in enumerate(self.read_container.get_long_reads()):
            command = f"singlem pipe --threads {threads} --sequences {long_read} --otu-table {self.intermediate_dir}/metagenome.longread_{i}_otu_table.csv".split()
            self.commands.append(command)
    
    def _create_shortread_commands(self):
        threads = max(self.threads // self.total_reads, 1)
        for i, short_read in enumerate(self.read_container.get_single_reads()):
            command = f"singlem pipe --threads {threads} --sequences {short_read} --otu-table {self.intermediate_dir}/metagenome.shortread_single_{i}_otu_table.csv".split()
            self.commands.append(command)
        
        for i, (short_read_1, short_read_2) in enumerate(self.read_container.get_paired_reads()):
            command = f"singlem pipe --threads {threads} --forward {short_read_1} --reverse {short_read_2} --otu-table {self.intermediate_dir}/metagenome.shortread_paired_{i}_otu_table.csv".split()
            self.commands.append(command)
        
        for i, interleaved_read in enumerate(self.read_container.get_interleaved_reads()):
            command = f"singlem pipe --threads {threads} --sequences {interleaved_read} --otu-table {self.intermediate_dir}/metagenome.shortread_interleaved_{i}_otu_table.csv".split()
            self.commands.append(command)



def run_singlem(
    read_container: ReadContainer,
    threads: int,
    log: str,
):
    output_dir = "data/singlem_out"
    singlem_container = SingleMContainer(threads, output_dir, read_container, log)
    singlem_container.run()

def valid_path(path: str) -> bool:
    return os.path.exists(path)

if __name__ == '__main__':
    # check if SINGLEM_METAPACKAGE_PATH environment variable is set and path is valid
    # if not then, error and exit
    os.environ["SINGLEM_METAPACKAGE_PATH"] = snakemake.params.package_path
    if "SINGLEM_METAPACKAGE_PATH" not in os.environ or not valid_path(os.environ["SINGLEM_METAPACKAGE_PATH"]):
        raise ValueError("SINGLEM_METAPACKAGE_PATH environment variable not set. Please set using 'aviary configure' or manually. Exiting.")

    long_reads = snakemake.config['long_reads']
    short_reads_1 = snakemake.config['short_reads_1']
    short_reads_2 = snakemake.config['short_reads_2']
    pplacer_threads = snakemake.threads
    log = snakemake.log[0]

    read_container = ReadContainer(long_reads, short_reads_1, short_reads_2)

    with open(log, "w") as logf: pass

    run_singlem(
        read_container,
        pplacer_threads,
        log,
    )
