#!/usr/bin/env python

#=======================================================================
# Author:
#
# Unit tests.
#
# Copyright
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License.
# If not, see <http://www.gnu.org/licenses/>.
#=======================================================================

import pytest
import os
import os.path
import subprocess
import shutil
import signal
import time
import unittest
import glob
import re
import gzip

data = os.path.join(os.path.dirname(__file__), 'data')

if os.environ.get("TEST_REQUEST_GPU", "0") == "1":
    request_gpu = "--request-gpu"
else:
    request_gpu = ""

singlem_metapackage = os.environ.get("SINGLEM_METAPACKAGE_PATH", "")
singlem_args = f"--skip-singlem false --singlem-metapackage-path {singlem_metapackage}" if singlem_metapackage else ""

def setup_output_dir(output_dir):
    try:
        shutil.rmtree(output_dir)
    except FileNotFoundError:
        pass
    os.makedirs(output_dir)

# We have a separate class for qsub tests (tests that run aviary with the CMR
# aqua snakemake profile) so that there is no need to run other expensive tests
# when running qsub tests.
@pytest.mark.qsub
class TestsQsub(unittest.TestCase):
    def test_short_read_recovery_queue_submission(self):
        output_dir = os.path.join("example", "test_short_read_recovery_queue_submission")
        setup_output_dir(output_dir)

        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"-n 32 -t 32 --local-cores 1 "
            f"--strict "
            f"--snakemake-profile aqua --cluster-retries 3 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertEqual(num_lines, 3)

    def test_short_read_recovery_queue_submission_gpus(self):
        output_dir = os.path.join("example", "test_short_read_recovery_queue_submission_gpus")
        setup_output_dir(output_dir)

        # Create inflated assembly file
        cmd = f"cat {data}/assembly.fasta > {output_dir}/assembly.fasta"
        multiplier = 100
        for i in range(multiplier):
            cmd += f" && awk '/^>/ {{print $0 \"{i}\"}} !/^>/ {{print $0}}' {data}/assembly.fasta >> {output_dir}/assembly.fasta"
        subprocess.run(cmd, shell=True, check=True)

        cmd = (
            f"aviary recover "
            f"--assembly {output_dir}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--skip-binners rosella metabat vamb "
            f"--extra-binners taxvamb comebin "
            f"--request-gpu "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            # The inflated/duplicated toy assembly this test bins against never
            # produces bins above ~48% completeness (no refinement is run, per
            # --refinery-max-iterations above), so the default 50%/5% filter
            # thresholds reject every bin and gtdbtk then fails on an empty
            # filtered_bins dir. Loosened here only for this smoke test.
            f"--min-completeness 20 "
            f"--max-contamination 20 "
            f"-n 32 -t 32 --local-cores 1 "
            f"--strict "
            f"--snakemake-profile aqua --cluster-retries 0 "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines >= 3)

@pytest.mark.expensive
class Tests(unittest.TestCase):
    def test_short_read_assembly(self):
        output_dir = os.path.join("example", "test_short_read_assembly")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

        with open(f"{output_dir}/aviary_out/data/final_contigs.fasta") as f:
            contigs = [c for c in f.read().strip().split('\n') if not c.startswith('>')]
            total_bp = sum(len(c) for c in contigs)

        self.assertTrue(total_bp > 1500000, "Assembly should be at least 1.5 million bp without host filtering")

    def test_short_read_assembly_host(self):
        output_dir = os.path.join("example", "test_short_read_assembly_host")
        setup_output_dir(output_dir)

        cmd = f"zcat {data}/GCA_000503915.1_ASM50391v1_genomic.fna.gz > {output_dir}/host_filter.fasta"
        subprocess.run(cmd, shell=True, check=True)

        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--host-filter {output_dir}/host_filter.fasta "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

        with open(f"{output_dir}/aviary_out/data/final_contigs.fasta") as f:
            contigs = [c for c in f.read().strip().split('\n') if not c.startswith('>')]
            total_bp = sum(len(c) for c in contigs)

        self.assertTrue(total_bp < 1000000, "Assembly should be smaller than 1 million bp after host filtering")

    def test_short_read_coassembly(self):
        output_dir = os.path.join("example", "test_short_read_coassembly")
        setup_output_dir(output_dir)

        cmd = f"cp {data}/wgsim.1.fq.gz {output_dir}/wgsimagain.1.fq.gz"
        cmd2 = f"cp {data}/wgsim.2.fq.gz {output_dir}/wgsimagain.2.fq.gz"
        subprocess.run(cmd + " && " + cmd2, shell=True, check=True)

        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz {output_dir}/wgsimagain.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz {output_dir}/wgsimagain.2.fq.gz "
            f"--coassemble yes "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

    def test_short_read_coassembly_skip_qc(self):
        output_dir = os.path.join("example", "test_short_read_coassembly_skip_qc")
        setup_output_dir(output_dir)

        cmd = f"cp {data}/wgsim.1.fq.gz {output_dir}/wgsimagain.1.fq.gz"
        cmd2 = f"cp {data}/wgsim.2.fq.gz {output_dir}/wgsimagain.2.fq.gz"
        subprocess.run(cmd + " && " + cmd2, shell=True, check=True)

        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz {output_dir}/wgsimagain.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz {output_dir}/wgsimagain.2.fq.gz "
            f"--coassemble yes --skip-qc "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

    def test_long_read_assembly_default(self):
        output_dir = os.path.join("example", "test_long_read_assembly")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--min-read-size 10 --min-mean-q 1 "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

    def test_long_read_assembly_flye(self):
        output_dir = os.path.join("example", "test_long_read_assembly_flye")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--long-read-assembler flye "
            f"--min-read-size 10 --min-mean-q 1 "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

    def test_long_read_assembly_no_short_reads(self):
        output_dir = os.path.join("example", "test_long_read_assembly_no_short_reads")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--min-read-size 10 --min-mean-q 1 "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

    def test_short_read_recovery(self):
        output_dir = os.path.join("example", "test_short_read_recovery")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"{singlem_args} "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 1)

        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

        if singlem_metapackage:
            self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/diversity"))
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraise.svg"))
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/benchmarks/singlem_appraise.benchmark.txt"))

    def test_long_read_recovery_split(self):
        output_dir = os.path.join("example", "test_long_read_recovery_split")
        setup_output_dir(output_dir)

        for i, size in enumerate([80000, 50000, 20000]):
            for end in [1, 2]:
                cmd = f"zcat {data}/wgsim.{end}.fq.gz | head -n {size} | gzip > {output_dir}/wgsim_{i}.{end}.fq.gz"
                subprocess.run(cmd, shell=True, check=True)

        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz {output_dir}/wgsim_0.1.fq.gz {output_dir}/wgsim_1.1.fq.gz {output_dir}/wgsim_2.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz {output_dir}/wgsim_0.2.fq.gz {output_dir}/wgsim_1.2.fq.gz {output_dir}/wgsim_2.2.fq.gz "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--coassemble no "
            f"--coverage-job-strategy always "
            f"--coverage-samples-per-job 2 "
            f"--min-read-size 10 --min-mean-q 1 "
            f"-n 32 -t 32 "
            f"{singlem_args} "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

        if singlem_metapackage:
            self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/diversity"))
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraise.svg"))

    def test_long_read_recovery_default(self):
        output_dir = os.path.join("example", "test_long_read_recovery")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--min-read-size 10 --min-mean-q 1 "
            f"-n 32 -t 32 "
            f"--strict "
            f"{singlem_args} "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

        if singlem_metapackage:
            self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/diversity"))
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraise.svg"))

    def test_long_read_only_recovery(self):
        output_dir = os.path.join("example", "test_long_read_only_recovery")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--min-read-size 10 --min-mean-q 1 "
            f"{singlem_args} "
            f"-n 32 -t 32 "
            f"--strict "
        )
        print(cmd)
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

        if singlem_metapackage:
            self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/diversity"))
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraise.svg"))

    def test_short_read_recovery_fast(self):
        output_dir = os.path.join("example", "test_short_read_recovery_fast")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella vamb metabat "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertEqual(num_lines, 3)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_semibin(self):
        suffix = "_gpu" if os.environ.get("TEST_REQUEST_GPU") == "1" else ""
        output_dir = os.path.join("example", f"test_short_read_recovery_semibin{suffix}")
        setup_output_dir(output_dir)

        cmd = f"ln -sr {data}/wgsim.1.fq.gz {output_dir}/wgsim2.1.fq.gz && ln -sr {data}/wgsim.2.fq.gz {output_dir}/wgsim2.2.fq.gz"
        subprocess.run(cmd, shell=True, check=True)

        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz {output_dir}/wgsim2.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz {output_dir}/wgsim2.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella vamb metabat "
            f"{request_gpu} "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertEqual(num_lines, 3)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
    
    def test_short_read_recovery_semibin_multi(self):
        suffix = "_gpu" if os.environ.get("TEST_REQUEST_GPU") == "1" else ""
        output_dir = os.path.join("example", f"test_short_read_recovery_semibin_multi{suffix}")
        setup_output_dir(output_dir)

        # Create a second assembly with renamed contig headers to simulate a separate sample.
        # Appending "_s2" to NODE names ensures contig names are unique across the two assemblies.
        assembly2 = os.path.join(output_dir, "assembly2.fasta")
        cmd = (
            f"awk '/^>/ {{print $0 \"_s2\"}} !/^>/ {{print $0}}' "
            f"{data}/assembly.fasta > {assembly2}"
        )
        subprocess.run(cmd, shell=True, check=True)

        # Symlink a second read set to simulate reads from the second sample.
        cmd = (
            f"ln -sr {data}/wgsim.1.fq.gz {output_dir}/wgsim2.1.fq.gz && "
            f"ln -sr {data}/wgsim.2.fq.gz {output_dir}/wgsim2.2.fq.gz"
        )
        subprocess.run(cmd, shell=True, check=True)

        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta {assembly2} "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz {output_dir}/wgsim2.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz {output_dir}/wgsim2.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella vamb metabat "
            f"--semibin-mode multi "
            f"{request_gpu} "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertGreater(num_lines, 1)

        # Two assemblies → two sample directories under samples/
        # (SemiBin2 also writes .fa and .csv files into the same directory, so filter dirs only)
        semibin_sample_dirs = [
            d for d in glob.glob(f"{output_dir}/aviary_out/data/semibin_bins/samples/*")
            if os.path.isdir(d)
        ]
        self.assertEqual(len(semibin_sample_dirs), 2)

        semibin_bins = glob.glob(f"{output_dir}/aviary_out/data/semibin_bins/output_bins/*.fa")
        self.assertGreater(len(semibin_bins), 0)

        # Bin filenames must be prefixed with their sample name
        sample_names = [os.path.basename(d) for d in semibin_sample_dirs]
        bin_names = [os.path.basename(b) for b in semibin_bins]
        self.assertTrue(
            any(
                bin_name.startswith(f"{sample_name}_")
                for sample_name in sample_names
                for bin_name in bin_names
            )
        )

        # Confirm the concatenated FASTA and per-sample BAMs were created
        self.assertTrue(os.path.isfile(
            f"{output_dir}/aviary_out/data/semibin_multi_prep/concatenated.fa"
        ))
        self.assertGreater(
            len(glob.glob(f"{output_dir}/aviary_out/data/semibin_multi_bams/*.bam")), 0
        )

        semibin_logs = glob.glob(
            f"{output_dir}/aviary_out/logs/semibin/*/attempt*.log"
        )
        self.assertGreater(len(semibin_logs), 0)
        self.assertTrue(
            any(
                "Performing multi-sample binning" in open(log_path).read()
                for log_path in semibin_logs
            ),
            "SemiBin log did not show multi-sample binning was run",
        )

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_vamb(self):
        output_dir = os.path.join("example", "test_short_read_recovery_vamb")
        setup_output_dir(output_dir)

        # Create inflated assembly file
        cmd = f"cat {data}/assembly.fasta > {output_dir}/assembly.fasta"
        multiplier = 100
        for i in range(multiplier):
            cmd += f" && awk '/^>/ {{print $0 \"{i}\"}} !/^>/ {{print $0}}' {data}/assembly.fasta >> {output_dir}/assembly.fasta"

        subprocess.run(cmd, shell=True, check=True)

        cmd = (
            f"aviary recover "
            f"--assembly {output_dir}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella semibin metabat "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 2)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_rosella(self):
        output_dir = os.path.join("example", "test_short_read_recovery_rosella")
        setup_output_dir(output_dir)

        # Create inflated assembly file
        cmd = f"cat {data}/assembly.fasta > {output_dir}/assembly.fasta"
        multiplier = 100
        for i in range(multiplier):
            cmd += f" && awk '/^>/ {{print $0 \"{i}\"}} !/^>/ {{print $0}}' {data}/assembly.fasta >> {output_dir}/assembly.fasta"

        subprocess.run(cmd, shell=True, check=True)

        cmd = (
            f"aviary recover "
            f"--assembly {output_dir}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners semibin metabat vamb "
            f"{request_gpu} "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 2)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_taxvamb(self):
        suffix = "_gpu" if os.environ.get("TEST_REQUEST_GPU") == "1" else ""
        output_dir = os.path.join("example", f"test_short_read_recovery_taxvamb{suffix}")
        logs_dir = os.path.join(output_dir, "aviary_out", "logs")
        logs_backup = logs_dir + ".bak"
        if os.path.exists(logs_dir):
            shutil.copytree(logs_dir, logs_backup)
        setup_output_dir(output_dir)
        if os.path.exists(logs_backup):
            os.makedirs(os.path.join(output_dir, "aviary_out"), exist_ok=True)
            shutil.copytree(logs_backup, logs_dir)
            shutil.rmtree(logs_backup)

        # Create inflated assembly file
        cmd = f"cat {data}/assembly.fasta > {output_dir}/assembly.fasta"
        multiplier = 100
        for i in range(multiplier):
            cmd += f" && awk '/^>/ {{print $0 \"{i}\"}} !/^>/ {{print $0}}' {data}/assembly.fasta >> {output_dir}/assembly.fasta"

        subprocess.run(cmd, shell=True, check=True)

        cmd = (
            f"aviary recover "
            f"--assembly {output_dir}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella semibin metabat vamb "
            f"--extra-binners taxvamb "
            f"{request_gpu} "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 2)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_comebin(self):
        suffix = "_gpu" if os.environ.get("TEST_REQUEST_GPU") == "1" else ""
        output_dir = os.path.join("example", f"test_short_read_recovery_comebin{suffix}")
        logs_dir = os.path.join(output_dir, "aviary_out", "logs")
        logs_backup = logs_dir + ".bak"
        if os.path.exists(logs_dir):
            shutil.copytree(logs_dir, logs_backup)
        setup_output_dir(output_dir)
        if os.path.exists(logs_backup):
            os.makedirs(os.path.join(output_dir, "aviary_out"), exist_ok=True)
            shutil.copytree(logs_backup, logs_dir)
            shutil.rmtree(logs_backup)

        # Create inflated assembly file
        cmd = f"cat {data}/assembly.fasta > {output_dir}/assembly.fasta"
        multiplier = 100
        for i in range(multiplier):
            cmd += f" && awk '/^>/ {{print $0 \"{i}\"}} !/^>/ {{print $0}}' {data}/assembly.fasta >> {output_dir}/assembly.fasta"

        subprocess.run(cmd, shell=True, check=True)

        cmd = (
            f"aviary recover "
            f"--assembly {output_dir}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella semibin metabat vamb "
            f"--extra-binners comebin "
            f"{request_gpu} "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 2)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_quickbin(self):
        output_dir = os.path.join("example", "test_short_read_recovery_quickbin")
        setup_output_dir(output_dir)

        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella semibin metabat vamb "
            f"--extra-binners quickbin "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 2)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_no_bins(self):
        output_dir = os.path.join("example", "test_short_read_recovery_no_bins")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/tiny_sample.1.fq "
            f"-2 {data}/tiny_sample.2.fq "
            f"--binning-only "
            f"--skip-binners rosella vamb semibin metabat1 "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertFalse(os.path.isfile(bin_info_path))

    def test_short_read_recovery_no_assembly_provided(self):
        output_dir = os.path.join("example", "test_short_read_recovery_no_assembly_provided")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella vamb metabat "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertEqual(num_lines, 3)

        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

    def test_short_read_complete(self):
        output_dir = os.path.join("example", "test_short_read_complete")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary complete "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"{request_gpu} "
            f"{singlem_args} "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 1)

        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

        if singlem_metapackage:
            self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/diversity"))
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraise.svg"))

        gtdbtk_path = f"{output_dir}/aviary_out/taxonomy/gtdbtk.bac120.summary.tsv"
        self.assertTrue(os.path.isfile(gtdbtk_path))
        with open(gtdbtk_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 1)

        eggnog_paths = glob.glob(os.path.join(output_dir, "aviary_out", "annotation", "eggnog", "*.annotations"))
        self.assertTrue(len(eggnog_paths) > 0)

        for eggnog_path in eggnog_paths:
            self.assertTrue(os.path.isfile(eggnog_path))
            with open(eggnog_path) as f:
                num_lines = sum(1 for _ in f)
            self.assertTrue(num_lines > 1)

    def test_error_integration(self):
        """Expect aviary_assemble to fail with tiny test data, then check logging.
        The small input cannot be assembled, so aviary_assemble will error first.
        We assert Snakemake reports the aviary_assemble error, our handler logs the
        rule failure with a concrete log path, and that an assemble log exists.
        """
        output_dir = os.path.join("example", "test_error_integration")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/tiny_sample_bad.1.fq "
            f"-2 {data}/tiny_sample.2.fq "
            f"--binning-only "
            f"--skip-qc --coassemble "
            f"-n 32 -t 32 "
            f"--strict "
        )

        try:
            proc = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            combined = (proc.stdout or b"").decode("utf-8", errors="replace") + (proc.stderr or b"").decode("utf-8", errors="replace")
        except subprocess.CalledProcessError as e:
            combined = (e.stdout or b"").decode("utf-8", errors="replace") + (e.stderr or b"").decode("utf-8", errors="replace")

        # Snakemake should report a rule error; our wrapper should emit log dumps
        self.assertTrue(("Error in rule assemble_short_reads" in combined) or ("RuleException in rule assemble_short_reads" in combined))
        self.assertIn("Rule failed: assemble_short_reads", combined)
        self.assertTrue(("===== BEGIN LOG (" in combined) and ("===== END LOG (" in combined))

        # Ensure an assemble log exists under the expected logs tree
        assemble_logs = glob.glob(os.path.join(output_dir, "aviary_out", "logs", "assemble_short_reads", "*", "attempt*.log"))
        self.assertTrue(len(assemble_logs) > 0)

        # One of the printed log paths should match an existing file
        printed_paths = re.findall(r"BEGIN LOG \([^)]*\):\s*(.*?)\s*=====", combined)
        if printed_paths:
            # Normalize whitespace and test for existence of at least one printed log
            self.assertTrue(any(os.path.exists(p.strip()) for p in printed_paths))


    def test_assembly_spades_multi_read_coassemble(self):
        """Verify SPAdes coassembly does not collapse distinct-but-identical-content readsets.

        assemble_short_reads.py's cat-based concatenation (aviary/modules/assembly/scripts/
        assemble_short_reads.py, ~line 49) only runs when --coassemble is set AND more than one
        distinct read-set path is supplied; a single-element list is treated like an
        uncoassembled run and concatenation is skipped entirely. Two readsets with identical
        content but distinct file paths are used here, so that a path-deduplication bug upstream
        (which would collapse the two-element list back to one and skip concatenation) is
        distinguishable from correct behaviour: the concatenated output must contain exactly
        twice the read count of a single input.
        """
        output_dir = os.path.join("example", "test_assembly_spades_multi_read_coassemble")
        setup_output_dir(output_dir)
        aviary_out = os.path.join(output_dir, "aviary_out")

        # Create copies of the test reads under distinct paths so aviary sees
        # two separate readsets and does not deduplicate them.
        reads2_1 = os.path.join(output_dir, "reads2.1.fq.gz")
        reads2_2 = os.path.join(output_dir, "reads2.2.fq.gz")
        shutil.copy(f"{data}/wgsim.1.fq.gz", reads2_1)
        shutil.copy(f"{data}/wgsim.2.fq.gz", reads2_2)

        final_assembly = os.path.join(aviary_out, "data", "final_contigs.fasta")
        concatenated_reads = os.path.join(aviary_out, "data", "short_reads.1.fastq.gz")
        concatenated_reads_2 = os.path.join(aviary_out, "data", "short_reads.2.fastq.gz")

        subprocess.run(
            f"aviary assemble "
            f"-o {aviary_out} "
            f"-1 {data}/wgsim.1.fq.gz {reads2_1} "
            f"-2 {data}/wgsim.2.fq.gz {reads2_2} "
            f"--skip-qc --coassemble "
            f"-n 32 -t 32 ",
            shell=True, check=True
        )

        self.assertTrue(
            os.path.isfile(final_assembly),
            "Final assembly missing after SPAdes multi-read coassemble with skip-qc"
        )
        self.assertTrue(
            os.path.isfile(concatenated_reads),
            "Concatenated forward reads file missing — assemble_short_reads.py did not concatenate"
        )
        self.assertTrue(
            os.path.isfile(concatenated_reads_2),
            "Concatenated reverse reads file missing — assemble_short_reads.py did not concatenate"
        )

        def gzip_line_count(path):
            with gzip.open(path, 'rt') as f:
                return sum(1 for _ in f)

        # The two input paths have identical content, so correct concatenation
        # doubles the line count. If the two distinct paths were deduplicated
        # back to one upstream, only a single copy would be present here.
        self.assertEqual(
            gzip_line_count(concatenated_reads),
            gzip_line_count(f"{data}/wgsim.1.fq.gz") * 2,
            "Forward concatenated reads are not double the single-readset line count — "
            "the two distinct paths were likely deduplicated back to one, skipping concatenation"
        )
        self.assertEqual(
            gzip_line_count(concatenated_reads_2),
            gzip_line_count(f"{data}/wgsim.2.fq.gz") * 2,
            "Reverse concatenated reads are not double the single-readset line count — "
            "the two distinct paths were likely deduplicated back to one, skipping concatenation"
        )


if __name__ == "__main__":
    unittest.main()
