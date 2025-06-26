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
import unittest

data = os.path.join(os.path.dirname(__file__), 'data')

if os.environ.get("TEST_REQUEST_GPU", "0") == "1":
    request_gpu = "--request-gpu"
else:
    request_gpu = ""

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

    def test_long_read_assembly(self):
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

    def test_short_read_recovery(self):
        output_dir = os.path.join("example", "test_short_read_recovery")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
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

        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/diversity"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv"))
        self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv") > 0)
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv"))
        self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv") > 0)
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraise.svg"))

    def test_long_read_recovery_split(self):
        output_dir = os.path.join("example", "test_long_read_recovery_split")
        setup_output_dir(output_dir)

        for i, size in enumerate([80000, 50000, 20000]):
            for end in [1, 2]:
                cmd = f"zcat {data}/wgsim.{end}.fq.gz | head -n {size} > {output_dir}/wgsim_{i}.{end}.fq.gz"
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
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/diversity"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv"))
        self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv") > 0)
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv"))
        self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv") > 0)
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraise.svg"))

    def test_long_read_recovery(self):
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
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

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
            f"-n 32 -t 32 "
            f"--strict "
        )
        print(cmd)
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

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

        semibin_log_path = f"{output_dir}/aviary_out/logs/semibin.log"
        self.assertTrue(os.path.isfile(semibin_log_path))
        with open(semibin_log_path) as f:
            log = f.read()
        self.assertTrue("Training model..." not in log)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_semibin(self):
        output_dir = os.path.join("example", "test_short_read_recovery_semibin")
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

        semibin_log_path = f"{output_dir}/aviary_out/logs/semibin.log"
        self.assertTrue(os.path.isfile(semibin_log_path))
        with open(semibin_log_path) as f:
            log = f.read()
        self.assertTrue("Training model..." in log)

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

    def test_short_read_recovery_taxvamb(self):
        output_dir = os.path.join("example", "test_short_read_recovery_taxvamb")
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
        output_dir = os.path.join("example", "test_short_read_recovery_comebin")
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

if __name__ == "__main__":
    unittest.main()
