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

import unittest
import os.path
import subprocess
import shutil

data = os.path.join(os.path.dirname(__file__), 'data')
path_to_conda = os.path.join(data,'.conda')

class Tests(unittest.TestCase):
    def setup_output_dir(self, output_dir):
        try:
            shutil.rmtree(output_dir)
        except FileNotFoundError:
            pass
        os.makedirs(output_dir)

    def test_short_read_assembly(self):
        output_dir = os.path.join("example", "test_short_read_assembly")
        self.setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--conda-prefix {path_to_conda} "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

    def test_long_read_assembly(self):
        output_dir = os.path.join("example", "test_long_read_assembly")
        self.setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--min-read-size 10 --min-mean-q 1 "
            f"--conda-prefix {path_to_conda} "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

    def test_short_read_recovery(self):
        output_dir = os.path.join("example", "test_short_read_recovery")
        self.setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--conda-prefix {path_to_conda} "
            f"-n 32 -t 32 "
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

    def test_long_read_recovery(self):
        output_dir = os.path.join("example", "test_long_read_recovery")
        self.setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--min-read-size 10 --min-mean-q 1 "
            f"--conda-prefix {path_to_conda} "
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

    def test_long_read_only_recovery(self):
        output_dir = os.path.join("example", "test_long_read_only_recovery")
        self.setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--min-read-size 10 --min-mean-q 1 "
            f"--conda-prefix {path_to_conda} "
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

    def test_short_read_recovery_fast(self):
        output_dir = os.path.join("example", "test_short_read_recovery_fast")
        self.setup_output_dir(output_dir)
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
            f"--conda-prefix {path_to_conda} "
            f"-n 32 -t 32 "
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
        self.setup_output_dir(output_dir)

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
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"--conda-prefix {path_to_conda} "
            f"-n 32 -t 32 "
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
        self.setup_output_dir(output_dir)

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
            f"--conda-prefix {path_to_conda} "
            f"-n 32 -t 32 "
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
        self.setup_output_dir(output_dir)

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
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"--conda-prefix {path_to_conda} "
            f"-n 32 -t 32 "
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
        self.setup_output_dir(output_dir)

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
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"--conda-prefix {path_to_conda} "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 2)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    # @unittest.skip("Skipping test due to queue submission")
    def test_short_read_recovery_queue_submission(self):
        output_dir = os.path.join("example", "test_short_read_recovery_queue_submission")
        self.setup_output_dir(output_dir)

        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--conda-prefix {path_to_conda} "
            f"-n 32 -t 32 --local-cores 1 "
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
        self.setup_output_dir(output_dir)

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
            f"--conda-prefix {path_to_conda} "
            f"-n 32 -t 32 --local-cores 1 "
            f"--snakemake-profile aqua --cluster-retries 0 "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines >= 3)

    def test_batch_recovery(self):
        output_dir = os.path.join("example", "test_batch_recovery")
        self.setup_output_dir(output_dir)
        cmd = (
            f"aviary batch "
            f"-o {output_dir}/aviary_out "
            f"-f {data}/example_batch.tsv "
            f"--conda-prefix {path_to_conda} "
            f"--skip-binners rosella vamb metabat "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"--min-read-size 10 --min-mean-q 1 "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/sample_1/data/final_contigs.fasta"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/sample_2/data/final_contigs.fasta"))

        bin_info_path_1 = f"{output_dir}/aviary_out/sample_1/bins/bin_info.tsv"
        bin_info_path_2 = f"{output_dir}/aviary_out/sample_2/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path_1))
        self.assertTrue(os.path.isfile(bin_info_path_2))
        with open(bin_info_path_1) as f:
            num_lines = sum(1 for _ in f)
        self.assertEqual(num_lines, 3)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out/aviary_cluster_ani_0.95"))
        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out/aviary_cluster_ani_0.97"))
        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out/aviary_cluster_ani_0.99"))

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out/aviary_cluster_ani_0.95/pangenomes"))
        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out/aviary_cluster_ani_0.97/pangenomes"))
        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out/aviary_cluster_ani_0.99/pangenomes"))


if __name__ == "__main__":
    unittest.main()
