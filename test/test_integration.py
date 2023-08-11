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
import tempfile
import os.path
import extern
import subprocess

data = os.path.join(os.path.dirname(__file__), 'data')
path_to_conda = os.path.join(data,'.conda')

class Tests(unittest.TestCase):
    def test_short_read_assembly(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"aviary assemble "
                f"-o {tmpdir}/aviary_out "
                f"-1 {data}/wgsim.1.fq.gz "
                f"-2 {data}/wgsim.2.fq.gz "
                f"--conda-prefix {path_to_conda} "
                f"-n 32 -t 32 --tmpdir {tmpdir} "
            )
            extern.run(cmd)

            self.assertTrue(os.path.isdir(f"{tmpdir}/aviary_out"))
            self.assertTrue(os.path.isfile(f"{tmpdir}/aviary_out/data/final_contigs.fasta"))
            self.assertTrue(os.path.islink(f"{tmpdir}/aviary_out/assembly/final_contigs.fasta"))

    def test_long_read_assembly(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"aviary assemble "
                f"-o {tmpdir}/aviary_out "
                f"-1 {data}/wgsim.1.fq.gz "
                f"-2 {data}/wgsim.2.fq.gz "
                f"-l {data}/pbsim.fq.gz "
                f"--longread-type ont "
                f"--conda-prefix {path_to_conda} "
                f"-n 32 -t 32 --tmpdir {tmpdir} "
            )
            extern.run(cmd)

            self.assertTrue(os.path.isdir(f"{tmpdir}/aviary_out"))
            self.assertTrue(os.path.isfile(f"{tmpdir}/aviary_out/data/final_contigs.fasta"))
            self.assertTrue(os.path.islink(f"{tmpdir}/aviary_out/assembly/final_contigs.fasta"))


    def test_short_read_recovery(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"aviary recover "
                f"-o {tmpdir}/aviary_out "
                f"-1 {data}/wgsim.1.fq.gz "
                f"-2 {data}/wgsim.2.fq.gz "
                f"--conda-prefix {path_to_conda} "
                f"-n 32 -t 32 --tmpdir {tmpdir} "
            )
            extern.run(cmd)

            self.assertTrue(os.path.isfile(f"{tmpdir}/aviary_out/bins/bin_info.tsv"))
            self.assertTrue(os.path.isfile(f"{tmpdir}/aviary_out/data/final_contigs.fasta"))
            self.assertTrue(os.path.islink(f"{tmpdir}/aviary_out/assembly/final_contigs.fasta"))


    def test_short_read_recovery_skip_binners(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"aviary recover "
                f"--assembly {data}/assembly.fasta "
                f"-o {tmpdir}/aviary_out "
                f"-1 {data}/wgsim.1.fq.gz "
                f"-2 {data}/wgsim.2.fq.gz "
                f"--skip-binners concoct rosella vamb metabat maxbin "
                f"--refinery-max-iterations 1 "
                f"--conda-prefix {path_to_conda} "
                f"-n 32 -t 32 --tmpdir {tmpdir} "
            )
            # output = subprocess.check_output(cmd, shell=True)
            extern.run(cmd)

            self.assertTrue(os.path.isfile(f"{tmpdir}/aviary_out/bins/bin_info.tsv"))
            self.assertFalse(os.path.isfile(f"{tmpdir}/aviary_out/data/final_contigs.fasta"))


if __name__ == "__main__":
    unittest.main()
