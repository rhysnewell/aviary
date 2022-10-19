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
import sys
import subprocess
import tempfile
import os.path

import extern

data = os.path.join(os.path.dirname(__file__), 'data')

class Tests(unittest.TestCase):
    def test_short_read_assembly(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = "aviary assemble -o {}/aviary_out -1 {}/wgsim.1.fq.gz -2 {}/wgsim.2.fq.gz".format(tmpdir, data, data)
            # print(cmd)
            extern.run(cmd)
            self.assertTrue(os.path.isdir("{}/aviary_out".format(tmpdir)))
            self.assertTrue(os.path.isfile("{}/aviary_out/data/final_contigs.fasta".format(tmpdir)))
            self.assertTrue(os.path.islink("{}/aviary_out/assembly/final_contigs.fasta".format(tmpdir)))


    def test_short_read_recovery(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = "aviary recover -o {}/aviary_out -1 {}/wgsim.1.fq.gz -2 {}/wgsim.2.fq.gz".format(tmpdir, data, data)
            # print(cmd)
            extern.run(cmd)
            self.assertTrue(os.path.isfile("{}/aviary_out/bins/bin_info.tsv".format(tmpdir)))
            self.assertTrue(os.path.isfile("{}/aviary_out/data/final_contigs.fasta".format(tmpdir)))
            self.assertTrue(os.path.islink("{}/aviary_out/assembly/final_contigs.fasta".format(tmpdir)))
        


if __name__ == "__main__":
	unittest.main()
