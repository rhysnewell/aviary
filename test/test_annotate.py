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

@pytest.mark.expensive
class Tests(unittest.TestCase):
    def test_default(self):
        output_dir = os.path.join("example", "test_annotate_default")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary annotate "
            f"-o {output_dir}/aviary_out "
            f"--genome-fasta-directory test/data/bins1"
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        
        #gtdbtk
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/gtdbtk/gtdbtk.bac120.summary.tsv"))
        self.assertGreater(os.path.getsize(f"{output_dir}/aviary_out/data/gtdbtk/gtdbtk.bac120.summary.tsv"), 0)
        # check there are 3 lines in the summary file
        with open(f"{output_dir}/aviary_out/data/gtdbtk/gtdbtk.bac120.summary.tsv") as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 3)

        # eggnog
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/eggnog/SemiBin_0.emapper.annotations"))
        self.assertGreater(os.path.getsize(f"{output_dir}/aviary_out/data/eggnog/SemiBin_0.emapper.annotations"), 0)
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/eggnog/SemiBin_1.emapper.annotations"))
        self.assertGreater(os.path.getsize(f"{output_dir}/aviary_out/data/eggnog/SemiBin_1.emapper.annotations"), 0)

if __name__ == "__main__":
    unittest.main()
