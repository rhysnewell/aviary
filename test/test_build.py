#!/usr/bin/env python3

import unittest
import os
import subprocess
from bird_tool_utils import in_tempdir

data = os.path.join(os.path.dirname(__file__), 'data')

class Tests(unittest.TestCase):
    def test_assembly_build_only(self):
        with in_tempdir():
            cmd = (
                f"aviary assemble "
                f"-o build_out/aviary_out "
                f"-1 {data}/wgsim.1.fq.gz "
                f"-2 {data}/wgsim.2.fq.gz "
                f"-n 1 -t 1 "
                f"--build "
            )
            subprocess.run(cmd, shell=True, check=True)

            self.assertTrue(os.path.isdir("build_out/aviary_out"))
            self.assertFalse(os.path.isfile("build_out/aviary_out/data/final_contigs.fasta"))
            self.assertFalse(os.path.islink("build_out/aviary_out/assembly/final_contigs.fasta"))


if __name__ == '__main__':
    unittest.main()
