#!/usr/bin/env python3

import unittest
import os
import tempfile
import extern

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
path_to_conda = os.path.join(path_to_data,'.conda')

FORWARD_READS = os.path.join(path_to_data, "wgsim.1.fq.gz")
REVERSE_READS = os.path.join(path_to_data, "wgsim.2.fq.gz")

class Tests(unittest.TestCase):
    def test_assemble_simple_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"aviary assemble "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("fastqc" in output)
            self.assertTrue("complete_qc_short" in output)
            self.assertTrue("assemble_short_reads" in output)
            self.assertTrue("complete_assembly_with_qc" in output)

            self.assertTrue("filter_illumina_ref" not in output)
            self.assertTrue("flye_assembly" not in output)

if __name__ == '__main__':
    unittest.main()