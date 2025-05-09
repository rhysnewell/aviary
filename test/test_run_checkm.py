#!/usr/bin/env python3

import unittest
import os
import tempfile
from aviary.modules.binning.scripts.run_checkm import checkm
from unittest.mock import patch, ANY
import subprocess
from pathlib import Path

def create_output(_, env, stdout, stderr):
    os.makedirs("output_folder")
    Path(os.path.join("output_folder", "quality_report.tsv")).touch()

class Tests(unittest.TestCase):
    def test_run_checkm(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            with patch.object(subprocess, "run", side_effect=create_output) as mock_subprocess:
                checkm2_db = os.path.join("checkm2_db")
                os.makedirs(checkm2_db)
                bin_folder = os.path.join("bin_folder")
                os.makedirs(bin_folder)
                Path(os.path.join(bin_folder, "bin_1.fna")).touch()

                checkm(checkm2_db, bin_folder, "fna", 1, "output_folder", "output_file", 1, "checkm.log")
                self.assertTrue(os.path.exists("output_file"))
                os.environ["CHECKM2DB"] = f"{checkm2_db}/uniref100.KO.1.dmnd"

                mock_subprocess.assert_called_once_with(
                    f"checkm2 predict -i {bin_folder}/ -x fna -o output_folder -t 1 --force".split(),
                    env=os.environ,
                    stdout=ANY,
                    stderr=subprocess.STDOUT
                    )

    def test_run_checkm_no_bins(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            with patch.object(subprocess, "run") as mock_subprocess:
                checkm2_db = os.path.join("checkm2_db")
                os.makedirs(checkm2_db)
                bin_folder = os.path.join("bin_folder")
                os.makedirs(bin_folder)

                checkm(checkm2_db, bin_folder, "fna", 1, "output_folder", "output_file", 1, "checkm.log")
                self.assertTrue(os.path.exists("output_file"))
                mock_subprocess.assert_not_called()

    def test_run_checkm_no_refinery(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            os.chdir(tmpdir)
            with patch.object(subprocess, "run") as mock_subprocess:
                checkm2_db = os.path.join("checkm2_db")
                os.makedirs(checkm2_db)
                bin_folder = os.path.join("bin_folder")
                os.makedirs(bin_folder)
                Path(os.path.join(bin_folder, "bin_1.fna")).touch()

                checkm(checkm2_db, bin_folder, "fna", 0, "output_folder", "output_file", 1, "checkm.log")
                self.assertTrue(os.path.exists("output_file"))
                mock_subprocess.assert_not_called()


if __name__ == '__main__':
    unittest.main()
