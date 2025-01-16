#!/usr/bin/env python3

import unittest
import os
import tempfile
import extern
from snakemake import load_configfile

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')
path_to_conda = os.path.join(path_to_data,'.conda')

FORWARD_READS = os.path.join(path_to_data, "wgsim.1.fq.gz")
REVERSE_READS = os.path.join(path_to_data, "wgsim.2.fq.gz")
ASSEMBLY = os.path.join(path_to_data, "assembly.fasta")

class Tests(unittest.TestCase):
    def test_recover_simple_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            # Binners
            self.assertTrue("prepare_binning_files" in output)
            self.assertTrue("get_bam_indices" in output)
            self.assertTrue("metabat_sens" in output)
            self.assertTrue("metabat_spec" in output)
            self.assertTrue("metabat_ssens" in output)
            self.assertTrue("metabat_sspec" in output)
            self.assertTrue("metabat2" in output)
            self.assertTrue("maxbin2" not in output)
            self.assertTrue("rosella" in output)
            self.assertTrue("semibin" in output)
            self.assertTrue("vamb" in output)
            self.assertTrue("concoct" not in output)
            self.assertTrue("das_tool" in output)

            # Refinery
            self.assertTrue("checkm_metabat2" in output)
            self.assertTrue("refine_metabat2" in output)
            self.assertTrue("checkm_rosella" in output)
            self.assertTrue("refine_rosella" in output)
            self.assertTrue("checkm_semibin" in output)
            self.assertTrue("refine_semibin" in output)
            self.assertTrue("checkm_das_tool" in output)
            self.assertTrue("refine_dastool" in output)

            # Extras
            self.assertTrue("checkm2" in output)
            self.assertTrue("gtdbtk" in output)
            self.assertTrue("get_abundances" in output)
            self.assertTrue("singlem_pipe_reads" in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("finalise_stats" in output)
            self.assertTrue("recover_mags" in output)

            # Unnecessary
            self.assertTrue("complete_assembly_with_qc" not in output)

    def test_recover_skip_binners(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test "
                f"--conda-prefix {path_to_conda} "
                f"--skip-binners metabat "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            # Binners
            self.assertTrue("prepare_binning_files" in output)
            self.assertTrue("get_bam_indices" in output)
            self.assertTrue("metabat_sens" not in output)
            self.assertTrue("metabat_spec" not in output)
            self.assertTrue("metabat_ssens" not in output)
            self.assertTrue("metabat_sspec" not in output)
            self.assertTrue("metabat2" not in output)
            self.assertTrue("maxbin2" not in output)
            self.assertTrue("rosella" in output)
            self.assertTrue("semibin" in output)
            self.assertTrue("vamb" in output)
            self.assertTrue("concoct" not in output)
            self.assertTrue("das_tool" in output)

            # Refinery
            self.assertTrue("checkm_metabat2" not in output)
            self.assertTrue("refine_metabat2" not in output)
            self.assertTrue("checkm_rosella" in output)
            self.assertTrue("refine_rosella" in output)
            self.assertTrue("checkm_semibin" in output)
            self.assertTrue("refine_semibin" in output)
            self.assertTrue("checkm_das_tool" in output)
            self.assertTrue("refine_dastool" in output)

            # Extras
            self.assertTrue("checkm2" in output)
            self.assertTrue("gtdbtk" in output)
            self.assertTrue("get_abundances" in output)
            self.assertTrue("singlem_pipe_reads" in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("finalise_stats" in output)
            self.assertTrue("recover_mags" in output)

            # Unnecessary
            self.assertTrue("complete_assembly_with_qc" not in output)

    def test_recover_no_singlem(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--skip-singlem "
                f"--assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            # Binners
            self.assertTrue("prepare_binning_files" in output)
            self.assertTrue("get_bam_indices" in output)
            self.assertTrue("metabat_sens" in output)
            self.assertTrue("metabat_spec" in output)
            self.assertTrue("metabat_ssens" in output)
            self.assertTrue("metabat_sspec" in output)
            self.assertTrue("metabat2" in output)
            self.assertTrue("maxbin2" not in output)
            self.assertTrue("rosella" in output)
            self.assertTrue("semibin" in output)
            self.assertTrue("vamb" in output)
            self.assertTrue("concoct" not in output)
            self.assertTrue("das_tool" in output)

            # Refinery
            self.assertTrue("checkm_metabat2" in output)
            self.assertTrue("refine_metabat2" in output)
            self.assertTrue("checkm_rosella" in output)
            self.assertTrue("refine_rosella" in output)
            self.assertTrue("checkm_semibin" in output)
            self.assertTrue("refine_semibin" in output)
            self.assertTrue("checkm_das_tool" in output)
            self.assertTrue("refine_dastool" in output)

            # Extras
            self.assertTrue("checkm2" in output)
            self.assertTrue("gtdbtk" in output)
            self.assertTrue("get_abundances" in output)
            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("singlem_appraise" not in output)
            self.assertTrue("finalise_stats" in output)
            self.assertTrue("recover_mags" in output)

            # Unnecessary
            self.assertTrue("complete_assembly_with_qc" not in output)

    def test_recover_no_abundances(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--skip-abundances "
                f"--assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            # Binners
            self.assertTrue("prepare_binning_files" in output)
            self.assertTrue("get_bam_indices" in output)
            self.assertTrue("metabat_sens" in output)
            self.assertTrue("metabat_spec" in output)
            self.assertTrue("metabat_ssens" in output)
            self.assertTrue("metabat_sspec" in output)
            self.assertTrue("metabat2" in output)
            self.assertTrue("maxbin2" not in output)
            self.assertTrue("rosella" in output)
            self.assertTrue("semibin" in output)
            self.assertTrue("vamb" in output)
            self.assertTrue("concoct" not in output)
            self.assertTrue("das_tool" in output)

            # Refinery
            self.assertTrue("checkm_metabat2" in output)
            self.assertTrue("refine_metabat2" in output)
            self.assertTrue("checkm_rosella" in output)
            self.assertTrue("refine_rosella" in output)
            self.assertTrue("checkm_semibin" in output)
            self.assertTrue("refine_semibin" in output)
            self.assertTrue("checkm_das_tool" in output)
            self.assertTrue("refine_dastool" in output)

            # Extras
            self.assertTrue("checkm2" in output)
            self.assertTrue("gtdbtk" in output)
            self.assertTrue("get_abundances" not in output)
            self.assertTrue("singlem_pipe_reads" in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("finalise_stats" in output)
            self.assertTrue("recover_mags" in output)

            # Unnecessary
            self.assertTrue("complete_assembly_with_qc" not in output)

    def test_recover_no_taxonomy(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--skip-taxonomy "
                f"--assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            # Binners
            self.assertTrue("prepare_binning_files" in output)
            self.assertTrue("get_bam_indices" in output)
            self.assertTrue("metabat_sens" in output)
            self.assertTrue("metabat_spec" in output)
            self.assertTrue("metabat_ssens" in output)
            self.assertTrue("metabat_sspec" in output)
            self.assertTrue("metabat2" in output)
            self.assertTrue("maxbin2" not in output)
            self.assertTrue("rosella" in output)
            self.assertTrue("semibin" in output)
            self.assertTrue("vamb" in output)
            self.assertTrue("concoct" not in output)
            self.assertTrue("das_tool" in output)

            # Refinery
            self.assertTrue("checkm_metabat2" in output)
            self.assertTrue("refine_metabat2" in output)
            self.assertTrue("checkm_rosella" in output)
            self.assertTrue("refine_rosella" in output)
            self.assertTrue("checkm_semibin" in output)
            self.assertTrue("refine_semibin" in output)
            self.assertTrue("checkm_das_tool" in output)
            self.assertTrue("refine_dastool" in output)

            # Extras
            self.assertTrue("checkm2" in output)
            self.assertTrue("gtdbtk" not in output)
            self.assertTrue("get_abundances" in output)
            self.assertTrue("singlem_pipe_reads" in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("finalise_stats" in output)
            self.assertTrue("recover_mags" in output)

            # Unnecessary
            self.assertTrue("complete_assembly_with_qc" not in output)

    def test_recover_binning_only(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--binning-only "
                f"--assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            # Binners
            self.assertTrue("prepare_binning_files" in output)
            self.assertTrue("get_bam_indices" in output)
            self.assertTrue("metabat_sens" in output)
            self.assertTrue("metabat_spec" in output)
            self.assertTrue("metabat_ssens" in output)
            self.assertTrue("metabat_sspec" in output)
            self.assertTrue("metabat2" in output)
            self.assertTrue("maxbin2" not in output)
            self.assertTrue("rosella" in output)
            self.assertTrue("semibin" in output)
            self.assertTrue("vamb" in output)
            self.assertTrue("concoct" not in output)
            self.assertTrue("das_tool" in output)

            # Refinery
            self.assertTrue("checkm_metabat2" in output)
            self.assertTrue("refine_metabat2" in output)
            self.assertTrue("checkm_rosella" in output)
            self.assertTrue("refine_rosella" in output)
            self.assertTrue("checkm_semibin" in output)
            self.assertTrue("refine_semibin" in output)
            self.assertTrue("checkm_das_tool" in output)
            self.assertTrue("refine_dastool" in output)

            # Extras
            self.assertTrue("checkm2" in output)
            self.assertTrue("gtdbtk" not in output)
            self.assertTrue("get_abundances" not in output)
            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("singlem_appraise" not in output)
            self.assertTrue("finalise_stats" in output)
            self.assertTrue("recover_mags" in output)

            # Unnecessary
            self.assertTrue("complete_assembly_with_qc" not in output)

    def test_recover_config(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--refinery-max-iterations 3 "
                f"--max-threads 8 "
                f"--assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test --tmpdir {tmpdir} "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-cmds \" --quiet\" "
            )
            extern.run(cmd)

            config_path = os.path.join(tmpdir, "test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)

            self.assertEqual(config["refinery_max_iterations"], 3)
            self.assertEqual(config["pplacer_threads"], 8)

    def test_recover_config_many_threads(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--max-threads 128 "
                f"--assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test --tmpdir {tmpdir} "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-cmds \" --quiet\" "
            )
            extern.run(cmd)

            config_path = os.path.join(tmpdir, "test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)

            self.assertEqual(config["refinery_max_iterations"], 5)
            self.assertEqual(config["pplacer_threads"], 8)

    def test_recover_config_many_pplacer_threads(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--max-threads 128 "
                f"--pplacer-threads 32 "
                f"--assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test --tmpdir {tmpdir} "
                f"--conda-prefix {path_to_conda} "
                f"--dryrun "
                f"--snakemake-cmds \" --quiet\" "
            )
            extern.run(cmd)

            config_path = os.path.join(tmpdir, "test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)

            self.assertEqual(config["refinery_max_iterations"], 5)
            self.assertEqual(config["pplacer_threads"], 32)

if __name__ == '__main__':
    unittest.main()
