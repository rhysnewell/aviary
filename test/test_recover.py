#!/usr/bin/env python3

import unittest
import os
import shutil
import subprocess
import tempfile
import extern
from snakemake.common.configfile import load_configfile

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

FORWARD_READS = os.path.join(path_to_data, "wgsim.1.fq.gz")
REVERSE_READS = os.path.join(path_to_data, "wgsim.2.fq.gz")
ASSEMBLY = os.path.join(path_to_data, "assembly.fasta")

class Tests(unittest.TestCase):
    def test_read_permissions(self):
        """Test with a non-readable fastq input fastq file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a file with no read permissions
            no_read_permissions = os.path.join(tmpdir, "no_read_permissions.fq.gz")
            with open(no_read_permissions, 'w') as f:
                f.write(">abc\nACGT\n")
            os.chmod(no_read_permissions, 0o200)

            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--assembly {ASSEMBLY} "
                f"-1 {no_read_permissions} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test "
            )
            try:
                output = extern.run(cmd)
                # If the command succeeds, it should not
                # raise an exception, so we fail the test.
                raise AssertionError("Command should have failed due to read permissions error.")
            except extern.ExternCalledProcessError as e:
                output = str(e)
                self.assertTrue(" Please check permissions." in output)

            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} "
                f"-2 {no_read_permissions} "
                f"--output {tmpdir}/test "
            )
            try:
                output = extern.run(cmd)
                # If the command succeeds, it should not
                # raise an exception, so we fail the test.
                raise AssertionError("Command should have failed due to read permissions error.")
            except extern.ExternCalledProcessError as e:
                output = str(e)
                self.assertTrue(" Please check permissions." in output)

            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary assemble "
                f"--interleaved {no_read_permissions} "
                f"--output {tmpdir}/test "
            )
            try:
                output = extern.run(cmd)
                # If the command succeeds, it should not
                # raise an exception, so we fail the test.
                raise AssertionError("Command should have failed due to read permissions error.")
            except extern.ExternCalledProcessError as e:
                output = str(e)
                self.assertTrue(" Please check permissions." in output)

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
                f"--dryrun --tmpdir {tmpdir} "
                #f"--snakemake-cmds \" --quiet\" "
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
            self.assertTrue("gtdbtk" in output)
            self.assertTrue("get_abundances" in output)
            self.assertFalse("singlem_pipe_reads" in output)
            self.assertFalse("singlem_appraise" in output)
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
                f"--skip-binners metabat "
                f"--dryrun --tmpdir {tmpdir} "
                #f"--snakemake-cmds \" --quiet\" "
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
            self.assertTrue("gtdbtk" in output)
            self.assertTrue("get_abundances" in output)
            self.assertFalse("singlem_pipe_reads" in output)
            self.assertFalse("singlem_appraise" in output)
            self.assertTrue("finalise_stats" in output)
            self.assertTrue("recover_mags" in output)

            # Unnecessary
            self.assertTrue("complete_assembly_with_qc" not in output)

    def test_recover_strobealign_aemb_auto_skips_bam_only_binners(self):
        """--short-read-mapper strobealign-aemb never produces a BAM file, so
        BAM-only binners (semibin/concoct/comebin/quickbin/metacat) must be
        auto-skipped with a warning rather than failing deep in the DAG.
        Coverage-table binners (rosella/metabat2/vamb) are unaffected."""
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
                f"--short-read-mapper strobealign-aemb "
                f"--output {tmpdir}/test "
                f"--dryrun --tmpdir {tmpdir} "
            )
            result = subprocess.run(
                ["bash", "-o", "pipefail", "-c", cmd],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(result.returncode, 0)
            stdout = result.stdout.decode()
            stderr = result.stderr.decode()

            self.assertIn("strobealign-aemb produces no BAM file", stderr)
            self.assertIn("semibin", stderr)
            self.assertIn("metacat", stderr)

            # concoct/comebin/quickbin are already skipped by default (see
            # test_recover_skip_binners), so semibin and metacat (a default
            # binner) are the rules this auto-skip newly removes from the DAG here.
            self.assertNotRegex(stdout, r"rule semibin:")
            self.assertNotRegex(stdout, r"rule metacat:")
            self.assertNotRegex(stdout, r"rule concoct:")
            self.assertNotRegex(stdout, r"rule comebin:")
            self.assertNotRegex(stdout, r"rule quickbin:")
            self.assertRegex(stdout, r"rule rosella:")
            self.assertRegex(stdout, r"rule metabat2:")
            self.assertRegex(stdout, r"rule vamb:")

            config_path = os.path.join(tmpdir, "test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)
            for binner in ("semibin", "concoct", "comebin", "quickbin"):
                self.assertIn(binner, config["skip_binners"])

    def test_recover_quickbin(self):
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
                f"--extra-binners quickbin "
                f"--output {tmpdir}/test "
                f"--dryrun --tmpdir {tmpdir} "
                # f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            self.assertTrue("quickbin" in output)

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
                f"--dryrun --tmpdir {tmpdir} "
                #f"--snakemake-cmds \" --quiet\" "
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
            self.assertTrue("gtdbtk" in output)
            self.assertTrue("get_abundances" in output)
            self.assertFalse("singlem_pipe_reads" in output)
            self.assertFalse("singlem_appraise" in output)
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
                f"--dryrun --tmpdir {tmpdir} "
                #f"--snakemake-cmds \" --quiet\" "
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
            self.assertTrue("gtdbtk" in output)
            self.assertTrue("get_abundances" not in output)
            self.assertFalse("singlem_pipe_reads" in output)
            self.assertFalse("singlem_appraise" in output)
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
                f"--dryrun --tmpdir {tmpdir} "
                #f"--snakemake-cmds \" --quiet\" "
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
            self.assertTrue("gtdbtk" not in output)
            self.assertTrue("get_abundances" in output)
            self.assertFalse("singlem_pipe_reads" in output)
            self.assertFalse("singlem_appraise" in output)
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
                f"--dryrun --tmpdir {tmpdir} "
                #f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)
            print(output)

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
            self.assertTrue("gtdbtk" not in output)
            self.assertTrue("get_abundances" not in output)
            self.assertFalse("singlem_pipe_reads" in output)
            self.assertFalse("singlem_appraise" in output)
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
                f"--min-completeness 62 "
                f"--max-contamination 3 "
                f"--max-threads 8 "
                f"--assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test --tmpdir {tmpdir} "
                f"--dryrun "
                #f"--snakemake-cmds \" --quiet\" "
            )
            extern.run(cmd)

            config_path = os.path.join(tmpdir, "test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)

            self.assertEqual(config["refinery_max_iterations"], 3)
            self.assertEqual(config["pplacer_threads"], 8)
            self.assertEqual(config["filter_bins_min_completeness"], 62.0)
            self.assertEqual(config["filter_bins_max_contamination"], 3.0)
            self.assertEqual(config["semibin_mode"], "single")

    def test_recover_semibin_mode_config(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--semibin-mode multi "
                f"--assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test --tmpdir {tmpdir} "
                f"--dryrun "
                #f"--snakemake-cmds \" --quiet\" "
            )
            extern.run(cmd)

            config_path = os.path.join(tmpdir, "test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)

            self.assertEqual(config["semibin_mode"], "multi")

    def test_recover_semibin_multi_single_assembly_warning(self):
        """A single assembly with --semibin-mode multi should warn but not fail."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--semibin-mode multi "
                f"--assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test --tmpdir {tmpdir} "
                f"--dryrun"
            )
            result = subprocess.run(
                ["bash", "-o", "pipefail", "-c", cmd],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            # Run must succeed — warning does not abort the pipeline
            self.assertEqual(result.returncode, 0)
            stderr = result.stderr.decode()
            self.assertIn("--semibin-mode multi", stderr)
            self.assertIn("only one assembly", stderr)

    def test_recover_multiple_assemblies_default_mode_fails_clearly(self):
        """Multiple assemblies without --semibin-mode multi must fail loudly,
        not silently bin only one assembly or silently ignore the others.

        Uses two distinct file paths (processor.py dedupes --assembly by
        os.path.abspath, so passing the same path twice would collapse to a
        single assembly and never trigger this check at all).
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            assembly_2 = os.path.join(tmpdir, "assembly_2.fasta")
            shutil.copyfile(ASSEMBLY, assembly_2)

            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--assembly {ASSEMBLY} {assembly_2} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test --tmpdir {tmpdir} "
                f"--dryrun"
            )
            result = subprocess.run(
                ["bash", "-o", "pipefail", "-c", cmd],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertNotEqual(result.returncode, 0)
            stderr = result.stderr.decode()
            self.assertIn("Multiple assemblies provided", stderr)
            self.assertIn("--semibin-mode", stderr)

    def test_recover_semibin_multi_no_assembly_fails_clearly(self):
        """--semibin-mode multi with no --assembly (reads-only) must fail
        loudly rather than reaching SemiBin2's concatenate_fasta with an
        empty input list."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--semibin-mode multi "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test --tmpdir {tmpdir} "
                f"--dryrun"
            )
            result = subprocess.run(
                ["bash", "-o", "pipefail", "-c", cmd],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertNotEqual(result.returncode, 0)
            stderr = result.stderr.decode()
            self.assertIn("--semibin-mode multi requires at least two assemblies", stderr)

    def test_recover_semibin_multi_dag_wiring(self):
        """Two assemblies + --semibin-mode multi must select the multi-sample
        rules/subcommand in the DAG, not just set the config value."""
        with tempfile.TemporaryDirectory() as tmpdir:
            assembly_2 = os.path.join(tmpdir, "assembly_2.fasta")
            shutil.copyfile(ASSEMBLY, assembly_2)

            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--semibin-mode multi "
                f"--assembly {ASSEMBLY} {assembly_2} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test --tmpdir {tmpdir} "
                f"--dryrun"
            )
            result = subprocess.run(
                ["bash", "-o", "pipefail", "-c", cmd],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(result.returncode, 0)
            output = result.stdout.decode() + result.stderr.decode()
            self.assertIn("semibin_multi_prepare", output)
            self.assertIn("semibin_multi_bams", output)

            config_path = os.path.join(tmpdir, "test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)
            self.assertEqual(config["semibin_mode"], "multi")
            self.assertEqual(len(config["fasta"]), 2)

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
                f"--dryrun "
                #f"--snakemake-cmds \" --quiet\" "
            )
            extern.run(cmd)

            config_path = os.path.join(tmpdir, "test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)

            self.assertEqual(config["refinery_max_iterations"], 5)
            self.assertEqual(config["pplacer_threads"], 8)
            self.assertEqual(config["filter_bins_min_completeness"], 50.0)
            self.assertEqual(config["filter_bins_max_contamination"], 5.0)

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
                f"--dryrun "
                #f"--snakemake-cmds \" --quiet\" "
            )
            extern.run(cmd)

            config_path = os.path.join(tmpdir, "test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)

            self.assertEqual(config["refinery_max_iterations"], 5)
            self.assertEqual(config["pplacer_threads"], 32)
            self.assertEqual(config["filter_bins_min_completeness"], 50.0)
            self.assertEqual(config["filter_bins_max_contamination"], 5.0)

    def test_mapper_config_defaults(self):
        # --short-read-mapper/--long-read-mapper default to None in argparse so
        # that "the user chose this" is distinguishable from "unset"; these pin
        # the resolved values that actually reach the workflow, which must stay
        # strobealign/rammap. The two model keys are written as the string
        # "none" rather than omitted, because the .smk rules interpolate them
        # unconditionally and a missing key is a KeyError at parse time.
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. CHECKM2DB=. EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test --tmpdir {tmpdir} --dryrun "
            )
            extern.run(cmd)
            config = load_configfile(os.path.join(tmpdir, "test", "config.yaml"))

            self.assertEqual(config["short_read_mapper"], "strobealign")
            self.assertEqual(config["long_read_mapper"], "rammap")
            self.assertEqual(config["long_read_mapper_model"], "none")
            self.assertEqual(config["minibwa_params"], "none")

    def test_mapper_config_with_explicit_selection(self):
        # bwa-mem resolves to a CoverM -p value verbatim; a --short-read-mapper
        # -model is folded into short_read_mapper rather than kept separate,
        # since that is the form CoverM and polish.py both consume.
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. CHECKM2DB=. EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--short-read-mapper minimap2 --short-read-mapper-model sr "
                f"--output {tmpdir}/test --tmpdir {tmpdir} --dryrun "
            )
            extern.run(cmd)
            config = load_configfile(os.path.join(tmpdir, "test", "config.yaml"))

            self.assertEqual(config["short_read_mapper"], "minimap2-sr")

    def test_minibwa_params_reach_the_config_verbatim(self):
        # The value starts with '-', so it has to survive argparse here and
        # shell-quoting in binning.smk/qc.smk on the way back out.
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. CHECKM2DB=. EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover "
                f"--assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--short-read-mapper minibwa --minibwa-params '-x lr' "
                f"--output {tmpdir}/test --tmpdir {tmpdir} --dryrun "
            )
            extern.run(cmd)
            config = load_configfile(os.path.join(tmpdir, "test", "config.yaml"))

            self.assertEqual(config["minibwa_params"], "-x lr")

    def test_metacat_options_reach_recovery_config(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. CHECKM2DB=. EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover --assembly {ASSEMBLY} "
                f"-1 {FORWARD_READS} -2 {REVERSE_READS} "
                "--metacat-coverage-params='--mapq 15' "
                "--metacat-cluster-params='--min-cluster 10' "
                f"--output {tmpdir}/test --tmpdir {tmpdir} --dryrun "
            )
            extern.run(cmd)
            config = load_configfile(os.path.join(tmpdir, "test", "config.yaml"))
            self.assertEqual(config["metacat_options"], {
                "coverage": "--mapq 15", "seed": "", "cluster": "--min-cluster 10"})
            self.assertNotIn("metacat", config["skip_binners"])

    def test_metacat_runs_by_default_and_can_be_skipped(self):
        for extra, expect_skipped in (("", False), ("--skip-binners metacat", True)):
            with tempfile.TemporaryDirectory() as tmpdir:
                extern.run(
                    f"GTDBTK_DATA_PATH=. CHECKM2DB=. EGGNOG_DATA_DIR=. "
                    f"METABULI_DB_PATH=. SINGLEM_METAPACKAGE_PATH=. "
                    f"aviary recover --assembly {ASSEMBLY} -1 {FORWARD_READS} -2 {REVERSE_READS} "
                    f"{extra} --output {tmpdir}/test --tmpdir {tmpdir} --dryrun")
                config = load_configfile(os.path.join(tmpdir, "test", "config.yaml"))
                self.assertEqual("metacat" in config["skip_binners"], expect_skipped)

    def test_metacat_options_with_metacat_skipped_fail_clearly(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            result = subprocess.run(
                f"GTDBTK_DATA_PATH=. CHECKM2DB=. EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover --assembly {ASSEMBLY} -1 {FORWARD_READS} -2 {REVERSE_READS} "
                f"--skip-binners metacat --metacat-cluster-params='--min-cluster 10' "
                f"--output {tmpdir}/test --tmpdir {tmpdir} --dryrun",
                shell=True, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("MetaCAT options were given but MetaCAT is skipped", result.stdout + result.stderr)

    SHORT_READS = f"-1 {FORWARD_READS} -2 {REVERSE_READS}"
    LONG_READS = f"-l {os.path.join(path_to_data, 'pbsim.fq.gz')} --longread-type ont"

    def _metacat_commands(self, extra, reads=SHORT_READS):
        """The MetaCAT commands Snakemake would run, from a dry run that prints shell commands."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output = extern.run(
                f"GTDBTK_DATA_PATH=. CHECKM2DB=. EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. SINGLEM_METAPACKAGE_PATH=. "
                f"aviary recover --assembly {ASSEMBLY} {reads} "
                f"{extra} --output {tmpdir}/test --tmpdir {tmpdir} --dryrun --snakemake-cmds ' -p'")
        return next(line for line in output.splitlines() if "MetaCAT coverage" in line)

    @staticmethod
    def _coverage_stage(command):
        return command[command.index("MetaCAT coverage"):command.index("MetaCAT seed")]

    def test_metacat_long_flag_only_for_long_read_only_runs(self):
        long_only = self._coverage_stage(self._metacat_commands("", self.LONG_READS))
        self.assertIn(" --long", long_only)
        short_only = self._coverage_stage(self._metacat_commands(""))
        self.assertNotIn("--long", short_only)
        hybrid = self._coverage_stage(self._metacat_commands("", f"{self.SHORT_READS} {self.LONG_READS}"))
        self.assertNotIn("--long", hybrid)

    def test_metacat_long_flag_is_not_duplicated_when_the_user_gives_it(self):
        for given in ("--long", "-l"):
            stage = self._coverage_stage(self._metacat_commands(
                f"--metacat-coverage-params='{given} --mapq 5'", self.LONG_READS))
            self.assertEqual(stage.split().count("--long") + stage.split().count("-l"), 1, stage)
            self.assertIn("--mapq 5", stage)

    def test_metacat_cpu_and_gpu_use_their_own_environment_and_engine(self):
        cpu = self._metacat_commands("")
        self.assertIn(" -e metacat bash", cpu)
        self.assertIn("--swdpgmm-engine cpu", cpu)
        gpu = self._metacat_commands("--request-gpu")
        self.assertIn(" -e metacat-gpu bash", gpu)
        self.assertIn("--swdpgmm-engine gpu", gpu)
        self.assertNotIn("--swdpgmm-engine cpu", gpu)

    def test_metacat_options_reach_the_config_through_complete(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            extern.run(
                f"GTDBTK_DATA_PATH=. CHECKM2DB=. EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. SINGLEM_METAPACKAGE_PATH=. "
                f"aviary complete {self.SHORT_READS} --metacat-cluster-params='--min-cluster 10' "
                f"--output {tmpdir}/test --tmpdir {tmpdir} --dryrun")
            config = load_configfile(os.path.join(tmpdir, "test", "config.yaml"))
            self.assertEqual(config["metacat_options"]["cluster"], "--min-cluster 10")
            self.assertNotIn("metacat", config["skip_binners"])

    def test_metacat_warns_about_one_coverage_mode_for_mixed_short_and_long_reads(self):
        def stderr_for(reads, extra=""):
            with tempfile.TemporaryDirectory() as tmpdir:
                return subprocess.run(
                    f"GTDBTK_DATA_PATH=. CHECKM2DB=. EGGNOG_DATA_DIR=. "
                    f"METABULI_DB_PATH=. SINGLEM_METAPACKAGE_PATH=. "
                    f"aviary recover --assembly {ASSEMBLY} {reads} {extra} "
                    f"--output {tmpdir}/test --tmpdir {tmpdir} --dryrun",
                    shell=True, capture_output=True, text=True).stderr
        warning = "MetaCAT uses one coverage mode"
        self.assertIn(warning, stderr_for(f"{self.SHORT_READS} {self.LONG_READS}"))
        self.assertNotIn(warning, stderr_for(self.SHORT_READS))
        self.assertNotIn(warning, stderr_for(f"{self.SHORT_READS} {self.LONG_READS}", "--skip-binners metacat"))

    def test_metacat_runs_all_three_stages_and_passes_options_to_each(self):
        command = self._metacat_commands(
            "--metacat-coverage-params='--mapq 15' --metacat-seed-params='--hmmsearch hmmsearch' "
            "--metacat-cluster-params='--min-cluster 10'")
        coverage, seed, cluster = (command.index(f"MetaCAT {stage}") for stage in ("coverage", "seed", "cluster"))
        self.assertLess(coverage, seed)
        self.assertLess(seed, cluster)
        self.assertIn("--mapq 15", command[coverage:seed])
        self.assertIn("--hmmsearch hmmsearch", command[seed:cluster])
        self.assertIn("--min-cluster 10", command[cluster:])
        self.assertIn("--swdpgmm-engine cpu", command[cluster:])
        self.assertIn("metacat_bins", command)

    def test_metacat_uses_aviarys_contig_and_bin_size_limits(self):
        # Like rosella, vamb, metabat and quickbin: --min-contig-size / --min-bin-size, not MetaCAT's own defaults.
        command = self._metacat_commands("--min-contig-size 2000 --min-bin-size 300000")
        coverage, seed, cluster = (command.index(f"MetaCAT {stage}") for stage in ("coverage", "seed", "cluster"))
        self.assertIn("--min-sequence-length 2000", command[coverage:seed])
        self.assertIn("--min-sequence-length 2000", command[cluster:])
        self.assertIn("--min-cluster 300000", command[cluster:])
        default = self._metacat_commands("")
        self.assertIn("--min-cluster 200000", default)
        self.assertIn("--min-sequence-length 1500", default)

    def test_user_metacat_options_come_after_aviarys_limits_so_they_win(self):
        command = self._metacat_commands("--metacat-cluster-params='--min-cluster 50000'")
        cluster = command[command.index("MetaCAT cluster"):]
        self.assertLess(cluster.index("--min-cluster 200000"), cluster.index("--min-cluster 50000"))


if __name__ == '__main__':
    unittest.main()
