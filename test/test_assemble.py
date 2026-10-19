#!/usr/bin/env python3

import unittest
import os
import tempfile
import importlib.util
from unittest import mock
import extern

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

FORWARD_READS = os.path.join(path_to_data, "wgsim.1.fq.gz")
REVERSE_READS = os.path.join(path_to_data, "wgsim.2.fq.gz")
LONG_READS = os.path.join(path_to_data, "pbsim.fq.gz")

class Tests(unittest.TestCase):
    def run_dryrun(self, extra_args, tmpdir):
        cmd = (
            f"GTDBTK_DATA_PATH=. "
            f"CHECKM2DB=. "
            f"EGGNOG_DATA_DIR=. "
            f"METABULI_DB_PATH=. "
            f"SINGLEM_METAPACKAGE_PATH=. "
            f"aviary assemble "
            f"{extra_args} "
            f"--output {tmpdir}/test "
            f"--dryrun --tmpdir {tmpdir} "
            #f"--snakemake-cmds \" --quiet\" "
        )
        return extern.run(cmd)
    def test_assemble_simple_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = self.run_dryrun(
                f"-1 {FORWARD_READS} -2 {REVERSE_READS}",
                tmpdir,
            )

            # General
            self.assertTrue("complete_assembly " not in output)
            self.assertTrue("complete_assembly_with_qc" in output)

            # Short-read
            self.assertTrue("assemble_short_reads" in output)
            self.assertTrue("complete_qc_short" in output)
            self.assertTrue("move_spades_assembly" in output)

            # Long-read
            self.assertTrue("long_read_assembly" not in output)
            self.assertTrue("polish_metagenome_flye" not in output)
            self.assertTrue("combine_long_only" not in output)
            self.assertTrue("skip_unicycler " not in output)
            self.assertTrue("skip_unicycler_with_qc" not in output)
            self.assertTrue("complete_qc_long" not in output)

            # Hybrid assembly
            self.assertTrue("generate_pilon_sort" not in output)
            self.assertTrue("polish_meta_pilon" not in output)
            self.assertTrue("polish_meta_racon_ill" not in output)
            self.assertTrue("get_high_cov_contigs" not in output)
            self.assertTrue("filter_illumina_assembly" not in output)
            self.assertTrue("\nspades_assembly " not in output)
            self.assertTrue("spades_assembly_coverage" not in output)
            self.assertTrue("metabat_binning_short" not in output)
            self.assertTrue("map_long_mega" not in output)
            self.assertTrue("pool_reads" not in output)
            self.assertTrue("get_read_pools" not in output)
            self.assertTrue("assemble_pools" not in output)
            self.assertTrue("combine_assemblies" not in output)
            self.assertTrue("reset_to_spades_assembly" not in output)
            self.assertTrue("remove_final_contigs" not in output)
            self.assertTrue("complete_qc_all" not in output)

            # QC
            self.assertTrue("qc_short_reads" in output)
            self.assertTrue("qc_long_reads" not in output)
            self.assertTrue("rastqc " in output)
            self.assertTrue("rastqc_long" not in output)
            self.assertTrue("nanoplot" not in output)
            self.assertTrue("metaquast" not in output)

    def test_assemble_simple_inputs_skip_qc(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = self.run_dryrun(
                f"-1 {FORWARD_READS} -2 {REVERSE_READS} --skip-qc",
                tmpdir,
            )

            # General
            self.assertTrue("complete_assembly " not in output)
            self.assertTrue("complete_assembly_with_qc" in output)

            # Short-read
            self.assertTrue("assemble_short_reads" in output)
            self.assertTrue("complete_qc_short" in output)
            self.assertTrue("move_spades_assembly" in output)

            # Long-read
            self.assertTrue("long_read_assembly" not in output)
            self.assertTrue("polish_metagenome_flye" not in output)
            self.assertTrue("combine_long_only" not in output)
            self.assertTrue("skip_unicycler " not in output)
            self.assertTrue("skip_unicycler_with_qc" not in output)
            self.assertTrue("complete_qc_long" not in output)

            # Hybrid assembly
            self.assertTrue("generate_pilon_sort" not in output)
            self.assertTrue("polish_meta_pilon" not in output)
            self.assertTrue("polish_meta_racon_ill" not in output)
            self.assertTrue("get_high_cov_contigs" not in output)
            self.assertTrue("filter_illumina_assembly" not in output)
            self.assertTrue("\nspades_assembly " not in output)
            self.assertTrue("spades_assembly_coverage" not in output)
            self.assertTrue("metabat_binning_short" not in output)
            self.assertTrue("map_long_mega" not in output)
            self.assertTrue("pool_reads" not in output)
            self.assertTrue("get_read_pools" not in output)
            self.assertTrue("assemble_pools" not in output)
            self.assertTrue("combine_assemblies" not in output)
            self.assertTrue("reset_to_spades_assembly" not in output)
            self.assertTrue("remove_final_contigs" not in output)
            self.assertTrue("complete_qc_all" not in output)

            # QC
            self.assertTrue("qc_short_reads" not in output)
            self.assertTrue("qc_long_reads" not in output)
            self.assertTrue("rastqc " in output)
            self.assertTrue("rastqc_long" not in output)
            self.assertTrue("nanoplot" not in output)
            self.assertTrue("metaquast" not in output)

    def test_assemble_hybrid(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = self.run_dryrun(
                f"-1 {FORWARD_READS} -2 {REVERSE_READS} -l {LONG_READS} "
                f"--longread-type ont",
                tmpdir,
            )
            print(output)

            # General
            self.assertTrue("complete_assembly " not in output)
            self.assertTrue("complete_assembly_with_qc" in output)

            # Short-read
            self.assertTrue("assemble_short_reads" not in output)
            self.assertTrue("complete_qc_short" not in output)
            self.assertTrue("move_spades_assembly" not in output)

            # Long-read
            self.assertTrue("long_read_assembly" in output)
            self.assertTrue("polish_metagenome_flye" in output)
            self.assertTrue("combine_long_only" not in output)
            self.assertTrue("skip_unicycler " not in output)
            self.assertTrue("skip_unicycler_with_qc" in output)
            self.assertTrue("complete_qc_long" not in output)

            # Hybrid assembly
            self.assertTrue("generate_pilon_sort" in output)
            self.assertTrue("polish_meta_pilon" in output)
            self.assertTrue("polish_meta_racon_ill" in output)
            self.assertTrue("get_high_cov_contigs" in output)
            self.assertTrue("filter_illumina_assembly" in output)
            self.assertTrue("\nspades_assembly " in output)
            self.assertTrue("spades_assembly_coverage" not in output)
            self.assertTrue("metabat_binning_short" not in output)
            self.assertTrue("map_long_mega" not in output)
            self.assertTrue("pool_reads" not in output)
            self.assertTrue("get_read_pools" not in output)
            self.assertTrue("assemble_pools" not in output)
            self.assertTrue("combine_assemblies" not in output)
            self.assertTrue("reset_to_spades_assembly" not in output)
            self.assertTrue("remove_final_contigs" not in output)
            self.assertTrue("complete_qc_all" in output)

            # QC
            self.assertTrue("qc_short_reads" in output)
            self.assertTrue("qc_long_reads" in output)
            self.assertTrue("rastqc " in output)
            self.assertTrue("rastqc_long" in output)
            self.assertTrue("nanoplot" in output)
            self.assertTrue("metaquast" not in output)

    def test_assemble_hybrid_skip_qc(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = self.run_dryrun(
                f"-1 {FORWARD_READS} -2 {REVERSE_READS} -l {LONG_READS} "
                f"--longread-type ont --skip-qc",
                tmpdir,
            )

            # General
            self.assertTrue("complete_assembly " not in output)
            self.assertTrue("complete_assembly_with_qc" in output)

            # Short-read
            self.assertTrue("assemble_short_reads" not in output)
            self.assertTrue("complete_qc_short" not in output)
            self.assertTrue("move_spades_assembly" not in output)

            # Long-read
            self.assertTrue("long_read_assembly" in output)
            self.assertTrue("polish_metagenome_flye" in output)
            self.assertTrue("combine_long_only" not in output)
            self.assertTrue("skip_unicycler " not in output)
            self.assertTrue("skip_unicycler_with_qc" in output)
            self.assertTrue("complete_qc_long" not in output)

            # Hybrid assembly
            self.assertTrue("generate_pilon_sort" in output)
            self.assertTrue("polish_meta_pilon" in output)
            self.assertTrue("polish_meta_racon_ill" in output)
            self.assertTrue("get_high_cov_contigs" in output)
            self.assertTrue("filter_illumina_assembly" in output)
            self.assertTrue("\nspades_assembly " in output)
            self.assertTrue("spades_assembly_coverage" not in output)
            self.assertTrue("metabat_binning_short" not in output)
            self.assertTrue("map_long_mega" not in output)
            self.assertTrue("pool_reads" not in output)
            self.assertTrue("get_read_pools" not in output)
            self.assertTrue("assemble_pools" not in output)
            self.assertTrue("combine_assemblies" not in output)
            self.assertTrue("reset_to_spades_assembly" not in output)
            self.assertTrue("remove_final_contigs" not in output)
            self.assertTrue("complete_qc_all" in output)

            # QC
            self.assertTrue("qc_short_reads" not in output)
            # Currently --skip-qc does not skip long read QC
            self.assertTrue("qc_long_reads" in output)
            self.assertTrue("rastqc " in output)
            self.assertTrue("rastqc_long" in output)
            self.assertTrue("nanoplot" in output)
            self.assertTrue("metaquast" not in output)

    def test_assemble_long_reads(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = self.run_dryrun(
                f"-l {LONG_READS} --longread-type ont",
                tmpdir,
            )

            # General
            self.assertTrue("complete_assembly " not in output)
            self.assertTrue("complete_assembly_with_qc" in output)

            # Short-read
            self.assertTrue("assemble_short_reads" not in output)
            self.assertTrue("complete_qc_short" not in output)
            self.assertTrue("move_spades_assembly" not in output)

            # Long-read
            self.assertTrue("long_read_assembly" in output)
            self.assertTrue("polish_metagenome_flye" in output)
            self.assertTrue("combine_long_only" in output)
            self.assertTrue("skip_unicycler " not in output)
            self.assertTrue("skip_unicycler_with_qc" not in output)
            self.assertTrue("complete_qc_long" in output)

            # Hybrid assembly
            self.assertTrue("generate_pilon_sort" not in output)
            self.assertTrue("polish_meta_pilon" not in output)
            self.assertTrue("polish_meta_racon_ill" not in output)
            self.assertTrue("get_high_cov_contigs" not in output)
            self.assertTrue("filter_illumina_assembly" not in output)
            self.assertTrue("\nspades_assembly " not in output)
            self.assertTrue("spades_assembly_coverage" not in output)
            self.assertTrue("metabat_binning_short" not in output)
            self.assertTrue("map_long_mega" not in output)
            self.assertTrue("pool_reads" not in output)
            self.assertTrue("get_read_pools" not in output)
            self.assertTrue("assemble_pools" not in output)
            self.assertTrue("combine_assemblies" not in output)
            self.assertTrue("reset_to_spades_assembly" not in output)
            self.assertTrue("remove_final_contigs" not in output)
            self.assertTrue("complete_qc_all" not in output)

            # QC
            self.assertTrue("qc_short_reads" not in output)
            self.assertTrue("qc_long_reads" in output)
            self.assertTrue("rastqc " not in output)
            self.assertTrue("rastqc_long" in output)
            self.assertTrue("nanoplot" in output)
            self.assertTrue("metaquast" not in output)

    def test_assemble_long_reads_skip_qc(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = self.run_dryrun(
                f"-l {LONG_READS} --longread-type ont --skip-qc",
                tmpdir,
            )

            # General
            self.assertTrue("complete_assembly " not in output)
            self.assertTrue("complete_assembly_with_qc" in output)

            # Short-read
            self.assertTrue("assemble_short_reads" not in output)
            self.assertTrue("complete_qc_short" not in output)
            self.assertTrue("move_spades_assembly" not in output)

            # Long-read
            # self.assertTrue("qc_long_reads" not in output) # This is present, but the script doesn't actually do any qc
            self.assertTrue("long_read_assembly" in output)
            self.assertTrue("polish_metagenome_flye" in output)
            self.assertTrue("combine_long_only" in output)
            self.assertTrue("skip_unicycler " not in output)
            self.assertTrue("skip_unicycler_with_qc" not in output)
            self.assertTrue("complete_qc_long" in output)

            # Hybrid assembly
            self.assertTrue("generate_pilon_sort" not in output)
            self.assertTrue("polish_meta_pilon" not in output)
            self.assertTrue("polish_meta_racon_ill" not in output)
            self.assertTrue("get_high_cov_contigs" not in output)
            self.assertTrue("filter_illumina_assembly" not in output)
            self.assertTrue("\nspades_assembly " not in output)
            self.assertTrue("spades_assembly_coverage" not in output)
            self.assertTrue("metabat_binning_short" not in output)
            self.assertTrue("map_long_mega" not in output)
            self.assertTrue("pool_reads" not in output)
            self.assertTrue("get_read_pools" not in output)
            self.assertTrue("assemble_pools" not in output)
            self.assertTrue("combine_assemblies" not in output)
            self.assertTrue("reset_to_spades_assembly" not in output)
            self.assertTrue("remove_final_contigs" not in output)
            self.assertTrue("complete_qc_all" not in output)

            # QC
            self.assertTrue("qc_short_reads" not in output)
            # Currently --skip-qc does not skip long read QC
            self.assertTrue("qc_long_reads" in output)
            self.assertTrue("rastqc " not in output)
            self.assertTrue("rastqc_long" in output)
            self.assertTrue("nanoplot" in output)
            self.assertTrue("metaquast" not in output)

    def test_assemble_long_reads_flye(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = self.run_dryrun(
                f"-l {LONG_READS} --longread-type ont --long-read-assembler flye",
                tmpdir,
            )

            self.assertTrue("long_read_assembly" in output)
            self.assertTrue("logs/flye_assembly.log" in output)
            self.assertTrue("logs/myloasm_assembly.log" not in output)

    def test_assemble_long_reads_myloasm(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = self.run_dryrun(
                f"-l {LONG_READS} --longread-type ont --long-read-assembler myloasm",
                tmpdir,
            )

            self.assertTrue("long_read_assembly" in output)
            self.assertTrue("logs/myloasm_assembly.log" in output)
            self.assertTrue("logs/flye_assembly.log" not in output)

    def test_myloasm_kmc_is_written_to_config(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            self.run_dryrun(f"-l {LONG_READS} --longread-type ont --myloasm-kmc --myloasm-kmc-ram 8", tmpdir)
            with open(f"{tmpdir}/test/config.yaml") as config:
                text = config.read()
            self.assertIn("myloasm_kmc: true", text)
            self.assertIn("myloasm_kmc_ram: 8", text)

    def test_myloasm_kmc_is_rejected_for_flye(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(extern.ExternCalledProcessError):
                self.run_dryrun(
                    f"-l {LONG_READS} --longread-type ont --long-read-assembler flye --myloasm-kmc", tmpdir)


class TestRunMyloasmCommand(unittest.TestCase):
    """run_myloasm.py turns --kmc/--kmc-ram into myloasm's own flags."""

    spec = importlib.util.spec_from_file_location(
        "run_myloasm", os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
                                    "aviary", "modules", "assembly", "scripts", "run_myloasm.py"))

    def _cmd(self, **kwargs):
        try:
            module = importlib.util.module_from_spec(self.spec)
            self.spec.loader.exec_module(module)
        except ImportError as error:  # Bio lives in the myloasm env
            self.skipTest(f"run_myloasm.py needs {error.name}")
        with tempfile.TemporaryDirectory() as tmpdir, \
                mock.patch.object(module.subprocess, "run") as run, \
                mock.patch.object(module, "write_assembly_info"):
            open(os.path.join(tmpdir, "assembly.fasta"), "w").close()
            module.run_myloasm("reads.fq.gz", tmpdir, "ont", 4, os.path.join(tmpdir, "log"), **kwargs)
            return run.call_args.args[0]

    def test_kmc_off_by_default(self):
        self.assertNotIn("--kmc", self._cmd())

    def test_kmc_with_ram_budget(self):
        cmd = self._cmd(kmc=True, kmc_ram=8)
        self.assertEqual(cmd[cmd.index("--kmc"):], ["--kmc", "--kmc-ram", "8"])

    def test_kmc_without_ram_uses_myloasm_default(self):
        self.assertNotIn("--kmc-ram", self._cmd(kmc=True))


if __name__ == '__main__':
    unittest.main()
