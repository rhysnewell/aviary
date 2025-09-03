#!/usr/bin/env python3

import unittest
import os
import tempfile
import extern

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

FORWARD_READS = os.path.join(path_to_data, "wgsim.1.fq.gz")
REVERSE_READS = os.path.join(path_to_data, "wgsim.2.fq.gz")
LONG_READS = os.path.join(path_to_data, "pbsim.fq.gz")

class Tests(unittest.TestCase):
    def test_assemble_simple_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary assemble "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            # General
            self.assertTrue("complete_assembly " not in output)
            self.assertTrue("complete_assembly_with_qc" in output)

            # Short-read
            self.assertTrue("assemble_short_reads" in output)
            self.assertTrue("complete_qc_short" in output)
            self.assertTrue("move_spades_assembly" in output)

            # Long-read
            self.assertTrue("flye_assembly" not in output)
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
            self.assertTrue("fastqc " in output)
            self.assertTrue("fastqc_long" not in output)
            self.assertTrue("nanoplot" not in output)
            self.assertTrue("metaquast" not in output)

    def test_assemble_simple_inputs_skip_qc(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary assemble "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--skip-qc "
                f"--output {tmpdir}/test "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            # General
            self.assertTrue("complete_assembly " not in output)
            self.assertTrue("complete_assembly_with_qc" in output)

            # Short-read
            self.assertTrue("assemble_short_reads" in output)
            self.assertTrue("complete_qc_short" in output)
            self.assertTrue("move_spades_assembly" in output)

            # Long-read
            self.assertTrue("flye_assembly" not in output)
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
            self.assertTrue("fastqc " in output)
            self.assertTrue("fastqc_long" not in output)
            self.assertTrue("nanoplot" not in output)
            self.assertTrue("metaquast" not in output)

    def test_assemble_hybrid(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary assemble "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"-l {LONG_READS} "
                f"--longread-type ont "
                f"--output {tmpdir}/test "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            # General
            self.assertTrue("complete_assembly " not in output)
            self.assertTrue("complete_assembly_with_qc" in output)

            # Short-read
            self.assertTrue("assemble_short_reads" not in output)
            self.assertTrue("complete_qc_short" not in output)
            self.assertTrue("move_spades_assembly" not in output)

            # Long-read
            self.assertTrue("flye_assembly" in output)
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
            self.assertTrue("fastqc " in output)
            self.assertTrue("fastqc_long" in output)
            self.assertTrue("nanoplot" in output)
            self.assertTrue("metaquast" not in output)

    def test_assemble_hybrid_skip_qc(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary assemble "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"-l {LONG_READS} "
                f"--longread-type ont "
                f"--skip-qc "
                f"--output {tmpdir}/test "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            # General
            self.assertTrue("complete_assembly " not in output)
            self.assertTrue("complete_assembly_with_qc" in output)

            # Short-read
            self.assertTrue("assemble_short_reads" not in output)
            self.assertTrue("complete_qc_short" not in output)
            self.assertTrue("move_spades_assembly" not in output)

            # Long-read
            self.assertTrue("flye_assembly" in output)
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
            self.assertTrue("fastqc " in output)
            self.assertTrue("fastqc_long" in output)
            self.assertTrue("nanoplot" in output)
            self.assertTrue("metaquast" not in output)

    def test_assemble_long_reads(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary assemble "
                f"-l {LONG_READS} "
                f"--longread-type ont "
                f"--output {tmpdir}/test "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            # General
            self.assertTrue("complete_assembly " not in output)
            self.assertTrue("complete_assembly_with_qc" in output)

            # Short-read
            self.assertTrue("assemble_short_reads" not in output)
            self.assertTrue("complete_qc_short" not in output)
            self.assertTrue("move_spades_assembly" not in output)

            # Long-read
            self.assertTrue("flye_assembly" in output)
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
            self.assertTrue("fastqc " not in output)
            self.assertTrue("fastqc_long" in output)
            self.assertTrue("nanoplot" in output)
            self.assertTrue("metaquast" not in output)

    def test_assemble_long_reads_skip_qc(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary assemble "
                f"-l {LONG_READS} "
                f"--longread-type ont "
                f"--skip-qc "
                f"--output {tmpdir}/test "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            # General
            self.assertTrue("complete_assembly " not in output)
            self.assertTrue("complete_assembly_with_qc" in output)

            # Short-read
            self.assertTrue("assemble_short_reads" not in output)
            self.assertTrue("complete_qc_short" not in output)
            self.assertTrue("move_spades_assembly" not in output)

            # Long-read
            self.assertTrue("flye_assembly" in output)
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
            self.assertTrue("fastqc " not in output)
            self.assertTrue("fastqc_long" in output)
            self.assertTrue("nanoplot" in output)
            self.assertTrue("metaquast" not in output)

if __name__ == '__main__':
    unittest.main()
