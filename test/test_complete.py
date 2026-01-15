#!/usr/bin/env python3

import unittest
import os
import tempfile
import extern
from snakemake.common.configfile import load_configfile

path_to_data = os.path.join(os.path.dirname(os.path.realpath(__file__)),'data')

FORWARD_READS = os.path.join(path_to_data, "wgsim.1.fq.gz")
REVERSE_READS = os.path.join(path_to_data, "wgsim.2.fq.gz")
LONG_READS = os.path.join(path_to_data, "pbsim.fq.gz")
ASSEMBLY = os.path.join(path_to_data, "assembly.fasta")

class Tests(unittest.TestCase):
    def test_complete_simple_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary complete "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            #-------------------------------------------------------------------
            # Assembly
            #-------------------------------------------------------------------
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

            #-------------------------------------------------------------------
            # Recovery
            #-------------------------------------------------------------------
            # Binners
            self.assertTrue("filter_contigs_by_size" in output)
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
            self.assertTrue("singlem_pipe_reads" in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("finalise_stats" in output)
            self.assertTrue("recover_mags" in output)

            #-------------------------------------------------------------------
            # Annotate
            #-------------------------------------------------------------------
            self.assertTrue("eggnog" in output)
            self.assertTrue("annotate" in output)

    def test_complete_assembly(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary complete "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--assembly {ASSEMBLY} "
                f"--output {tmpdir}/test "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            #-------------------------------------------------------------------
            # Assembly
            #-------------------------------------------------------------------
            # General
            self.assertTrue("complete_assembly " not in output)
            self.assertTrue("complete_assembly_with_qc" not in output)

            # Short-read
            self.assertTrue("assemble_short_reads" not in output)
            self.assertTrue("complete_qc_short" not in output)
            self.assertTrue("move_spades_assembly" not in output)

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
            self.assertTrue("fastqc " not in output)
            self.assertTrue("fastqc_long" not in output)
            self.assertTrue("nanoplot" not in output)
            self.assertTrue("metaquast" not in output)

            #-------------------------------------------------------------------
            # Recovery
            #-------------------------------------------------------------------
            # Binners
            self.assertTrue("filter_contigs_by_size" in output)
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
            self.assertTrue("singlem_pipe_reads" in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("finalise_stats" in output)
            self.assertTrue("recover_mags" in output)

            #-------------------------------------------------------------------
            # Annotate
            #-------------------------------------------------------------------
            self.assertTrue("eggnog" in output)
            self.assertTrue("annotate" in output)

    def test_complete_hybrid(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary complete "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"-l {LONG_READS} "
                f"--longread-type ont "
                f"--output {tmpdir}/test "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            #-------------------------------------------------------------------
            # Assembly
            #-------------------------------------------------------------------
            # General
            self.assertTrue("complete_assembly " not in output)
            self.assertTrue("complete_assembly_with_qc" not in output)

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

            #-------------------------------------------------------------------
            # Recovery
            #-------------------------------------------------------------------
            # Binners
            self.assertTrue("filter_contigs_by_size" in output)
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
            self.assertTrue("singlem_pipe_reads" in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("finalise_stats" in output)
            self.assertTrue("recover_mags" in output)

            #-------------------------------------------------------------------
            # Annotate
            #-------------------------------------------------------------------
            self.assertTrue("eggnog" in output)
            self.assertTrue("annotate" in output)

    def test_complete_long(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary complete "
                f"-l {LONG_READS} "
                f"--longread-type ont "
                f"--output {tmpdir}/test "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            #-------------------------------------------------------------------
            # Assembly
            #-------------------------------------------------------------------
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

            #-------------------------------------------------------------------
            # Recovery
            #-------------------------------------------------------------------
            # Binners
            self.assertTrue("filter_contigs_by_size" in output)
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
            self.assertTrue("singlem_pipe_reads" in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("finalise_stats" in output)
            self.assertTrue("recover_mags" in output)

            #-------------------------------------------------------------------
            # Annotate
            #-------------------------------------------------------------------
            self.assertTrue("eggnog" in output)
            self.assertTrue("annotate" in output)

    def test_complete_skip_binners(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary complete "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test "
                f"--skip-binners metabat "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            #-------------------------------------------------------------------
            # Assembly
            #-------------------------------------------------------------------
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

            #-------------------------------------------------------------------
            # Recovery
            #-------------------------------------------------------------------
            # Binners
            self.assertTrue("filter_contigs_by_size" in output)
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
            self.assertTrue("singlem_pipe_reads" in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("finalise_stats" in output)
            self.assertTrue("recover_mags" in output)

            #-------------------------------------------------------------------
            # Annotate
            #-------------------------------------------------------------------
            self.assertTrue("eggnog" in output)
            self.assertTrue("annotate" in output)

    def test_complete_no_singlem(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary complete "
                f"--skip-singlem "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            #-------------------------------------------------------------------
            # Assembly
            #-------------------------------------------------------------------
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

            #-------------------------------------------------------------------
            # Recovery
            #-------------------------------------------------------------------
            # Binners
            self.assertTrue("filter_contigs_by_size" in output)
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
            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("singlem_appraise" not in output)
            self.assertTrue("finalise_stats" in output)
            self.assertTrue("recover_mags" in output)

            #-------------------------------------------------------------------
            # Annotate
            #-------------------------------------------------------------------
            self.assertTrue("eggnog" in output)
            self.assertTrue("annotate" in output)

    def test_complete_no_abundances(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary complete "
                f"--skip-abundances "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            #-------------------------------------------------------------------
            # Assembly
            #-------------------------------------------------------------------
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

            #-------------------------------------------------------------------
            # Recovery
            #-------------------------------------------------------------------
            # Binners
            self.assertTrue("filter_contigs_by_size" in output)
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
            self.assertTrue("singlem_pipe_reads" in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("finalise_stats" in output)
            self.assertTrue("recover_mags" in output)

            #-------------------------------------------------------------------
            # Annotate
            #-------------------------------------------------------------------
            self.assertTrue("eggnog" in output)
            self.assertTrue("annotate" in output)

    def test_complete_no_taxonomy(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary complete "
                f"--skip-taxonomy "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            #-------------------------------------------------------------------
            # Assembly
            #-------------------------------------------------------------------
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

            #-------------------------------------------------------------------
            # Recovery
            #-------------------------------------------------------------------
            # Binners
            self.assertTrue("filter_contigs_by_size" in output)
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
            self.assertTrue("singlem_pipe_reads" in output)
            self.assertTrue("singlem_appraise" in output)
            self.assertTrue("finalise_stats" in output)
            self.assertTrue("recover_mags" in output)

            #-------------------------------------------------------------------
            # Annotate
            #-------------------------------------------------------------------
            self.assertTrue("eggnog" in output)
            self.assertTrue("annotate" in output)

    def test_complete_binning_only(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary complete "
                f"--binning-only "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test "
                f"--dryrun --tmpdir {tmpdir} "
                f"--snakemake-cmds \" --quiet\" "
            )
            output = extern.run(cmd)

            #-------------------------------------------------------------------
            # Assembly
            #-------------------------------------------------------------------
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

            #-------------------------------------------------------------------
            # Recovery
            #-------------------------------------------------------------------
            # Binners
            self.assertTrue("filter_contigs_by_size" in output)
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
            self.assertTrue("singlem_pipe_reads" not in output)
            self.assertTrue("singlem_appraise" not in output)
            self.assertTrue("finalise_stats" in output)
            self.assertTrue("recover_mags" in output)

            #-------------------------------------------------------------------
            # Annotate
            #-------------------------------------------------------------------
            self.assertTrue("eggnog" in output)
            self.assertTrue("annotate" in output)

    def test_complete_config(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary complete "
                f"--refinery-max-iterations 3 "
                f"--max-threads 8 "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test --tmpdir {tmpdir} "
                f"--dryrun "
                f"--snakemake-cmds \" --quiet\" "
            )
            extern.run(cmd)

            config_path = os.path.join(tmpdir, "test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)

            self.assertEqual(config["refinery_max_iterations"], 3)
            self.assertEqual(config["pplacer_threads"], 8)

    def test_complete_config_many_threads(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary complete "
                f"--max-threads 128 "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test --tmpdir {tmpdir} "
                f"--dryrun "
                f"--snakemake-cmds \" --quiet\" "
            )
            extern.run(cmd)

            config_path = os.path.join(tmpdir, "test", "config.yaml")
            self.assertTrue(os.path.exists(config_path))
            config = load_configfile(config_path)

            self.assertEqual(config["refinery_max_iterations"], 5)
            self.assertEqual(config["pplacer_threads"], 8)

    def test_complete_config_many_pplacer_threads(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cmd = (
                f"GTDBTK_DATA_PATH=. "
                f"CHECKM2DB=. "
                f"EGGNOG_DATA_DIR=. "
                f"METABULI_DB_PATH=. "
                f"SINGLEM_METAPACKAGE_PATH=. "
                f"aviary complete "
                f"--max-threads 128 "
                f"--pplacer-threads 32 "
                f"-1 {FORWARD_READS} "
                f"-2 {REVERSE_READS} "
                f"--output {tmpdir}/test --tmpdir {tmpdir} "
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
