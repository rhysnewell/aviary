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
import glob
import random
import re

data = os.path.join(os.path.dirname(__file__), 'data')

if os.environ.get("TEST_REQUEST_GPU", "0") == "1":
    request_gpu = "--request-gpu"
else:
    request_gpu = ""

singlem_metapackage = os.environ.get("SINGLEM_METAPACKAGE_PATH", "")
singlem_args = f"--skip-singlem false --singlem-metapackage-path {singlem_metapackage}" if singlem_metapackage else ""

def setup_output_dir(output_dir):
    try:
        shutil.rmtree(output_dir)
    except FileNotFoundError:
        pass
    os.makedirs(output_dir)


FASTA_LINE_WIDTH = 60
INFLATED_ASSEMBLY_COPIES = 100
INFLATED_ASSEMBLY_DIVERGENCE = 0.01
INFLATED_ASSEMBLY_SEED = 42

# Which bases a given base may mutate into. Anything not listed (N, IUPAC
# ambiguity codes) is left alone rather than silently turned into a real base.
_SUBSTITUTIONS = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}


def _read_fasta(path: str) -> list[tuple[str, str]]:
    """Read a FASTA file into (header, sequence) pairs, header including '>'."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    parts: list[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(parts)))
                header = line
                parts = []
            else:
                parts.append(line)
    if header is not None:
        records.append((header, "".join(parts)))
    return records


def _mutate(sequence: str, divergence: float, rng: random.Random) -> str:
    """Substitute `divergence` of the bases in `sequence` at random positions.

    Always substitutes at least one base, so that short contigs -- where the
    rate rounds down to zero -- do not come back byte-identical to their
    source, which is the very thing this is here to avoid.
    """
    bases = bytearray(sequence, "ascii")
    if not bases:
        return sequence
    count = min(len(bases), max(1, int(len(bases) * divergence)))
    for position in rng.sample(range(len(bases)), count):
        replacements = _SUBSTITUTIONS.get(chr(bases[position]))
        if replacements:
            bases[position] = ord(rng.choice(replacements))
    return bases.decode("ascii")


def _write_record(handle, header: str, sequence: str) -> None:
    handle.write(header + "\n")
    for start in range(0, len(sequence), FASTA_LINE_WIDTH):
        handle.write(sequence[start:start + FASTA_LINE_WIDTH] + "\n")


def write_inflated_assembly(
    source: str,
    destination: str,
    copies: int = INFLATED_ASSEMBLY_COPIES,
    divergence: float = INFLATED_ASSEMBLY_DIVERGENCE,
    seed: int = INFLATED_ASSEMBLY_SEED,
) -> None:
    """Write `source` plus `copies` slightly-mutated copies of it to `destination`.

    The binning tests need far more contigs than the toy assembly provides, so
    they inflate it. The copies used to be byte-identical (only the headers
    differed), which handed COMEBin's clustering step thousands of co-located
    points. hnswlib's graph construction cannot reliably traverse those, and
    the run fails partway with "Cannot return the results in a contiguous 2D
    array. Probably ef or M is too small" -- intermittently, because COMEBin's
    embedding step is unseeded, so whether a given run lands in a traversable
    configuration is chance.

    Substituting `divergence` of the bases in each copy gives every contig its
    own kmer composition, while keeping copies of one contig similar enough to
    remain a coherent, binnable group (strain-level divergence). It also makes
    read mapping less ambiguous: reads map to the closest variant rather than
    at random across identical twins, which steadies the coverage signal.

    Seeded, so the inflated assembly is byte-identical from run to run.
    """
    records = _read_fasta(source)
    rng = random.Random(seed)
    with open(destination, "w") as handle:
        for header, sequence in records:
            _write_record(handle, header, sequence)
        for copy_index in range(copies):
            for header, sequence in records:
                _write_record(
                    handle,
                    f"{header}{copy_index}",
                    _mutate(sequence, divergence, rng),
                )


def assert_mapper_logged(testcase, output_dir, rule_log_dir, expected_flag):
    """Assert a CoverM call in this rule used the expected -p mapper.

    The coverage scripts print the full command they run, so this checks which
    aligner was actually invoked rather than inferring it from the output. That
    distinction matters: a short-read mapper handed long reads, or a mapper
    that failed to align anything, still produces a well-formed coverage table
    -- of zeroes. Asserting on the table alone would pass.
    """
    pattern = f"{output_dir}/aviary_out/logs/{rule_log_dir}/*/attempt*.log"
    matches = glob.glob(pattern)
    testcase.assertTrue(matches, f"no log files matching {pattern}")
    text = "".join(open(path).read() for path in matches)
    testcase.assertIn(
        expected_flag, text,
        f"expected {expected_flag!r} in {rule_log_dir} log, got:\n{text[:2000]}")


def assert_coverage_has_depths(testcase, path):
    """Assert a CoverM contig-coverage table holds real, non-zero depths."""
    testcase.assertTrue(os.path.isfile(path), f"missing coverage table {path}")
    with open(path) as handle:
        header = [column.strip() for column in handle.readline().split("\t")]
        depths = [float(line.split("\t")[2]) for line in handle if line.strip()]
    testcase.assertIn("totalAvgDepth", header)
    testcase.assertTrue(len(depths) > 1, f"only {len(depths)} contigs in {path}")
    testcase.assertTrue(any(depth > 0 for depth in depths),
                        f"all depths zero in {path} -- did any reads map?")

# We have a separate class for qsub tests (tests that run aviary with the CMR
# aqua snakemake profile) so that there is no need to run other expensive tests
# when running qsub tests.
@pytest.mark.qsub
class TestsQsub(unittest.TestCase):
    def test_short_read_recovery_queue_submission(self):
        output_dir = os.path.join("example", "test_short_read_recovery_queue_submission")
        setup_output_dir(output_dir)

        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"-n 32 -t 32 --local-cores 1 "
            f"--strict "
            f"--snakemake-profile aqua --cluster-retries 3 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertEqual(num_lines, 3)

    def test_short_read_recovery_queue_submission_gpus(self):
        output_dir = os.path.join("example", "test_short_read_recovery_queue_submission_gpus")
        setup_output_dir(output_dir)

        # Create inflated assembly file
        write_inflated_assembly(f"{data}/assembly.fasta", f"{output_dir}/assembly.fasta")

        cmd = (
            f"aviary recover "
            f"--assembly {output_dir}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--skip-binners rosella metabat vamb "
            f"--extra-binners taxvamb comebin "
            f"--request-gpu "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            # The inflated/duplicated toy assembly this test bins against never
            # produces bins above ~48% completeness (no refinement is run, per
            # --refinery-max-iterations above), so the default 50%/5% filter
            # thresholds reject every bin and gtdbtk then fails on an empty
            # filtered_bins dir. Loosened here only for this smoke test.
            f"--min-completeness 20 "
            f"--max-contamination 20 "
            f"-n 32 -t 32 --local-cores 1 "
            f"--strict "
            # Retries because the GPU binners here can fail transiently on this
            # synthetic assembly: COMEBin's embeddings are unseeded, and a run
            # whose embedding space happens to hold tight duplicate clusters
            # trips hnswlib ("Cannot return the results in a contigious 2D
            # array"). write_inflated_assembly's divergence makes that far less
            # likely; this rides out the remainder rather than failing the whole
            # suite. The binner rules wipe their output dir first, so a retry is
            # a clean re-run. Matches the non-GPU sibling test, which uses 3.
            f"--snakemake-profile aqua --cluster-retries 2 "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines >= 3)

@pytest.mark.expensive
class Tests(unittest.TestCase):
    def test_short_read_assembly(self):
        output_dir = os.path.join("example", "test_short_read_assembly")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

        with open(f"{output_dir}/aviary_out/data/final_contigs.fasta") as f:
            contigs = [c for c in f.read().strip().split('\n') if not c.startswith('>')]
            total_bp = sum(len(c) for c in contigs)

        self.assertTrue(total_bp > 1500000, "Assembly should be at least 1.5 million bp without host filtering")

    def test_short_read_assembly_host(self):
        output_dir = os.path.join("example", "test_short_read_assembly_host")
        setup_output_dir(output_dir)

        cmd = f"zcat {data}/GCA_000503915.1_ASM50391v1_genomic.fna.gz > {output_dir}/host_filter.fasta"
        subprocess.run(cmd, shell=True, check=True)

        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--host-filter {output_dir}/host_filter.fasta "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

        with open(f"{output_dir}/aviary_out/data/final_contigs.fasta") as f:
            contigs = [c for c in f.read().strip().split('\n') if not c.startswith('>')]
            total_bp = sum(len(c) for c in contigs)

        self.assertTrue(total_bp < 1000000, "Assembly should be smaller than 1 million bp after host filtering")

    def test_short_read_coassembly(self):
        output_dir = os.path.join("example", "test_short_read_coassembly")
        setup_output_dir(output_dir)

        cmd = f"cp {data}/wgsim.1.fq.gz {output_dir}/wgsimagain.1.fq.gz"
        cmd2 = f"cp {data}/wgsim.2.fq.gz {output_dir}/wgsimagain.2.fq.gz"
        subprocess.run(cmd + " && " + cmd2, shell=True, check=True)

        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz {output_dir}/wgsimagain.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz {output_dir}/wgsimagain.2.fq.gz "
            f"--coassemble yes "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

    def test_short_read_coassembly_skip_qc(self):
        output_dir = os.path.join("example", "test_short_read_coassembly_skip_qc")
        setup_output_dir(output_dir)

        cmd = f"cp {data}/wgsim.1.fq.gz {output_dir}/wgsimagain.1.fq.gz"
        cmd2 = f"cp {data}/wgsim.2.fq.gz {output_dir}/wgsimagain.2.fq.gz"
        subprocess.run(cmd + " && " + cmd2, shell=True, check=True)

        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz {output_dir}/wgsimagain.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz {output_dir}/wgsimagain.2.fq.gz "
            f"--coassemble yes --skip-qc "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

    def test_long_read_assembly_default(self):
        output_dir = os.path.join("example", "test_long_read_assembly")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--min-read-size 10 --min-mean-q 1 "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

    def test_long_read_assembly_flye(self):
        output_dir = os.path.join("example", "test_long_read_assembly_flye")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--long-read-assembler flye "
            f"--min-read-size 10 --min-mean-q 1 "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

    def test_long_read_assembly_no_short_reads(self):
        output_dir = os.path.join("example", "test_long_read_assembly_no_short_reads")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--min-read-size 10 --min-mean-q 1 "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

    def test_short_read_recovery(self):
        output_dir = os.path.join("example", "test_short_read_recovery")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"{singlem_args} "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 1)

        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

        if singlem_metapackage:
            self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/diversity"))
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraise.svg"))
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/benchmarks/singlem_appraise.benchmark.txt"))

    def test_long_read_recovery_split(self):
        output_dir = os.path.join("example", "test_long_read_recovery_split")
        setup_output_dir(output_dir)

        for i, size in enumerate([80000, 50000, 20000]):
            for end in [1, 2]:
                cmd = f"zcat {data}/wgsim.{end}.fq.gz | head -n {size} | gzip > {output_dir}/wgsim_{i}.{end}.fq.gz"
                subprocess.run(cmd, shell=True, check=True)

        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz {output_dir}/wgsim_0.1.fq.gz {output_dir}/wgsim_1.1.fq.gz {output_dir}/wgsim_2.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz {output_dir}/wgsim_0.2.fq.gz {output_dir}/wgsim_1.2.fq.gz {output_dir}/wgsim_2.2.fq.gz "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--coassemble no "
            f"--coverage-job-strategy always "
            f"--coverage-samples-per-job 2 "
            f"--min-read-size 10 --min-mean-q 1 "
            # This tiny synthetic dataset only ever produces bins in the
            # 15-27% completeness range -- this test is checking that bins
            # within threshold get picked up, not that the (default 50%)
            # threshold itself works, so loosen it to match what the data
            # can actually produce.
            f"--min-completeness 20 --max-contamination 20 "
            f"-n 32 -t 32 "
            f"{singlem_args} "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

        if singlem_metapackage:
            self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/diversity"))
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraise.svg"))

    def test_long_read_recovery_default(self):
        output_dir = os.path.join("example", "test_long_read_recovery")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--min-read-size 10 --min-mean-q 1 "
            # This tiny synthetic dataset only ever produces bins in the
            # 15-27% completeness range -- this test is checking that bins
            # within threshold get picked up, not that the (default 50%)
            # threshold itself works, so loosen it to match what the data
            # can actually produce.
            f"--min-completeness 20 --max-contamination 20 "
            f"-n 32 -t 32 "
            f"--strict "
            f"{singlem_args} "
        )
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

        if singlem_metapackage:
            self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/diversity"))
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraise.svg"))

    def test_long_read_only_recovery(self):
        output_dir = os.path.join("example", "test_long_read_only_recovery")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--min-read-size 10 --min-mean-q 1 "
            f"{singlem_args} "
            f"-n 32 -t 32 "
            f"--strict "
        )
        print(cmd)
        subprocess.run(cmd, shell=True, check=True)

        self.assertTrue(os.path.isdir(f"{output_dir}/aviary_out"))
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

        if singlem_metapackage:
            self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/diversity"))
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraise.svg"))

    def test_short_read_recovery_fast(self):
        output_dir = os.path.join("example", "test_short_read_recovery_fast")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella vamb metabat "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertEqual(num_lines, 3)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_semibin(self):
        suffix = "_gpu" if os.environ.get("TEST_REQUEST_GPU") == "1" else ""
        output_dir = os.path.join("example", f"test_short_read_recovery_semibin{suffix}")
        setup_output_dir(output_dir)

        cmd = f"ln -sr {data}/wgsim.1.fq.gz {output_dir}/wgsim2.1.fq.gz && ln -sr {data}/wgsim.2.fq.gz {output_dir}/wgsim2.2.fq.gz"
        subprocess.run(cmd, shell=True, check=True)

        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz {output_dir}/wgsim2.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz {output_dir}/wgsim2.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella vamb metabat "
            f"{request_gpu} "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertEqual(num_lines, 3)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_semibin_multi(self):
        """Real end-to-end SemiBin2 multi_easy_bin run: two assemblies, two
        read sets, --semibin-mode multi.

        Deliberately uses two assemblies with IDENTICAL (colliding) contig
        names -- a copy of assembly.fasta, not a pre-uniquified variant --
        because that is what a real co-assembly scenario looks like (contigs
        named NODE_1, k141_1, etc. repeat across samples). This is the case
        that actually exercises SemiBin2 concatenate_fasta's per-sample
        prefixing and filter_contigs_by_size's prefix-stripping; a
        pre-uniquified second assembly would pass even if that name
        reconciliation were broken, and prove nothing.
        """
        suffix = "_gpu" if os.environ.get("TEST_REQUEST_GPU") == "1" else ""
        output_dir = os.path.join("example", f"test_short_read_recovery_semibin_multi{suffix}")
        setup_output_dir(output_dir)

        cmd = (
            f"cp {data}/assembly.fasta {output_dir}/assembly2.fasta && "
            f"ln -sr {data}/wgsim.1.fq.gz {output_dir}/wgsim2.1.fq.gz && "
            f"ln -sr {data}/wgsim.2.fq.gz {output_dir}/wgsim2.2.fq.gz"
        )
        subprocess.run(cmd, shell=True, check=True)

        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta {output_dir}/assembly2.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz {output_dir}/wgsim2.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz {output_dir}/wgsim2.2.fq.gz "
            f"--semibin-mode multi "
            f"--binning-only "
            f"--skip-binners rosella vamb metabat "
            f"{request_gpu} "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 1)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_vamb(self):
        output_dir = os.path.join("example", "test_short_read_recovery_vamb")
        setup_output_dir(output_dir)

        # Create inflated assembly file
        write_inflated_assembly(f"{data}/assembly.fasta", f"{output_dir}/assembly.fasta")

        cmd = (
            f"aviary recover "
            f"--assembly {output_dir}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella semibin metabat "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 2)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_rosella(self):
        output_dir = os.path.join("example", "test_short_read_recovery_rosella")
        setup_output_dir(output_dir)

        # Create inflated assembly file
        write_inflated_assembly(f"{data}/assembly.fasta", f"{output_dir}/assembly.fasta")

        cmd = (
            f"aviary recover "
            f"--assembly {output_dir}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners semibin metabat vamb "
            f"{request_gpu} "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 2)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_taxvamb(self):
        suffix = "_gpu" if os.environ.get("TEST_REQUEST_GPU") == "1" else ""
        output_dir = os.path.join("example", f"test_short_read_recovery_taxvamb{suffix}")
        logs_dir = os.path.join(output_dir, "aviary_out", "logs")
        logs_backup = logs_dir + ".bak"
        if os.path.exists(logs_dir):
            shutil.copytree(logs_dir, logs_backup)
        setup_output_dir(output_dir)
        if os.path.exists(logs_backup):
            os.makedirs(os.path.join(output_dir, "aviary_out"), exist_ok=True)
            shutil.copytree(logs_backup, logs_dir)
            shutil.rmtree(logs_backup)

        # Create inflated assembly file
        write_inflated_assembly(f"{data}/assembly.fasta", f"{output_dir}/assembly.fasta")

        cmd = (
            f"aviary recover "
            f"--assembly {output_dir}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella semibin metabat vamb "
            f"--extra-binners taxvamb "
            f"{request_gpu} "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 2)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_metabuli_unclassified(self):
        """Run the real pipeline over contigs Metabuli cannot classify.

        convert_metabuli used to crash outright on Metabuli's unclassified
        rows -- they carry a trailing tab, giving 9 tab-fields against the
        8 of classified rows, so pd.read_csv(header=None) aborted with
        "Expected 8 fields, saw 9". Any real run containing a single
        unclassified read died.

        No existing test caught it. Every metabuli-producing test
        (queue_submission_gpus, taxvamb, taxvamb_gpu) bins the same toy
        assembly, in which all 7070 contigs classify -- 0 unclassified rows --
        so they pass identically against broken and fixed code.

        This test forces the condition: the assembly is the usual test contigs
        (so binning still has something real to work with) plus synthetic
        random-sequence contigs, which are absent from Metabuli's GTDB
        database and so come back is_classified=0.

        It asserts that precondition explicitly. Without that check the test
        would silently rot into exactly the vacuous shape of the others the
        moment the database or input changed.
        """
        suffix = "_gpu" if os.environ.get("TEST_REQUEST_GPU") == "1" else ""
        output_dir = os.path.join(
            "example", f"test_short_read_recovery_metabuli_unclassified{suffix}")
        setup_output_dir(output_dir)

        assembly = os.path.join(output_dir, "assembly.fasta")

        # Inflate the real assembly the same way the taxvamb test does, so
        # there is enough signal for binning to produce bins at all.
        cmd = f"cat {data}/assembly.fasta > {assembly}"
        for i in range(100):
            cmd += (f" && awk '/^>/ {{print $0 \"{i}\"}} !/^>/ {{print $0}}' "
                    f"{data}/assembly.fasta >> {assembly}")
        subprocess.run(cmd, shell=True, check=True)

        # Append unclassifiable contigs. Fixed seed so the input -- and hence
        # whether Metabuli classifies them -- is reproducible between runs.
        rng = random.Random(1234)
        with open(assembly, "a") as fh:
            for n in range(20):
                seq = "".join(rng.choice("ACGT") for _ in range(20000))
                fh.write(f">synthetic_unclassifiable_{n}\n")
                for j in range(0, len(seq), 80):
                    fh.write(seq[j:j + 80] + "\n")

        cmd = (
            f"aviary recover "
            f"--assembly {assembly} "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella semibin metabat vamb "
            f"--extra-binners taxvamb "   # taxvamb pulls in metabuli + convert_metabuli
            f"{request_gpu} "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        # PRECONDITION: the run must actually contain unclassified rows,
        # otherwise this test proves nothing about the bug it exists for.
        classifications = os.path.join(
            output_dir, "aviary_out", "data", "metabuli_taxonomy",
            "tax_classifications.tsv")
        self.assertTrue(os.path.isfile(classifications),
                        "Metabuli produced no tax_classifications.tsv")
        with open(classifications) as f:
            rows = [l.split("\t") for l in f if not l.startswith("#")]
        unclassified = [r for r in rows if r and r[0] == "0"]
        self.assertGreater(
            len(unclassified), 0,
            "No unclassified (is_classified=0) rows were produced, so this "
            "test is not exercising the trailing-tab path it exists for. The "
            "synthetic contigs may now be classifiable, or the assembly "
            "changed.")

        # THE ACTUAL ASSERTION: conversion survived those rows. Before the fix
        # this file is missing because convert_metabuli raised ParserError.
        taxonomy = os.path.join(
            output_dir, "aviary_out", "data", "metabuli_taxonomy", "taxonomy.tsv")
        self.assertTrue(
            os.path.isfile(taxonomy),
            "convert_metabuli produced no taxonomy.tsv despite Metabuli "
            "having run -- it crashed on the unclassified rows.")

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))

    def test_short_read_recovery_comebin(self):
        suffix = "_gpu" if os.environ.get("TEST_REQUEST_GPU") == "1" else ""
        output_dir = os.path.join("example", f"test_short_read_recovery_comebin{suffix}")
        logs_dir = os.path.join(output_dir, "aviary_out", "logs")
        logs_backup = logs_dir + ".bak"
        if os.path.exists(logs_dir):
            shutil.copytree(logs_dir, logs_backup)
        setup_output_dir(output_dir)
        if os.path.exists(logs_backup):
            os.makedirs(os.path.join(output_dir, "aviary_out"), exist_ok=True)
            shutil.copytree(logs_backup, logs_dir)
            shutil.rmtree(logs_backup)

        # Create inflated assembly file
        write_inflated_assembly(f"{data}/assembly.fasta", f"{output_dir}/assembly.fasta")

        cmd = (
            f"aviary recover "
            f"--assembly {output_dir}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella semibin metabat vamb "
            f"--extra-binners comebin "
            f"{request_gpu} "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 2)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_quickbin(self):
        output_dir = os.path.join("example", "test_short_read_recovery_quickbin")
        setup_output_dir(output_dir)

        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella semibin metabat vamb "
            f"--extra-binners quickbin "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 2)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_concoct(self):
        output_dir = os.path.join("example", "test_short_read_recovery_concoct")
        setup_output_dir(output_dir)

        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella semibin metabat vamb "
            f"--extra-binners concoct "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 2)

        self.assertFalse(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))

    def test_short_read_recovery_default_mapper(self):
        # CoverM >=0.7.0 maps with strobealign, which is the path real
        # short-read runs now take. Reads are .fq.gz deliberately: strobealign
        # 0.17.0 has an mmap bug on gzipped fastq.
        output_dir = os.path.join("example", "test_short_read_recovery_default_mapper")
        setup_output_dir(output_dir)

        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella semibin vamb quickbin "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        short_cov_path = f"{output_dir}/aviary_out/data/coverm.cov"
        self.assertTrue(os.path.isfile(short_cov_path))
        with open(short_cov_path) as f:
            header = f.readline().split("\t")
            depths = [float(line.split("\t")[2]) for line in f if line.strip()]
        self.assertIn("totalAvgDepth", [c.strip() for c in header])
        # A mapper failure yields a table of zeroes rather than no table, so
        # assert reads actually mapped.
        self.assertTrue(len(depths) > 1)
        self.assertTrue(any(d > 0 for d in depths))

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 1)

    # ------------------------------------------------------------------
    # Read-mapper selection (--short-read-mapper / --long-read-mapper).
    #
    # CoverM >=0.7.0 defaults short reads to strobealign; aviary defaults long
    # reads to rammap. Both are overridable. Each test asserts the LOGGED
    # command as well as the coverage values, because the failure mode here is
    # silent: a mapper that maps nothing still emits a well-formed table.
    #
    # Coverage tests use --binning-only --skip-qc and a single binner so they
    # exercise get_coverage.py and little else.
    # ------------------------------------------------------------------

    def _run_short_read_coverage(self, name, mapper_flag, expected_p):
        output_dir = os.path.join("example", name)
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"{mapper_flag} "
            f"--binning-only "
            f"--skip-binners rosella semibin vamb quickbin "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
            # A rule (e.g. das_tool's diamond call) can stall on transient
            # storage contention when many benchmarking jobs run concurrently
            # on shared /mnt/weka. Retrying re-runs just the failed rule
            # rather than failing the whole test outright.
            f"--cluster-retries 2 "
        )
        subprocess.run(cmd, shell=True, check=True)
        assert_coverage_has_depths(self, f"{output_dir}/aviary_out/data/coverm.cov")
        assert_mapper_logged(self, output_dir, "coverm_prepare", f"-p {expected_p}")

    def _run_long_read_coverage(self, name, mapper_flag, expected_p):
        output_dir = os.path.join("example", name)
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--min-read-size 10 --min-mean-q 1 "
            f"{mapper_flag} "
            f"--binning-only "
            f"--skip-binners rosella semibin vamb quickbin "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)
        assert_coverage_has_depths(self, f"{output_dir}/aviary_out/data/coverm.cov")
        assert_mapper_logged(self, output_dir, "coverm_prepare", f"-p {expected_p}")

    def test_short_read_coverage_strobealign_default(self):
        # No flag: proves the CoverM default is what aviary actually uses.
        self._run_short_read_coverage(
            "test_short_read_coverage_strobealign_default", "", "strobealign")

    def test_short_read_coverage_minimap2(self):
        # Proves the whole chain is connected -- CLI, processor.py, config,
        # Snakefile, script. If any link broke, the default would silently win.
        self._run_short_read_coverage(
            "test_short_read_coverage_minimap2",
            "--short-read-mapper minimap2", "minimap2-sr")

    def test_short_read_coverage_rammap(self):
        self._run_short_read_coverage(
            "test_short_read_coverage_rammap",
            "--short-read-mapper rammap", "rammap-sr")

    def test_short_read_coverage_minibwa(self):
        self._run_short_read_coverage(
            "test_short_read_coverage_minibwa",
            "--short-read-mapper minibwa", "minibwa")

    def test_short_read_coverage_bwa_mem(self):
        # bwa-mem/bwa-mem2 are the only short-read mappers CoverM does not pull
        # in itself, so this also proves the binaries are actually present in
        # the coverm pixi env rather than merely declared.
        self._run_short_read_coverage(
            "test_short_read_coverage_bwa_mem",
            "--short-read-mapper bwa-mem", "bwa-mem")

    def test_short_read_coverage_bwa_mem2(self):
        self._run_short_read_coverage(
            "test_short_read_coverage_bwa_mem2",
            "--short-read-mapper bwa-mem2", "bwa-mem2")

    def test_short_read_polishing_with_bwa_mem(self):
        # The one path where bwa-mem is not just another CoverM -p value: racon
        # needs PAF, bwa mem only emits SAM, so polish.py indexes the reference
        # and pipes bwa mem through paftools.js sam2paf. Nothing else covers
        # that pipeline end to end, and it runs in the polishing env, which is a
        # different env from the coverage tests above.
        output_dir = os.path.join("example", "test_short_read_polishing_bwa_mem")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-l {data}/pbsim.fq.gz "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--longread-type ccs "
            f"--short-read-mapper bwa-mem "
            f"--min-read-size 10 --min-mean-q 1 "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)
        self.assertTrue(os.path.isfile(
            f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        assert_mapper_logged(self, output_dir, "polish_meta_racon_ill", "bwa mem")
        assert_mapper_logged(self, output_dir, "polish_meta_racon_ill", "sam2paf")

    def test_short_read_polishing_with_bwa_mem2(self):
        # Same pipeline as bwa-mem above (index, then pipe through
        # paftools.js sam2paf), but bwa-mem2 is a distinct binary that CoverM
        # does not pull in itself -- this also proves it is actually present
        # in the polishing pixi env, not just the coverm one.
        output_dir = os.path.join("example", "test_short_read_polishing_bwa_mem2")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-l {data}/pbsim.fq.gz "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--longread-type ccs "
            f"--short-read-mapper bwa-mem2 "
            f"--min-read-size 10 --min-mean-q 1 "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)
        self.assertTrue(os.path.isfile(
            f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        assert_mapper_logged(self, output_dir, "polish_meta_racon_ill", "bwa-mem2 mem")
        assert_mapper_logged(self, output_dir, "polish_meta_racon_ill", "sam2paf")

    def test_short_read_polishing_with_minimap2(self):
        # minimap2/rammap emit PAF directly (no sam2paf conversion needed),
        # but on this non-interleaved paired-read path they do not reliably
        # preserve the /1 or /2 mate suffix in their PAF query names, which
        # has to be reconciled against the still-suffixed original reads for
        # racon's seqkit read-extraction step to find anything at all.
        output_dir = os.path.join("example", "test_short_read_polishing_minimap2")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-l {data}/pbsim.fq.gz "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--longread-type ccs "
            f"--short-read-mapper minimap2 "
            f"--min-read-size 10 --min-mean-q 1 "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)
        self.assertTrue(os.path.isfile(
            f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        assert_mapper_logged(self, output_dir, "polish_meta_racon_ill", "minimap2 -x sr")

    def test_short_read_polishing_with_rammap(self):
        # Same mate-suffix concern as minimap2 above -- rammap is a separate
        # binary sharing minimap2's preset semantics, so this is not covered
        # by the minimap2 test.
        output_dir = os.path.join("example", "test_short_read_polishing_rammap")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-l {data}/pbsim.fq.gz "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--longread-type ccs "
            f"--short-read-mapper rammap "
            f"--min-read-size 10 --min-mean-q 1 "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)
        self.assertTrue(os.path.isfile(
            f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        assert_mapper_logged(self, output_dir, "polish_meta_racon_ill", "rammap -x sr")

    def test_short_read_polishing_with_minibwa(self):
        # minibwa's `map` subcommand requires a prebuilt on-disk index, unlike
        # every other short-read PAF mapper here, which index inline from the
        # FASTA path in one invocation -- polish.py has to build that index
        # itself before calling `minibwa map`.
        output_dir = os.path.join("example", "test_short_read_polishing_minibwa")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-l {data}/pbsim.fq.gz "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--longread-type ccs "
            f"--short-read-mapper minibwa "
            f"--min-read-size 10 --min-mean-q 1 "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)
        self.assertTrue(os.path.isfile(
            f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        assert_mapper_logged(self, output_dir, "polish_meta_racon_ill", "minibwa map")

    def test_long_read_polishing_with_minibwa(self):
        # minibwa's -x preset belongs to its own `map` subcommand, not a
        # top-level flag like minimap2's -- calling it the minimap2 way
        # (`minibwa -x lr ...`) fails with "unknown command '-x'". This
        # exercises the working `minibwa map -f -x lr` invocation end to end
        # via the long-read racon polishing rule (polish_metagenome_flye),
        # which is the code path that previously refused this combination.
        output_dir = os.path.join("example", "test_long_read_polishing_minibwa")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-l {data}/pbsim.fq.gz "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--longread-type ccs "
            f"--long-read-mapper minibwa "
            f"--min-read-size 10 --min-mean-q 1 "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)
        self.assertTrue(os.path.isfile(
            f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        assert_mapper_logged(self, output_dir, "polish_metagenome_flye", "minibwa map -f -x lr")

    def test_long_read_coverage_rammap_default(self):
        # rammap is the new long-read default and nothing else exercises it.
        # --longread-type ont with no --long-read-mapper-model override
        # resolves to the "lr-hq" model, not the legacy "ont"/map-ont preset
        # -- see LONG_READ_TYPE_TO_MODEL in aviary/__init__.py.
        self._run_long_read_coverage(
            "test_long_read_coverage_rammap_default", "", "rammap-lr-hq")

    def test_long_read_coverage_minimap2(self):
        self._run_long_read_coverage(
            "test_long_read_coverage_minimap2",
            "--long-read-mapper minimap2", "minimap2-lr-hq")

    def test_short_read_abundances_use_selected_mapper(self):
        # get_abundances.py runs only when --binning-only is NOT set. Short-read
        # abundances moved from a hardcoded minimap2-sr to the selected mapper,
        # so this pins the new behaviour.
        output_dir = os.path.join("example", "test_short_read_abundances_mapper")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--skip-binners rosella semibin vamb quickbin "
            f"--skip-qc --skip-taxonomy --skip-singlem "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)
        self.assertTrue(os.path.isfile(
            f"{output_dir}/aviary_out/data/coverm_abundances.tsv"))
        assert_mapper_logged(self, output_dir, "coverm_abundances", "-p strobealign")

    def test_long_read_abundances_use_selected_mapper(self):
        # Also exercises the run_coverm guard on a real run: long reads must
        # reach CoverM with a long-read mapper or the script raises.
        output_dir = os.path.join("example", "test_long_read_abundances_mapper")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ont "
            f"--min-read-size 10 --min-mean-q 1 "
            f"--skip-binners rosella semibin vamb quickbin "
            f"--skip-qc --skip-taxonomy --skip-singlem "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)
        self.assertTrue(os.path.isfile(
            f"{output_dir}/aviary_out/data/coverm_abundances.tsv"))
        # --longread-type ont resolves to the "lr-hq" model by default, not
        # the legacy "ont"/map-ont preset -- see LONG_READ_TYPE_TO_MODEL in
        # aviary/__init__.py.
        assert_mapper_logged(self, output_dir, "coverm_abundances", "-p rammap-lr-hq")

    def test_fraction_recovered_uses_selected_mapper(self):
        # read_fraction_recovered is not part of the default DAG -- nothing in
        # a normal `recover` run requests its output. It is reached through
        # `assembly_quality`, which lists it as an input, so it has to be
        # requested explicitly rather than reached by a normal run.
        output_dir = os.path.join("example", "test_fraction_recovered_mapper")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            # Target the output file rather than the assembly_quality rule:
            # that rule also needs www/assembly_stats.txt, which both
            # assembly_size and complete_assembly_with_qc can produce, and
            # snakemake refuses the ambiguity. Requesting the file directly
            # exercises read_fraction_recovered without touching that.
            f"-w www/fraction_recovered/short_fraction_recovered.tsv "
            f"--skip-qc "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)
        self.assertTrue(os.path.isfile(
            f"{output_dir}/aviary_out/www/fraction_recovered/short_fraction_recovered.tsv"))
        assert_mapper_logged(self, output_dir, "fraction_recovered", "-p strobealign")

    def test_long_read_polishing_uses_selected_mapper(self):
        # polish.py generates racon's PAF itself rather than via CoverM, so the
        # command shape differs per aligner. This is a full long-read assembly
        # (no --assembly shortcut) because polishing only runs on one aviary
        # has assembled, and it is the only test that would catch rammap
        # emitting PAF racon cannot consume.
        output_dir = os.path.join("example", "test_long_read_polishing_mapper")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary assemble "
            f"-o {output_dir}/aviary_out "
            f"-l {data}/pbsim.fq.gz "
            f"--longread-type ccs "
            f"--min-read-size 10 --min-mean-q 1 "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)
        self.assertTrue(os.path.isfile(
            f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        # ccs routes through racon (ONT would use medaka instead), so a PAF was
        # generated with the selected long-read mapper.
        assert_mapper_logged(self, output_dir, "polish_metagenome_flye", "rammap")

    def test_short_read_recovery_no_bins(self):
        output_dir = os.path.join("example", "test_short_read_recovery_no_bins")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"--assembly {data}/assembly.fasta "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/tiny_sample.1.fq "
            f"-2 {data}/tiny_sample.2.fq "
            f"--binning-only "
            f"--skip-binners rosella vamb semibin metabat1 "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertFalse(os.path.isfile(bin_info_path))

    def test_short_read_recovery_no_assembly_provided(self):
        output_dir = os.path.join("example", "test_short_read_recovery_no_assembly_provided")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"--binning-only "
            f"--skip-binners rosella vamb metabat "
            f"--skip-qc "
            f"--refinery-max-iterations 0 "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertEqual(num_lines, 3)

        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

    def test_short_read_complete(self):
        output_dir = os.path.join("example", "test_short_read_complete")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary complete "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/wgsim.1.fq.gz "
            f"-2 {data}/wgsim.2.fq.gz "
            f"{request_gpu} "
            f"{singlem_args} "
            f"-n 32 -t 32 "
            f"--strict "
        )
        subprocess.run(cmd, shell=True, check=True)

        bin_info_path = f"{output_dir}/aviary_out/bins/bin_info.tsv"
        self.assertTrue(os.path.isfile(bin_info_path))
        with open(bin_info_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 1)

        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/assembly/final_contigs.fasta"))

        if singlem_metapackage:
            self.assertTrue(os.path.islink(f"{output_dir}/aviary_out/diversity"))
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/metagenome.combined_otu_table.csv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv"))
            self.assertTrue(os.path.getsize(f"{output_dir}/aviary_out/diversity/singlem_appraisal.tsv") > 0)
            self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/diversity/singlem_appraise.svg"))

        gtdbtk_path = f"{output_dir}/aviary_out/taxonomy/gtdbtk.bac120.summary.tsv"
        self.assertTrue(os.path.isfile(gtdbtk_path))
        with open(gtdbtk_path) as f:
            num_lines = sum(1 for _ in f)
        self.assertTrue(num_lines > 1)

        eggnog_paths = glob.glob(os.path.join(output_dir, "aviary_out", "annotation", "eggnog", "*.annotations"))
        self.assertTrue(len(eggnog_paths) > 0)

        for eggnog_path in eggnog_paths:
            self.assertTrue(os.path.isfile(eggnog_path))
            with open(eggnog_path) as f:
                num_lines = sum(1 for _ in f)
            self.assertTrue(num_lines > 1)

    def test_error_integration(self):
        """Expect aviary_assemble to fail with tiny test data, then check logging.
        The small input cannot be assembled, so aviary_assemble will error first.
        We assert Snakemake reports the aviary_assemble error, our handler logs the
        rule failure with a concrete log path, and that an assemble log exists.
        """
        output_dir = os.path.join("example", "test_error_integration")
        setup_output_dir(output_dir)
        cmd = (
            f"aviary recover "
            f"-o {output_dir}/aviary_out "
            f"-1 {data}/tiny_sample_bad.1.fq "
            f"-2 {data}/tiny_sample.2.fq "
            f"--binning-only "
            f"--skip-qc --coassemble "
            f"-n 32 -t 32 "
            f"--strict "
        )

        try:
            proc = subprocess.run(cmd, shell=True, check=True, capture_output=True)
            combined = (proc.stdout or b"").decode("utf-8", errors="replace") + (proc.stderr or b"").decode("utf-8", errors="replace")
        except subprocess.CalledProcessError as e:
            combined = (e.stdout or b"").decode("utf-8", errors="replace") + (e.stderr or b"").decode("utf-8", errors="replace")

        # Snakemake should report a rule error; our wrapper should emit log dumps
        self.assertTrue(("Error in rule assemble_short_reads" in combined) or ("RuleException in rule assemble_short_reads" in combined))
        self.assertIn("Rule failed: assemble_short_reads", combined)
        self.assertTrue(("===== BEGIN LOG (" in combined) and ("===== END LOG (" in combined))

        # Ensure an assemble log exists under the expected logs tree
        assemble_logs = glob.glob(os.path.join(output_dir, "aviary_out", "logs", "assemble_short_reads", "*", "attempt*.log"))
        self.assertTrue(len(assemble_logs) > 0)

        # One of the printed log paths should match an existing file
        printed_paths = re.findall(r"BEGIN LOG \([^)]*\):\s*(.*?)\s*=====", combined)
        if printed_paths:
            # Normalize whitespace and test for existence of at least one printed log
            self.assertTrue(any(os.path.exists(p.strip()) for p in printed_paths))

if __name__ == "__main__":
    unittest.main()
