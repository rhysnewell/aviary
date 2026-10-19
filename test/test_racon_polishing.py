#!/usr/bin/env python
"""End-to-end tests for racon short-read polishing (polish_meta_racon_ill)
with paired reads whose mates differ in length, as they do after fastp.

racon treats reads as single-end and looks them up by ID. Interleaved mates
share an ID, so it paired one mate's overlap with the other's sequence:
"unequal lengths in sequence and overlap file" (isovic/racon#68), or silent
wrong-mate polishing when the mates are the same length. The racon tests in
test_integration.py use untrimmed 150 bp wgsim reads and only ever hit the
silent case.

Marked expensive: they run in the real polishing pixi env. Run them via
`pixi run run-tests-at-cmr` or `pytest --run-expensive`.
"""

import glob
import gzip
import os
import random
import shutil
import subprocess
import unittest

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO_ROOT, "test", "data")
POLISH = os.path.join(REPO_ROOT, "aviary", "modules", "assembly", "scripts", "polish.py")
PIXI_RUN = f"pixi run --manifest-path {os.path.join(REPO_ROOT, 'aviary', 'pixi.toml')}"
THREADS = 16


def _fastq_records(path):
    with gzip.open(path, "rt") as handle:
        while header := handle.readline().rstrip("\n"):
            seq, _, qual = (handle.readline().rstrip("\n") for _ in range(3))
            yield header[1:], seq, qual


def write_paired_reads(out_dir, style, trim=True):
    """wgsim pairs with "suffix" (r/1, r/2) or "illumina" (same ID, CASAVA
    description) headers; mates trimmed to unequal lengths when `trim`.
    Returns (interleaved, r1, r2) paths."""
    rng = random.Random(7)
    paths = [os.path.join(out_dir, f"{style}.{name}.fq.gz") for name in ("interleaved", "R1", "R2")]
    handles = [gzip.open(path, "wt") for path in paths]
    mates = zip(_fastq_records(os.path.join(DATA, "wgsim.1.fq.gz")),
                _fastq_records(os.path.join(DATA, "wgsim.2.fq.gz")))
    for pair in mates:
        for mate, (read_id, seq, qual) in enumerate(pair, start=1):
            if style == "illumina":
                read_id = f"{read_id[:-2].replace('_', ':')} {mate}:N:0:CT+GT"
            cut = rng.randint(0, 25) if trim else 0
            record = f"@{read_id}\n{seq[:len(seq) - cut]}\n+\n{qual[:len(qual) - cut]}\n"
            handles[0].write(record)
            handles[mate].write(record)
    for handle in handles:
        handle.close()
    return paths


def setup_output_dir(output_dir):
    shutil.rmtree(output_dir, ignore_errors=True)
    os.makedirs(output_dir)


@pytest.mark.expensive
class TestRaconPolishingPairedReads(unittest.TestCase):

    def _polish(self, work_dir, reads_1, reads_2, mapper):
        os.makedirs(os.path.join(work_dir, "data"))
        cmd = (
            f"{PIXI_RUN} -e polishing python {POLISH} "
            f"--short-reads-1 {reads_1} --short-reads-2 {reads_2} "
            f"--reference {DATA}/assembly.fasta --output-dir data/polishing "
            f"--output-prefix racon_ill --output-fasta out.fa --rounds 1 "
            f"--long-read-type ont --short-read-mapper {mapper} --illumina True "
            f"--max-cov 200 --threads {THREADS} --coassemble True --log polish.log"
        )
        result = subprocess.run(cmd, shell=True, cwd=work_dir, capture_output=True, text=True)
        log_path = os.path.join(work_dir, "polish.log")
        log = open(log_path).read() if os.path.exists(log_path) else ""
        self.assertEqual(result.returncode, 0,
                         f"{mapper}:\n{result.stderr[-3000:]}\n--- polish.log ---\n{log[-3000:]}")
        self.assertIn("racon return: 0", log)
        with open(os.path.join(work_dir, "out.fa")) as handle:
            self.assertTrue(any(line.startswith(">") for line in handle), "empty assembly")
        with gzip.open(os.path.join(work_dir, "data/polishing/reads.racon_ill.0.fastq.gz"), "rt") as handle:
            ids = [line.split()[0] for i, line in enumerate(handle) if i % 4 == 0]
        self.assertTrue(ids, "no reads were extracted for racon")
        self.assertEqual(len(ids), len(set(ids)), "racon was given two reads with one ID")

    def _check(self, name, style, mappers, trim=True, interleaved=True):
        output_dir = os.path.join(REPO_ROOT, "example", name)
        setup_output_dir(output_dir)
        interleaved_reads, r1, r2 = write_paired_reads(output_dir, style, trim)
        reads = (interleaved_reads, "none") if interleaved else (r1, r2)
        for mapper in mappers:
            with self.subTest(mapper=mapper):
                self._polish(os.path.join(output_dir, mapper), *reads, mapper)

    def test_interleaved_illumina_headers_every_mapper(self):
        # The reported failure: QC'd interleaved reads, Illumina headers.
        self._check("test_racon_interleaved_illumina", "illumina",
                    ["strobealign", "minimap2-sr", "rammap-sr", "minibwa", "bwa-mem"])

    def test_interleaved_suffixed_headers(self):
        self._check("test_racon_interleaved_suffix", "suffix", ["strobealign", "minimap2-sr"])

    def test_interleaved_equal_length_mates(self):
        # Never failed -- polished with the wrong mate. Only the ID check sees it.
        self._check("test_racon_interleaved_equal", "illumina", ["strobealign"], trim=False)

    def test_separate_r1_r2(self):
        self._check("test_racon_separate", "illumina", ["strobealign", "bwa-mem"], interleaved=False)

    def test_aviary_assemble_hybrid_with_trimmed_reads(self):
        # Whole pipeline as users run it: R1/R2 in, qc_short_reads interleaves
        # them, polish_meta_racon_ill polishes the interleaved file.
        output_dir = os.path.join(REPO_ROOT, "example", "test_racon_aviary_assemble")
        setup_output_dir(output_dir)
        _, r1, r2 = write_paired_reads(output_dir, "illumina")
        subprocess.run(
            f"aviary assemble -o {output_dir}/aviary_out -l {DATA}/pbsim.fq.gz "
            f"-1 {r1} -2 {r2} --longread-type ccs --min-read-size 10 --min-mean-q 1 "
            f"-n {THREADS} -t {THREADS}",
            shell=True, check=True)
        self.assertTrue(os.path.isfile(f"{output_dir}/aviary_out/data/final_contigs.fasta"))
        logs = glob.glob(f"{output_dir}/aviary_out/logs/polish_meta_racon_ill/*/attempt*.log")
        self.assertTrue(logs, "polish_meta_racon_ill never ran")
        text = "".join(open(path).read() for path in logs)
        self.assertIn("racon return: 0", text)


if __name__ == '__main__':
    unittest.main()
