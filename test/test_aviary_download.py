#!/usr/bin/env python

#=======================================================================
# Tests for `aviary configure --download`, the database download workflow.
#
# Aviary can download and install the GTDB, EggNOG, SingleM, CheckM2 and
# Metabuli databases via the `--download` flag. Each database is installed
# into the directory pointed to by its associated environment variable, and
# the snakemake `download_databases` workflow touches a
# `<db_dir>.<db>.download.done` marker on success.
#
# This file contains:
#   1. A fast static regression test guarding the fix that lets the download
#      workflow run without reads (download_databases must be in
#      SUBCOMMANDS_WITHOUT_READS).
#   2. An @expensive dry-run test that drives the real CLI/snakemake path and
#      asserts the reads-validation guard is NOT triggered (no download).
#   3. A @download parametrized test that performs the real download of each
#      database and verifies the installed artifacts.
#   4. A @download end-to-end test that downloads ALL databases and then runs a
#      full `aviary complete` walkthrough on the synthetic reads in test/data,
#      proving the downloaded databases work through the whole pipeline.
#
# The real-download tests (3, 4) are marked `download`, NOT `expensive`, so they
# are excluded from the normal expensive suite (e.g. run_tests_at_cmr, which
# collects `-m expensive`) and only run when pytest is given `--run-download`.
# They require:
#   - the aviary pixi environments to be built (`aviary build`); each download
#     runs inside the relevant env (e.g. `pixi run -e checkm2 checkm2 ...`).
#   - network access and substantial disk/time (GTDB in particular is ~100 GB).
# The dry-run test (2) is marked `expensive` (it builds the snakemake DAG but
# downloads nothing) and runs with `--run-expensive`.
#
# RESUMABILITY: databases are downloaded into a persistent directory (DB_BASE,
# set via AVIARY_TEST_DB_DIR -- point it at a roomy disk like /scratch). After a
# database downloads AND passes its works-check, a `.<db>.verified.done`
# checkpoint is written there. If a run is killed by OOM/walltime/etc. (not a
# code problem), simply re-run: already-verified databases are skipped and only
# the interrupted/remaining ones are attempted. Delete a checkpoint to force a
# re-download of that database.
#
# This can be run standalone with `python test_aviary_download.py` or via pixi:
#                   pixi run run-test-aviary-download
#
# Distributed under the GNU General Public License v3 or later.
#=======================================================================

import os
import glob
import tempfile
import subprocess

import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PROCESSOR_PY = os.path.join(REPO_ROOT, "aviary", "modules", "processor.py")
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

# Threads/cores for the full walkthrough; override with AVIARY_TEST_CORES.
WALKTHROUGH_CORES = os.environ.get("AVIARY_TEST_CORES", "16")

# Where databases are downloaded. GTDB alone is ~100 GB so use a roomy disk.
# Override with AVIARY_TEST_DB_DIR; defaults to ~/aviary_test_dbs locally.
# On HPC: export AVIARY_TEST_DB_DIR=/scratch/herholdt/aviary_test_dbs
DB_BASE = os.environ.get(
    "AVIARY_TEST_DB_DIR",
    os.path.join(os.path.expanduser("~"), "aviary_test_dbs"),
)

# Temporary directory for aviary intermediate files. Must be set or aviary
# will prompt interactively (which fails in pytest).
# On HPC: export AVIARY_TEST_TMPDIR=/scratch/herholdt
_tmpdir = os.environ.get("AVIARY_TEST_TMPDIR", os.path.expanduser("~"))
os.environ.setdefault("TMPDIR", _tmpdir)

# (db key passed to --download, env var that points to its install dir,
#  list of glob patterns -- relative to the db dir -- that must match a
#  non-empty file/dir after a successful download)
DATABASE_CASES = [
    pytest.param("checkm2", "CHECKM2DB", ["*.dmnd"], id="checkm2"),
    pytest.param("eggnog", "EGGNOG_DATA_DIR", ["eggnog.db"], id="eggnog"),
    pytest.param("singlem", "SINGLEM_METAPACKAGE_PATH", [""], id="singlem"),
    pytest.param("metabuli", "METABULI_DB_PATH", ["*"], id="metabuli"),
    pytest.param("gtdb", "GTDBTK_DATA_PATH", ["*"], id="gtdb"),
]

# db key -> env var that points to its install directory
DB_ENV_VARS = {db: env for db, env, _ in (p.values for p in DATABASE_CASES)}


def _done_marker(db_dir, db):
    # Matches annotation.smk: get_db_done_file -> f"{get_db_dir(db)}.{db}.download.done"
    return f"{db_dir}.{db}.download.done"


def _verified_marker(db):
    # Our own checkpoint: written only after download AND the works-check pass.
    return os.path.join(DB_BASE, f".{db}.verified.done")


def _run_aviary(args, cwd=REPO_ROOT):
    return subprocess.run(["aviary", *args], cwd=cwd, check=True)


def _run_build(cwd=REPO_ROOT):
    return subprocess.run(["pixi", "run", "aviary", "build"], cwd=cwd, check=True)


def _verify_artifacts(db, db_dir, artifact_globs):
    """The 'simple call to see if it works' check: the install dir holds real,
    non-empty content matching the expected artifact patterns. Raises on failure
    so no checkpoint is written and a restart will retry this database."""
    entries = os.listdir(db_dir)
    assert entries, f"{db} install dir is empty after download: {db_dir}"
    for pattern in artifact_globs:
        if pattern == "":
            assert any(
                os.path.getsize(os.path.join(db_dir, e)) > 0
                if os.path.isfile(os.path.join(db_dir, e))
                else os.listdir(os.path.join(db_dir, e))
                for e in entries
            ), f"{db} produced only empty artifacts in {db_dir}"
            continue
        matches = glob.glob(os.path.join(db_dir, pattern))
        assert matches, f"{db} expected artifact matching '{pattern}' in {db_dir}"
        assert any(
            (os.path.isfile(m) and os.path.getsize(m) > 0) or
            (os.path.isdir(m) and os.listdir(m))
            for m in matches
        ), f"{db} artifact '{pattern}' exists but is empty in {db_dir}"


def _ensure_db(db, env_var, artifact_globs, cores, monkeypatch):
    """Resumable download + verify for one database.

    Always points `env_var` at the persistent install dir (so downstream steps
    find it). If a verified checkpoint already exists, returns immediately
    without re-downloading. Otherwise downloads (only when the snakemake
    download marker is absent -- re-running the download into a populated GTDB
    dir would error), runs the works-check, then writes the checkpoint.

    A crash (OOM/walltime) before the checkpoint is written leaves the (already
    finished) databases checkpointed and only the interrupted one to retry.

    Returns True if it had to download/verify, False if resumed from checkpoint.
    """
    db_dir = os.path.join(DB_BASE, db)
    os.makedirs(db_dir, exist_ok=True)
    # aviary `configure` resolves every db path on startup and prompts
    # interactively for any unset one (EOFError under pytest). Point all
    # non-target dbs at placeholder dirs so only the db under test matters.
    for _db, _env in DB_ENV_VARS.items():
        monkeypatch.setenv(_env, os.path.join(DB_BASE, _db))
    monkeypatch.setenv(env_var, db_dir)

    if os.path.isfile(_verified_marker(db)):
        return False  # already downloaded and verified in a previous run

    if not os.path.isfile(_done_marker(db_dir, db)):
        _run_aviary([
            "configure", "--download", db,
            "-o", os.path.join(DB_BASE, "_configure_out"),
            "-n", "1", "-t", str(cores),
        ])
        assert os.path.isfile(_done_marker(db_dir, db)), \
            f"download .done marker missing for {db}: {_done_marker(db_dir, db)}"

    _verify_artifacts(db, db_dir, artifact_globs)
    open(_verified_marker(db), "w").close()  # checkpoint: safe to resume past this
    return True


# --------------------------------------------------------------------------- #
# 1. Fast static regression test (no download, always runs)
# --------------------------------------------------------------------------- #
def test_download_databases_is_reads_exempt():
    """download_databases must be exempt from the reads check, otherwise
    `aviary configure --download` fails with
    'both long_reads and short_reads_1 are set to none'."""
    with open(PROCESSOR_PY) as f:
        contents = f.read()

    # Find the SUBCOMMANDS_WITHOUT_READS assignment and confirm membership.
    line = next(
        (l for l in contents.splitlines() if l.strip().startswith("SUBCOMMANDS_WITHOUT_READS")),
        None,
    )
    assert line is not None, "SUBCOMMANDS_WITHOUT_READS not found in processor.py"
    assert "download_databases" in line, (
        "download_databases must be in SUBCOMMANDS_WITHOUT_READS so database-only "
        "downloads skip the reads-validation guard"
    )


# --------------------------------------------------------------------------- #
# 2. Dry-run: exercises the real CLI + snakemake DAG without downloading.
#    Before the fix this raised the reads error during workflow construction.
# --------------------------------------------------------------------------- #
@pytest.mark.expensive
def test_configure_download_dryrun_no_reads_error(tmp_path, monkeypatch):
    db_dir = tmp_path / "checkm2"
    db_dir.mkdir()
    monkeypatch.setenv("CHECKM2DB", str(db_dir))

    proc = subprocess.run(
        [
            "aviary", "configure",
            "--download", "checkm2",
            "-o", str(tmp_path / "aviary_out"),
            "-n", "1", "-t", "1",
            "--dryrun",
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    combined = proc.stdout + proc.stderr
    assert "both long_reads and short_reads_1 are set to" not in combined, (
        f"reads-validation guard was triggered during a database-only download:\n{combined}"
    )
    assert proc.returncode == 0, f"dry-run failed:\n{combined}"


@pytest.fixture(scope="session", autouse=False)
def aviary_build():
    """Ensure all pixi sub-environments are built before download tests run."""
    _run_build()


# --------------------------------------------------------------------------- #
# 3. Real download + verification, one database per parametrization.
#    Gated by `download` marker / --run-download (NOT --run-expensive), so the
#    expensive suite never triggers a multi-GB download.
#    Run a single light one with: pytest --run-download -k checkm2
# --------------------------------------------------------------------------- #
@pytest.mark.download
@pytest.mark.parametrize("db, env_var, artifact_globs", DATABASE_CASES)
def test_configure_download_installs_database(db, env_var, artifact_globs, monkeypatch, aviary_build):
    if os.path.isfile(_verified_marker(db)):
        pytest.skip(
            f"{db} already downloaded and verified (checkpoint {_verified_marker(db)}); "
            f"delete it to force a re-download"
        )
    _ensure_db(db, env_var, artifact_globs, cores=8, monkeypatch=monkeypatch)
    assert os.path.isfile(_verified_marker(db))


# --------------------------------------------------------------------------- #
# 4. End-to-end: download ALL databases, then run a full `aviary complete`
#    walkthrough on the synthetic short reads (test/data/wgsim.{1,2}.fq.gz) and
#    verify the assemble -> recover -> annotate pipeline produced real outputs.
#    This proves the freshly-downloaded databases are actually usable.
#
#    Heavy: downloads every database (GTDB alone is ~100 GB) and runs the full
#    pipeline (GTDB-Tk needs ~128 GB+ RAM). Gated behind --run-download; intended
#    to be submitted to a large node. Tune cores with AVIARY_TEST_CORES.
# --------------------------------------------------------------------------- #
@pytest.mark.download
def test_download_then_complete_walkthrough(tmp_path, monkeypatch, aviary_build):
    # 1) download + verify every database into the persistent DB_BASE. Each is
    #    checkpointed, so a restart after OOM/walltime skips finished databases.
    for db, env_var, artifact_globs in (p.values for p in DATABASE_CASES):
        _ensure_db(db, env_var, artifact_globs, cores=WALKTHROUGH_CORES, monkeypatch=monkeypatch)

    # 2) run the full assemble -> recover -> annotate pipeline on synthetic reads
    out = tmp_path / "aviary_out"
    _run_aviary([
        "complete",
        "-o", str(out),
        "-1", os.path.join(DATA_DIR, "wgsim.1.fq.gz"),
        "-2", os.path.join(DATA_DIR, "wgsim.2.fq.gz"),
        "-n", WALKTHROUGH_CORES, "-t", WALKTHROUGH_CORES,
    ])

    # 3) verify each stage produced real output (mirrors test_short_read_complete)
    # assembly
    assert os.path.isfile(out / "data" / "final_contigs.fasta"), "assembly missing"
    assert os.path.getsize(out / "data" / "final_contigs.fasta") > 0

    # recovery / binning
    bin_info = out / "bins" / "bin_info.tsv"
    assert os.path.isfile(bin_info), "binning bin_info.tsv missing"
    with open(bin_info) as f:
        assert sum(1 for _ in f) > 1, "no bins recovered"

    # annotation: GTDB-Tk taxonomy (needs the downloaded GTDB db)
    gtdbtk = out / "taxonomy" / "gtdbtk.bac120.summary.tsv"
    assert os.path.isfile(gtdbtk), "GTDB-Tk summary missing (gtdb db not used?)"
    with open(gtdbtk) as f:
        assert sum(1 for _ in f) > 1, "GTDB-Tk produced no classifications"

    # annotation: EggNOG (needs the downloaded eggnog db)
    eggnog = glob.glob(str(out / "annotation" / "eggnog" / "*.annotations"))
    assert eggnog, "no EggNOG annotation files produced"
    for path in eggnog:
        assert os.path.getsize(path) > 0, f"empty EggNOG annotation: {path}"


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
