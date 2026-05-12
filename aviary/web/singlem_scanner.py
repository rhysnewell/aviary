"""
Background SingleM scanner for Aviary web interface.

At server launch, scans all aviary output directories and writes
singlem_profile.tsv (a singlem condense taxonomic profile) for each.
Already-processed directories are skipped (idempotent).

Priority order per output directory:
  1. Pipeline OTU table already exists (data/singlem_out/metagenome.combined_otu_table.csv)
     → run singlem condense directly on it (fast, seconds per sample).
  2. Bin FASTA files exist (recover/bins/ or bins/)
     → run singlem pipe on bins then condense (slower, minutes per sample).
"""
import os
import logging
import shutil
import subprocess
import threading
from pathlib import Path

logger = logging.getLogger(__name__)

PROFILE_FILENAME   = "singlem_profile.tsv"
OTU_TMP_FILENAME   = "singlem_otu_table_tmp.tsv"
PIPELINE_OTU_TABLE = "data/singlem_out/metagenome.combined_otu_table.csv"
BIN_EXTENSIONS     = {".fna", ".fa", ".fasta"}


def _find_bins(output_dir: Path):
    """Return list of bin FASTA files in an aviary output directory."""
    for subdir in (output_dir / "recover" / "bins", output_dir / "bins"):
        if subdir.is_dir():
            bins = [p for p in subdir.iterdir()
                    if p.suffix in BIN_EXTENSIONS and p.is_file()]
            if bins:
                return bins
    return []


def _run_condense(otu_table: Path, profile: Path,
                  metapackage: str, singlem_bin: str) -> bool:
    """
    Run singlem condense on an existing OTU table and write a profile.
    Returns True on success.
    """
    condense_cmd = [
        singlem_bin, "condense",
        "--input-otu-tables", str(otu_table),
        "--metapackage", metapackage,
        "--taxonomic-profile", str(profile),
    ]
    try:
        r = subprocess.run(condense_cmd, capture_output=True, text=True, timeout=int(os.environ.get("SINGLEM_CONDENSE_TIMEOUT", 600)))
        if r.returncode != 0:
            logger.warning(f"[singlem] condense failed ({profile.parent.name}): {r.stderr[-400:]}")
            return False
        return True
    except subprocess.TimeoutExpired:
        logger.warning(f"[singlem] condense timed out: {profile.parent.name}")
        return False


def run_singlem(output_dir: Path, metapackage: str, singlem_bin: str) -> bool:
    """
    Produce singlem_profile.tsv for output_dir. Returns True if successful.

    Fast path: condense the pipeline OTU table that aviary already wrote.
    Slow path: run singlem pipe on bin FASTAs, then condense.
    """
    profile = output_dir / PROFILE_FILENAME
    if profile.exists():
        return True   # already done — skip

    # ── Fast path: pipeline OTU table already exists ──────────────────────
    pipeline_otu = output_dir / PIPELINE_OTU_TABLE
    if pipeline_otu.exists() and pipeline_otu.stat().st_size > 0:
        logger.info(f"[singlem] condense (pipeline OTU table) → {output_dir.name}")
        if _run_condense(pipeline_otu, profile, metapackage, singlem_bin):
            logger.info(f"[singlem] done (fast path): {output_dir.name}")
            return True
        # condense failed — fall through to bin-based path

    # ── Slow path: run pipe on bins then condense ──────────────────────────
    bins = _find_bins(output_dir)
    if not bins:
        logger.debug(f"[singlem] no pipeline OTU table or bins found: {output_dir.name}")
        return False

    otu_tmp = output_dir / OTU_TMP_FILENAME
    try:
        pipe_cmd = [
            singlem_bin, "pipe",
            "--genome-fasta-files", *[str(b) for b in bins],
            "--metapackage", metapackage,
            "--otu-table", str(otu_tmp),
        ]
        logger.info(f"[singlem] pipe (bins) → {output_dir.name}")
        r = subprocess.run(pipe_cmd, capture_output=True, text=True, timeout=int(os.environ.get("SINGLEM_PIPE_TIMEOUT", 1800)))
        if r.returncode != 0:
            logger.warning(f"[singlem] pipe failed ({output_dir.name}): {r.stderr[-400:]}")
            return False

        if not otu_tmp.exists() or otu_tmp.stat().st_size == 0:
            logger.info(f"[singlem] no marker genes found in bins: {output_dir.name}")
            profile.write_text("sample\tcoverage\ttaxonomy\n")
            return True

        logger.info(f"[singlem] condense (bins) → {output_dir.name}")
        if _run_condense(otu_tmp, profile, metapackage, singlem_bin):
            logger.info(f"[singlem] done (bin path): {output_dir.name}")
            return True
        return False

    except subprocess.TimeoutExpired:
        logger.warning(f"[singlem] pipe timed out: {output_dir.name}")
        return False
    except Exception as e:
        logger.error(f"[singlem] error ({output_dir.name}): {e}")
        return False
    finally:
        otu_tmp.unlink(missing_ok=True)


def _scan_worker(output_dirs, metapackage, singlem_bin):
    """Worker: sequentially processes all output dirs at low OS priority."""
    try:
        os.nice(10)
    except (AttributeError, OSError):
        pass

    pending = [Path(d) for d in output_dirs
               if not (Path(d) / PROFILE_FILENAME).exists()]

    if not pending:
        print("[singlem] All directories already have profiles — nothing to do.")
        return

    n = len(pending)
    print(f"[singlem] Scanning {n} director{'y' if n == 1 else 'ies'} in background...")
    done = 0
    for d in pending:
        if run_singlem(d, metapackage, singlem_bin):
            done += 1
    print(f"[singlem] Scan complete — {done}/{n} profiles generated.")


def start_background_scan(output_dirs, metapackage, singlem_bin=None):
    """
    Launch a background daemon thread to run singlem across all output_dirs.
    Returns the thread (or None if prerequisites are missing).
    """
    if not singlem_bin:
        singlem_bin = shutil.which("singlem")
    if not singlem_bin:
        print("[singlem] WARNING: 'singlem' not found in PATH — skipping scan.")
        print("[singlem]   Activate the singlem pixi environment before starting the server.")
        return None
    if not metapackage:
        print("[singlem] WARNING: no metapackage path — skipping scan.")
        print("[singlem]   Set SINGLEM_METAPACKAGE_PATH or pass --metapackage.")
        return None

    t = threading.Thread(
        target=_scan_worker,
        args=(list(output_dirs), metapackage, singlem_bin),
        daemon=True,
        name="singlem-scanner",
    )
    t.start()
    return t
