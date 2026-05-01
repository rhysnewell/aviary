#!/usr/bin/env python3
"""
Aviary Web Interface
Live log viewer for Snakemake workflows.

Usage:
    python aviary/web/server.py --output-dir /path/to/aviary/output [--port 8090]
    python -m aviary.web.server --output-dir /path/to/aviary/output
"""

import os
import re
import json
import time
import math
import bisect
import socket
import random as _random
import threading
import argparse
import logging
from pathlib import Path
from datetime import datetime
from flask import Flask, request, Response, jsonify

logger = logging.getLogger(__name__)

# singlem_scanner.py writes this file after running singlem condense.
# This is the correct input for the phylo tree — NOT the raw OTU table.
SINGLEM_PROFILE_FILENAME  = "singlem_profile.tsv"
GTDBTK_NEWICK_SUBPATH    = Path("data") / "gtdbtk" / "classify" / "gtdbtk.bac120.classify.tree.6.tree"
GTDBTK_SUMMARY_SUBPATH   = Path("data") / "gtdbtk" / "classify" / "gtdbtk.bac120.summary.tsv"

# GFA stats cache: (path_str, mtime, version) -> result dict
# Bump _GFA_CACHE_VERSION whenever _parse_gfa output schema changes
_GFA_CACHE_VERSION = 2
_gfa_cache      = {}
_gfa_cache_lock = threading.Lock()

# ---------------------------------------------------------------------------
# Log parsing
# ---------------------------------------------------------------------------

TIMESTAMP_RE  = re.compile(r'^\[(\w+ \w+ +\d+ \d+:\d+:\d+ \d+)\]$')
RULE_RE       = re.compile(r'^(local)?rule (\w+):$')
CHECKPOINT_RE = re.compile(r'^checkpoint (\w+):$')
JOBID_RE      = re.compile(r'^\s+jobid: (\d+)$')
SUBMITTED_RE  = re.compile(r"^Submitted job (\d+) with external jobid '([^']+)'")
FINISHED_RE   = re.compile(r'^Finished jobid: (\d+) \(Rule: ([\w.\-]+)\)$')
PROGRESS_RE   = re.compile(r'^(\d+) of (\d+) steps \((\d+)%\) done$')
HOST_RE       = re.compile(r'^host: (.+)$')
JOB_STATS_RE  = re.compile(r'^([A-Za-z_]\S*)\s+(\d+)$')
THREADS_RE    = re.compile(r'^\s+threads: (\d+)$')
RESOURCES_RE  = re.compile(r'^\s+resources: (.+)$')
BENCHMARK_RE  = re.compile(r'^\s+benchmark: (.+)$')
LOG_PATH_RE   = re.compile(r'log_path=([^\s,]+)')
MEM_MB_RE     = re.compile(r'mem_mb=(\d+)')
RUNTIME_RE    = re.compile(r'runtime=(\d+)')
INPUT_RE      = re.compile(r'^\s+input: (.+)$')
OUTPUT_RE     = re.compile(r'^\s+output: (.+)$')
ERROR_RE      = re.compile(r'^Error in rule (\w+):')
WFLOW_ERROR_RE = re.compile(r'^WorkflowError:')


def _parse_ts(ts_str):
    if not ts_str:
        return None
    try:
        return datetime.strptime(ts_str.strip(), "%a %b %d %H:%M:%S %Y")
    except ValueError:
        return None


def _fmt_duration(start_str, end_str):
    s = _parse_ts(start_str)
    e = _parse_ts(end_str)
    if s is None or e is None:
        return None
    total = int((e - s).total_seconds())
    if total < 60:
        return f"{total}s"
    if total < 3600:
        return f"{total // 60}m {total % 60}s"
    h = total // 3600
    m = (total % 3600) // 60
    return f"{h}h {m}m"


def parse_snakemake_log(log_path):
    result = {
        "host": None,
        "log_file": str(log_path),
        "log_modified": None,
        "total_steps": 0,
        "completed_steps": 0,
        "progress_pct": 0,
        "job_stats": {},
        "jobs": {},
        "rule_to_jobid": {},
        "errors": [],
        "is_complete": False,
        "workflow_failed": False,
        "started_at": None,
    }

    p = Path(log_path)
    if not p.exists():
        return result

    result["log_modified"] = p.stat().st_mtime
    try:
        content = p.read_text(errors="replace")
    except OSError:
        return result

    lines = content.splitlines()
    current_ts = None
    in_job_stats = False
    i = 0

    while i < len(lines):
        line = lines[i]

        m = HOST_RE.match(line)
        if m:
            result["host"] = m.group(1).strip()
            i += 1
            continue

        if line.strip() == "Job stats:":
            in_job_stats = True
            i += 1
            continue

        if in_job_stats:
            if line.startswith("---") or line.strip() in ("job", "count", ""):
                if line.strip() == "":
                    in_job_stats = False
                i += 1
                continue
            m = JOB_STATS_RE.match(line.strip())
            if m:
                name, count = m.group(1), int(m.group(2))
                if name == "total":
                    result["total_steps"] = count
                else:
                    result["job_stats"][name] = count
                i += 1
                continue
            in_job_stats = False

        m = TIMESTAMP_RE.match(line)
        if m:
            current_ts = m.group(1)
            if result["started_at"] is None:
                result["started_at"] = current_ts
            i += 1
            continue

        # Rule or checkpoint
        rm = RULE_RE.match(line)
        cm = CHECKPOINT_RE.match(line)
        if rm or cm:
            if cm:
                is_local = False
                rule_name = cm.group(1)
            else:
                is_local = bool(rm.group(1))
                rule_name = rm.group(2)

            job = {
                "rule": rule_name,
                "is_local": is_local,
                "start_time": current_ts,
                "start_ts": _parse_ts(current_ts).isoformat() if _parse_ts(current_ts) else None,
                "end_time": None,
                "end_ts": None,
                "duration": None,
                "status": "running",
                "jobid": None,
                "external_jobid": None,
                "threads": 1 if is_local else None,
                "mem_mb": None,
                "runtime_limit": None,
                "log_path": None,
                "benchmark": None,
                "input": None,
                "output": None,
            }

            i += 1
            while i < len(lines) and (lines[i].startswith("    ") or lines[i].startswith("\t")):
                dl = lines[i]
                jm = JOBID_RE.match(dl)
                if jm:
                    job["jobid"] = int(jm.group(1))
                tm = THREADS_RE.match(dl)
                if tm:
                    job["threads"] = int(tm.group(1))
                resm = RESOURCES_RE.match(dl)
                if resm:
                    res = resm.group(1)
                    lm = LOG_PATH_RE.search(res)
                    if lm:
                        job["log_path"] = lm.group(1)
                    mm = MEM_MB_RE.search(res)
                    if mm:
                        job["mem_mb"] = int(mm.group(1))
                    rtm = RUNTIME_RE.search(res)
                    if rtm:
                        job["runtime_limit"] = int(rtm.group(1))
                bm = BENCHMARK_RE.match(dl)
                if bm:
                    job["benchmark"] = bm.group(1).strip()
                im = INPUT_RE.match(dl)
                if im:
                    job["input"] = im.group(1).strip()
                om = OUTPUT_RE.match(dl)
                if om:
                    job["output"] = om.group(1).strip()
                i += 1

            jid = job["jobid"]
            if jid is not None:
                result["jobs"][jid] = job
                result["rule_to_jobid"][rule_name] = jid
            else:
                synth = -(len(result["jobs"]) + 1)
                job["jobid"] = synth
                result["jobs"][synth] = job
                result["rule_to_jobid"][rule_name] = synth
            continue

        m = SUBMITTED_RE.match(line)
        if m:
            jid = int(m.group(1))
            if jid in result["jobs"]:
                result["jobs"][jid]["external_jobid"] = m.group(2)
            i += 1
            continue

        m = FINISHED_RE.match(line)
        if m:
            jid = int(m.group(1))
            rule = m.group(2)
            target = result["jobs"].get(jid) or result["jobs"].get(result["rule_to_jobid"].get(rule))
            if target:
                target["status"] = "completed"
                target["end_time"] = current_ts
                target["end_ts"] = _parse_ts(current_ts).isoformat() if _parse_ts(current_ts) else None
                target["duration"] = _fmt_duration(target["start_time"], current_ts)
            i += 1
            continue

        m = PROGRESS_RE.match(line)
        if m:
            result["completed_steps"] = int(m.group(1))
            result["total_steps"] = int(m.group(2))
            result["progress_pct"] = int(m.group(3))
            i += 1
            continue

        m = ERROR_RE.match(line)
        if m:
            rule = m.group(1)
            result["errors"].append({"rule": rule, "timestamp": current_ts})
            jid = result["rule_to_jobid"].get(rule)
            if jid is not None and jid in result["jobs"]:
                result["jobs"][jid]["status"] = "failed"
            result["workflow_failed"] = True
            i += 1
            continue

        if WFLOW_ERROR_RE.match(line):
            result["workflow_failed"] = True

        i += 1

    if result["total_steps"] > 0 and result["completed_steps"] >= result["total_steps"]:
        result["is_complete"] = True
        for job in result["jobs"].values():
            if job["status"] == "running":
                job["status"] = "completed"

    jobs_list = sorted(result["jobs"].values(), key=lambda x: x.get("start_ts") or "")
    result["jobs_list"] = jobs_list

    statuses = [j["status"] for j in result["jobs"].values()]
    result["counts"] = {
        "running":   statuses.count("running"),
        "completed": statuses.count("completed"),
        "failed":    statuses.count("failed"),
        "total":     len(statuses),
    }

    return result


def _quality_tier(completeness, contamination):
    try:
        c, x = float(completeness), float(contamination)
        if c >= 90 and x <= 5:  return "high"
        if c >= 50 and x <= 10: return "medium"
    except (ValueError, TypeError):
        pass
    return "low"


def _parse_assembly_stats(stats_path):
    """Parse www/assembly_stats.txt into a dict."""
    p = Path(stats_path)
    if not p.exists():
        return None
    try:
        text = p.read_text()
    except OSError:
        return None
    result = {}
    patterns = {
        "total_length_mb": r"Main genome scaffold sequence total:\s+([\d.]+) Mb",
        "contig_count":    r"Main genome scaffold total:\s+([\d,]+)",
        "n50":             r"Main genome scaffold N/L50:\s+[\d,]+/([\d,]+)",
        "l50":             r"Main genome scaffold N/L50:\s+([\d,]+)/",
        "max_length":      r"Max scaffold length:\s+(\S+)",
        "scaffolds_50kb":  r"Number of scaffolds > 50 KB:\s+([\d,]+)",
    }
    for key, pat in patterns.items():
        m = re.search(pat, text, re.MULTILINE)
        if m:
            val = m.group(1).replace(",", "")
            try:
                result[key] = float(val) if "." in val else int(val)
            except ValueError:
                result[key] = val
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith("A\tC\tG"):
            try:
                vals = lines[i + 1].split("\t")
                result["gc_content"] = round(float(vals[7]) * 100, 1)
            except (IndexError, ValueError):
                pass
            break
    return result or None


def _parse_bin_info(bin_info_path):
    """Parse bins/bin_info.tsv into list of dicts with numeric coercion."""
    p = Path(bin_info_path)
    if not p.exists():
        return None
    try:
        lines = p.read_text().splitlines()
    except OSError:
        return None
    if len(lines) < 2:
        return []
    headers = lines[0].split("\t")
    numeric = {"Completeness", "Contamination", "Coding_Density", "Contig_N50",
               "Genome_Size", "GC_Content", "Total_Coding_Sequences", "Total_Contigs",
               "Max_Contig_Length", "Average_Gene_Length", "Circular contigs",
               "Circular bp", "Circular fraction"}
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        row = {}
        for i, h in enumerate(headers):
            val = parts[i] if i < len(parts) else ""
            if h in numeric:
                try:
                    row[h] = float(val)
                except (ValueError, TypeError):
                    row[h] = None
            else:
                row[h] = val
        rows.append(row)
    return rows


def find_latest_log(output_dir):
    log_dir = Path(output_dir) / ".snakemake" / "log"
    if not log_dir.exists():
        return None
    logs = sorted(log_dir.glob("*.snakemake.log"), key=lambda p: p.stat().st_mtime, reverse=True)
    return logs[0] if logs else None


def list_logs(output_dir):
    log_dir = Path(output_dir) / ".snakemake" / "log"
    if not log_dir.exists():
        return []
    return sorted(log_dir.glob("*.snakemake.log"), key=lambda p: p.stat().st_mtime, reverse=True)


def read_job_log(output_dir, log_path_rel, max_bytes=200_000):
    p = Path(output_dir) / log_path_rel
    if not p.exists():
        # Job may have been retried — search for any attempt log in the rule dir
        rule_dir = p.parent.parent
        candidates = sorted(rule_dir.glob("*/attempt*.log"), key=lambda x: x.stat().st_mtime, reverse=True)
        if candidates:
            p = candidates[0]
        else:
            return None, f"Log file not found: {p}"
    try:
        size = p.stat().st_size
        with open(p, errors="replace") as f:
            if size > max_bytes:
                f.seek(max(0, size - max_bytes))
                return f"[... showing last {max_bytes // 1024}KB of {size // 1024}KB ...]\n\n" + f.read(), None
            return f.read(), None
    except OSError as e:
        return None, str(e)


def read_benchmark(output_dir, benchmark_path):
    p = Path(output_dir) / benchmark_path
    if not p.exists():
        return None
    try:
        lines = p.read_text().splitlines()
        if len(lines) < 2:
            return None
        return dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    except OSError:
        return None


# ---------------------------------------------------------------------------
# HTTP server (no external dependencies)
# ---------------------------------------------------------------------------

def _parse_gfa(path_str):
    """
    Stream-parse a GFA file returning histogram stats + reservoir scatter sample.
    Uses fixed log-scale bins so no per-segment storage is needed — safe for
    very large graphs (8 M+ segments).
    """
    LEN_BINS   = 60          # log10-spaced length bins
    DEPTH_BINS = 60          # linear depth bins 0-200x
    DEPTH_CAP  = 200.0
    SCATTER_N  = 3000        # reservoir sample size

    log_min, log_max = 2.0, 7.0   # 100 bp → 10 Mbp

    len_counts   = [0] * LEN_BINS
    depth_counts = [0] * DEPTH_BINS
    depth_over   = 0

    res_len   = []   # reservoir for scatter
    res_dep   = []
    res_k     = 0    # total segments seen (for reservoir)

    total_seg   = 0
    total_len   = 0
    total_dep   = 0.0
    max_len     = 0
    min_len     = 10**9
    total_links = 0

    # Coverage threshold counters
    n_cov_1x  = 0;  bp_cov_1x  = 0
    n_cov_5x  = 0;  bp_cov_5x  = 0
    n_cov_10x = 0;  bp_cov_10x = 0
    n_low_cov = 0   # contigs with depth < 1×

    # Binning length threshold counters
    n_1kbp = 0;  bp_1kbp = 0
    n_2kbp = 0;  bp_2kbp = 0

    # Dead-end tracking: nodes seen in L-lines; disabled for very large graphs
    _DEAD_LIMIT   = 4_000_000
    _linked_nodes = set()
    _dead_track   = True

    # For approximate N50 via histogram (length × count cumsum from large end)
    # We track a finer length histogram in parallel
    N50_BINS   = 200
    n50_counts    = [0]   * N50_BINS   # count of segments per bin
    n50_totals    = [0]   * N50_BINS   # total bp per bin
    n50_weighted  = [0.0] * N50_BINS   # depth-weighted bp per bin (for cov-weighted N50)

    with open(path_str, 'r', buffering=1 << 20) as fh:
        for line in fh:
            c = line[0] if line else ''

            if c == 'S':
                # Fast field extraction: find tab positions without full split
                try:
                    t1 = line.index('\t')
                    t2 = line.index('\t', t1 + 1)
                    t3 = line.find('\t', t2 + 1)
                    # Sequence length from field width (avoids allocating the string)
                    seq_len = (t3 - t2 - 1) if t3 != -1 else (len(line.rstrip('\n')) - t2 - 1)
                    # Handle '*' placeholder
                    if line[t2 + 1] == '*':
                        seq_len = 0

                    # Parse coverage from tags or megahit-style segment name
                    depth = 0.0
                    name  = line[t1 + 1:t2]
                    if t3 != -1:
                        rest = line[t3 + 1:]
                        dp_i = rest.find('DP:f:')
                        if dp_i != -1:
                            dp_s = dp_i + 5
                            dp_e = rest.find('\t', dp_s)
                            try:
                                depth = float(rest[dp_s:dp_e] if dp_e != -1 else rest[dp_s:].rstrip())
                            except ValueError:
                                pass
                        elif 'KC:i:' in rest:
                            # fallback: k-mer count tag
                            kc_i = rest.find('KC:i:') + 5
                            kc_e = rest.find('\t', kc_i)
                            try:
                                depth = float(rest[kc_i:kc_e] if kc_e != -1 else rest[kc_i:].rstrip())
                            except ValueError:
                                pass
                        if 'LN:i:' in rest:
                            ln_i = rest.find('LN:i:') + 5
                            ln_e = rest.find('\t', ln_i)
                            try:
                                seq_len = int(rest[ln_i:ln_e] if ln_e != -1 else rest[ln_i:].rstrip())
                            except ValueError:
                                pass
                    # Megahit encodes coverage in the node name: NODE_N_length_L_cov_C_ID_I
                    if depth == 0.0 and '_cov_' in name:
                        cov_i = name.find('_cov_') + 5
                        cov_e = name.find('_', cov_i)
                        try:
                            depth = float(name[cov_i:cov_e] if cov_e != -1 else name[cov_i:])
                        except ValueError:
                            pass
                except ValueError:
                    continue

                total_seg += 1
                total_len += seq_len
                total_dep += depth
                if seq_len > max_len:
                    max_len = seq_len
                if seq_len < min_len:
                    min_len = seq_len

                # Coverage thresholds
                if depth < 1.0:
                    n_low_cov += 1
                if depth >= 1.0:
                    n_cov_1x  += 1;  bp_cov_1x  += seq_len
                if depth >= 5.0:
                    n_cov_5x  += 1;  bp_cov_5x  += seq_len
                if depth >= 10.0:
                    n_cov_10x += 1;  bp_cov_10x += seq_len

                # Binning length thresholds
                if seq_len >= 1000:
                    n_1kbp += 1;  bp_1kbp += seq_len
                if seq_len >= 2000:
                    n_2kbp += 1;  bp_2kbp += seq_len

                # Length histogram (log10)
                if seq_len > 0:
                    log_l = math.log10(seq_len)
                    idx = int((log_l - log_min) / (log_max - log_min) * LEN_BINS)
                    idx = max(0, min(idx, LEN_BINS - 1))
                    len_counts[idx] += 1

                    # N50 histogram (finer bins)
                    n50_idx = int((log_l - log_min) / (log_max - log_min) * N50_BINS)
                    n50_idx = max(0, min(n50_idx, N50_BINS - 1))
                    n50_counts[n50_idx]   += 1
                    n50_totals[n50_idx]   += seq_len
                    n50_weighted[n50_idx] += seq_len * depth

                # Depth histogram (linear 0-200)
                if depth >= DEPTH_CAP:
                    depth_over += 1
                else:
                    d_idx = int(depth / DEPTH_CAP * DEPTH_BINS)
                    depth_counts[max(0, min(d_idx, DEPTH_BINS - 1))] += 1

                # Reservoir sampling for scatter
                res_k += 1
                if len(res_len) < SCATTER_N:
                    res_len.append(seq_len)
                    res_dep.append(round(depth, 2))
                else:
                    j = _random.randint(0, res_k - 1)
                    if j < SCATTER_N:
                        res_len[j] = seq_len
                        res_dep[j] = round(depth, 2)

            elif c == 'L':
                total_links += 1
                if _dead_track:
                    try:
                        t1 = line.index('\t')
                        t2 = line.index('\t', t1 + 1)
                        t3 = line.index('\t', t2 + 1)
                        t4 = line.index('\t', t3 + 1)
                        _linked_nodes.add(line[t1 + 1:t2])   # from_name
                        _linked_nodes.add(line[t3 + 1:t4])   # to_name
                        if len(_linked_nodes) > _DEAD_LIMIT:
                            _dead_track = False
                            _linked_nodes.clear()
                    except ValueError:
                        pass

    if total_seg == 0:
        return {"error": "No segments found in GFA"}

    if min_len == 10**9:
        min_len = 0

    # Approximate N50 from the N50 histogram (iterate large→small bins)
    half_total = total_len / 2
    cumsum = 0
    n50_approx = 0
    for i in range(N50_BINS - 1, -1, -1):
        cumsum += n50_totals[i]
        if cumsum >= half_total:
            frac_lo = log_min + i / N50_BINS * (log_max - log_min)
            frac_hi = log_min + (i + 1) / N50_BINS * (log_max - log_min)
            n50_approx = int(10 ** ((frac_lo + frac_hi) / 2))
            break

    # Coverage-weighted N50: weight each contig by length × depth
    total_weighted = sum(n50_weighted)
    half_weighted  = total_weighted / 2
    cumw = 0
    n50_cov_weighted = 0
    for i in range(N50_BINS - 1, -1, -1):
        cumw += n50_weighted[i]
        if cumw >= half_weighted:
            frac_lo = log_min + i / N50_BINS * (log_max - log_min)
            frac_hi = log_min + (i + 1) / N50_BINS * (log_max - log_min)
            n50_cov_weighted = int(10 ** ((frac_lo + frac_hi) / 2))
            break

    # Build axis labels for length histogram
    len_edges = [
        round(10 ** (log_min + i / LEN_BINS * (log_max - log_min)))
        for i in range(LEN_BINS + 1)
    ]
    # Build axis labels for depth histogram
    depth_edges = [round(DEPTH_CAP * i / DEPTH_BINS, 1) for i in range(DEPTH_BINS + 1)]

    dead_ends = (total_seg - len(_linked_nodes)) if _dead_track else None

    return {
        "stats": {
            "total_segments": total_seg,
            "total_length":   total_len,
            "n50":            n50_approx,
            "n50_cov_weighted": n50_cov_weighted,
            "max_length":     max_len,
            "min_length":     min_len,
            "total_links":    total_links,
            "avg_links_per_seg": round(2 * total_links / total_seg, 2) if total_seg else 0,
            "avg_depth":      round(total_dep / total_seg, 2) if total_seg else 0,
            # Coverage thresholds
            "n_cov_1x":  n_cov_1x,   "bp_cov_1x":  bp_cov_1x,
            "n_cov_5x":  n_cov_5x,   "bp_cov_5x":  bp_cov_5x,
            "n_cov_10x": n_cov_10x,  "bp_cov_10x": bp_cov_10x,
            "n_low_cov": n_low_cov,
            # Binning length thresholds
            "n_1kbp":  n_1kbp,  "bp_1kbp":  bp_1kbp,
            "n_2kbp":  n_2kbp,  "bp_2kbp":  bp_2kbp,
            # Graph topology
            "dead_ends": dead_ends,
        },
        "length_hist":  {"edges": len_edges,   "counts": len_counts},
        "depth_hist":   {"edges": depth_edges, "counts": depth_counts, "overflow": depth_over, "cap": DEPTH_CAP},
        "scatter":      [[res_len[i], res_dep[i]] for i in range(len(res_len))],
    }


LOGO_PATH = Path(__file__).parent / "Aviary_logo_nobg.png"
LANDING_TEMPLATE_PATH = Path(__file__).parent / "templates" / "landing.html"
TEMPLATE_PATH = Path(__file__).parent / "templates" / "index.html"
DOCS_TEMPLATE_PATH = Path(__file__).parent / "templates" / "docs.html"
EXPORT_TEMPLATE_PATH = Path(__file__).parent / "templates" / "export.html"
GRAPH_TEMPLATE_PATH  = Path(__file__).parent / "templates" / "graph.html"
GUIDE_TEMPLATE_PATH    = Path(__file__).parent / "templates" / "web_guide.html"
ASSEMBLY_TEMPLATE_PATH = Path(__file__).parent / "templates" / "assembly_graph.html"
PHYLO_TEMPLATE_PATH    = Path(__file__).parent / "templates" / "phylo.html"
DEFAULT_OUTPUT_DIR = None
PROJECT_ROOT = None


def _build_taxonomy_tree(root_dir):
    """
    Build a GTDB-Tk taxonomy tree from all bin_info.tsv files under root_dir.
    Returns a D3-compatible nested tree plus a flat bin list.
    """
    RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

    root_p = Path(root_dir)
    dirs = find_all_output_dirs_raw(root_dir)

    # nodes: path_key -> {name, full, rank, children_keys set, bin_count, hq, mq, lq}
    nodes = {}
    ROOT_KEY = "__root__"
    nodes[ROOT_KEY] = {"name": "Life", "full": "Life", "rank": "root",
                       "children": set(), "count": 0, "hq": 0, "mq": 0, "lq": 0}
    bins_list = []

    for d in dirs:
        sample_dir = Path(d)
        try:
            rel = sample_dir.relative_to(root_p)
        except ValueError:
            rel = sample_dir
        parts = rel.parts
        if len(parts) >= 2:
            assembler, sample = parts[-2], parts[-1]
        elif len(parts) == 1:
            assembler, sample = "default", parts[0]
        else:
            continue

        # Determine batch date the same way api_summary does
        logs = list_logs(d)
        if logs:
            log_name = logs[0].name.replace(".snakemake.log", "")
            try:
                dt = datetime.strptime(log_name[:15], "%Y-%m-%dT%H%M%S")
                batch_key = dt.strftime("%Y-%m-%d")
            except ValueError:
                batch_key = log_name[:10]
        else:
            try:
                batch_key = datetime.fromtimestamp(sample_dir.stat().st_mtime).strftime("%Y-%m-%d")
            except OSError:
                batch_key = "Unknown"

        bins_dir = sample_dir / "bins"
        bin_rows = _parse_bin_info(bins_dir / "bin_info.tsv")
        if not bin_rows:
            bins_dir = sample_dir / "recover" / "bins"
            bin_rows = _parse_bin_info(bins_dir / "bin_info.tsv")
        if not bin_rows:
            continue

        for b in bin_rows:
            classification = (b.get("classification") or "").strip()
            if not classification or classification in ("N/A", "none", "None"):
                continue

            comp  = b.get("Completeness")
            cont  = b.get("Contamination")
            qual  = _quality_tier(comp, cont)

            taxa = [t.strip() for t in classification.split(";")]

            parent_key = ROOT_KEY
            for i, taxon in enumerate(taxa):
                if not taxon or taxon in ("d__", "p__", "c__", "o__", "f__", "g__", "s__"):
                    break
                # Build a unique path key by joining taxa up to this level
                node_key = ";".join(taxa[:i + 1])
                if node_key not in nodes:
                    rank = RANKS[i] if i < len(RANKS) else "unknown"
                    # Strip rank prefix (d__, p__, etc.)
                    display = taxon[3:] if len(taxon) > 3 and taxon[1:3] == "__" else taxon
                    nodes[node_key] = {
                        "name":     display,
                        "full":     taxon,
                        "rank":     rank,
                        "children": set(),
                        "count":    0,
                        "hq": 0, "mq": 0, "lq": 0,
                    }
                nodes[node_key]["count"] += 1
                if qual == "high":   nodes[node_key]["hq"] += 1
                elif qual == "medium": nodes[node_key]["mq"] += 1
                else:               nodes[node_key]["lq"] += 1
                nodes[parent_key]["children"].add(node_key)
                parent_key = node_key

            bins_list.append({
                "name":           b.get("Name", ""),
                "batch":          batch_key,
                "sample":         sample,
                "assembler":      assembler,
                "classification": classification,
                "completeness":   comp,
                "contamination":  cont,
                "quality":        qual,
                "genome_size":    b.get("Genome_Size"),
                "gc_content":     b.get("GC_Content"),
            })

    def _to_dict(key):
        n = nodes[key]
        result = {
            "name":  n["name"],
            "full":  n["full"],
            "rank":  n["rank"],
            "count": n["count"],
            "hq":    n["hq"],
            "mq":    n["mq"],
            "lq":    n["lq"],
        }
        children = [_to_dict(ck) for ck in sorted(n["children"])]
        if children:
            result["children"] = children
        return result

    tree = _to_dict(ROOT_KEY)
    return {"tree": tree, "bins": bins_list, "total": len(bins_list)}


def find_all_output_dirs_raw(root):
    """Internal alias used by _build_taxonomy_tree before find_all_output_dirs is defined."""
    root = Path(root)
    results = []
    try:
        for log_dir in root.rglob(".snakemake/log"):
            if log_dir.is_dir():
                results.append(str(log_dir.parent.parent))
    except PermissionError:
        pass
    return sorted(results)


def find_all_output_dirs(root):
    """Recursively find all directories containing .snakemake/log/ under root."""
    root = Path(root)
    results = []
    try:
        for log_dir in root.rglob(".snakemake/log"):
            if log_dir.is_dir():
                results.append(str(log_dir.parent.parent))
    except PermissionError:
        pass
    return sorted(results)


# ---------------------------------------------------------------------------
# SingleM OTU table parsing (reads aviary's existing pipeline output)
# ---------------------------------------------------------------------------

def _parse_condensed_profile(path):
    """
    Parse a singlem condensed profile (output of `singlem condense`).
    Columns: sample, coverage, taxonomy
    Each row is one taxon — already deduplicated and abundance-normalised by singlem.
    Returns list of {taxonomy, coverage} dicts ready to build the tree.
    """
    rows = []
    try:
        with open(path) as f:
            lines = [l.rstrip("\n") for l in f if l.strip()]
        if len(lines) < 2:
            return []
        header = lines[0].split("\t")
        cov_idx = next((i for i, h in enumerate(header) if h == "coverage"), 1)
        tax_idx = next((i for i, h in enumerate(header) if h == "taxonomy"), 2)
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) <= max(cov_idx, tax_idx):
                continue
            try:
                cov = float(parts[cov_idx])
                tax = parts[tax_idx].strip()
                if tax and cov > 0:
                    rows.append({"taxonomy": tax, "coverage": cov})
            except ValueError:
                continue
    except OSError:
        pass
    return rows


def _build_singlem_tree(root_dir):
    """
    Build a taxonomy tree from singlem condensed profiles (singlem_profile.tsv).
    Written by singlem_scanner.py via `singlem condense` — one row per taxon,
    already deduplicated and abundance-normalised. Do NOT use the raw OTU table.
    Coverage values are normalised to percentages across all samples.
    """
    RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    root_p = Path(root_dir)
    dirs   = find_all_output_dirs_raw(root_dir)

    total_dirs    = len(dirs)
    profiled_dirs = sum(1 for d in dirs if (Path(d) / SINGLEM_PROFILE_FILENAME).exists())

    nodes = {}
    ROOT_KEY = "__root__"
    nodes[ROOT_KEY] = {
        "name": "Life", "full": "Life", "rank": "root",
        "children": set(), "coverage": 0.0, "count": 0,
    }

    profile_rows = []   # flat list returned to the client for per-sample filtering

    for d in dirs:
        profile_path = Path(d) / SINGLEM_PROFILE_FILENAME
        if not profile_path.exists():
            continue

        # Derive sample / assembler / batch from directory path (same logic as elsewhere)
        sample_dir = Path(d)
        try:
            rel = sample_dir.relative_to(root_p)
        except ValueError:
            rel = sample_dir
        parts = rel.parts
        if len(parts) >= 2:
            assembler, sample = parts[-2], parts[-1]
        elif len(parts) == 1:
            assembler, sample = "default", parts[0]
        else:
            continue

        logs = list_logs(d)
        if logs:
            log_name = logs[0].name.replace(".snakemake.log", "")
            try:
                dt = datetime.strptime(log_name[:15], "%Y-%m-%dT%H%M%S")
                batch_key = dt.strftime("%Y-%m-%d")
            except ValueError:
                batch_key = log_name[:10]
        else:
            try:
                batch_key = datetime.fromtimestamp(
                    sample_dir.stat().st_mtime).strftime("%Y-%m-%d")
            except OSError:
                batch_key = "Unknown"

        rows = _parse_condensed_profile(str(profile_path))
        if not rows:
            continue

        # Normalise this sample's coverage to relative abundance (sum → 1.0)
        # so that multi-sample aggregation gives mean RA, not depth-weighted sums.
        sample_total = sum(r["coverage"] for r in rows) or 1.0
        norm_rows = [{"taxonomy": r["taxonomy"],
                      "coverage": r["coverage"] / sample_total}
                     for r in rows]

        for row in norm_rows:
            tax = row["taxonomy"].strip()
            cov = row["coverage"]

            # Strip "Root; " prefix if present
            if tax.lower().startswith("root;"):
                tax = tax[tax.index(";") + 1:].strip().lstrip("; ")

            raw_taxa = [t.strip() for t in tax.split("; ") if t.strip()]
            # Truncate at the first empty rank level (e.g. "c__", "o__")
            # so unclassified levels don't create phantom nodes in the tree.
            taxa = []
            for t in raw_taxa:
                if len(t) == 3 and t[1:3] == "__":
                    break   # e.g. "c__" — unclassified from here
                taxa.append(t)
            if not taxa:
                continue

            nodes[ROOT_KEY]["coverage"] += cov

            parent_key = ROOT_KEY
            for i, taxon in enumerate(taxa):
                node_key = "; ".join(taxa[: i + 1])
                if node_key not in nodes:
                    rank    = RANKS[i] if i < len(RANKS) else "unknown"
                    display = taxon[3:] if len(taxon) > 3 and taxon[1:3] == "__" else taxon
                    nodes[node_key] = {
                        "name":     display,
                        "full":     taxon,
                        "rank":     rank,
                        "children": set(),
                        "coverage": 0.0,
                        "samples":  set(),
                    }
                nodes[node_key]["coverage"] += cov
                nodes[node_key]["samples"].add(sample)
                nodes[parent_key]["children"].add(node_key)
                parent_key = node_key

            profile_rows.append({
                "taxonomy":  tax,
                "coverage":  cov,   # already normalised to RA fraction
                "sample":    sample,
                "assembler": assembler,
                "batch":     batch_key,
            })

    if not profile_rows:
        return None

    n_samples = profiled_dirs or 1
    # coverage_pct = mean relative abundance across all samples that have profiles.
    # Each sample was normalised to sum-to-1 before accumulation, so dividing by
    # n_samples gives the mean RA fraction; multiply by 100 for percentage.
    for n in nodes.values():
        n["coverage_pct"] = round(n["coverage"] / n_samples * 100, 2)
        n["sample_count"] = len(n.get("samples", set()))

    def _to_dict(key):
        n      = nodes[key]
        result = {
            "name":         n["name"],
            "full":         n["full"],
            "rank":         n["rank"],
            "_key":         key if key != ROOT_KEY else "",
            "sample_count": n["sample_count"],
            "coverage":     round(n["coverage"], 6),
            "coverage_pct": n["coverage_pct"],
        }
        children = [_to_dict(ck) for ck in sorted(n["children"])]
        if children:
            result["children"] = children
        return result

    return {
        "data_source":  "singlem",
        "tree":         _to_dict(ROOT_KEY),
        "profile_rows": profile_rows,
        "total":        len(profile_rows),
        "scan_status":  {"total": total_dirs, "done": profiled_dirs},
    }


# ---------------------------------------------------------------------------
# GTDB-Tk phylogenetic tree (Newick) — parse, prune, annotate
# ---------------------------------------------------------------------------

def _is_reference_genome(name):
    """Return True if the leaf name is a GTDB reference genome accession."""
    return not name or name.startswith("GB_GCA_") or name.startswith("RS_GCF_")


def _parse_newick(text):
    """
    Parse a Newick string into nested dicts.
    Node format: {name, length, support, children:[...]}
    Handles GTDB-Tk format: quoted labels like '100.0:o__CG2-30-54-11',
    float bootstrap values on internal nodes, and branch lengths.
    """
    text = text.strip().rstrip(";").strip()
    pos = [0]

    def _node():
        children = []
        if pos[0] < len(text) and text[pos[0]] == "(":
            pos[0] += 1
            while True:
                children.append(_node())
                if pos[0] >= len(text) or text[pos[0]] in (")", ";"):
                    if pos[0] < len(text) and text[pos[0]] == ")":
                        pos[0] += 1
                    break
                if text[pos[0]] == ",":
                    pos[0] += 1

        # Parse label (possibly single-quoted)
        name = ""
        if pos[0] < len(text) and text[pos[0]] == "'":
            pos[0] += 1
            while pos[0] < len(text) and text[pos[0]] != "'":
                name += text[pos[0]]
                pos[0] += 1
            if pos[0] < len(text):
                pos[0] += 1
        else:
            while pos[0] < len(text) and text[pos[0]] not in (",", ")", ":", ";"):
                name += text[pos[0]]
                pos[0] += 1

        # Parse branch length
        length = 0.0
        if pos[0] < len(text) and text[pos[0]] == ":":
            pos[0] += 1
            ls = ""
            while pos[0] < len(text) and text[pos[0]] not in (",", ")", ";"):
                ls += text[pos[0]]
                pos[0] += 1
            try:
                length = float(ls)
            except ValueError:
                pass

        # Decode bootstrap from internal node label
        support = None
        node_name = name
        if children and name:
            if ":" in name:
                # e.g. '100.0:o__CG2-30-54-11'
                parts = name.split(":", 1)
                try:
                    support = float(parts[0])
                    node_name = parts[1]
                except ValueError:
                    pass
            else:
                try:
                    support = float(name)
                    node_name = ""
                except ValueError:
                    pass

        return {"name": node_name, "length": length, "support": support, "children": children}

    return _node()


def _annotate_has_user(node):
    """Mark each node with has_user=True if it has any user MAG leaf descendants."""
    if not node["children"]:
        node["is_user"]  = not _is_reference_genome(node["name"])
        node["has_user"] = node["is_user"]
        return
    for c in node["children"]:
        _annotate_has_user(c)
    node["is_user"]  = False
    node["has_user"] = any(c["has_user"] for c in node["children"])


def _count_leaves(node):
    if not node["children"]:
        return 1
    return sum(_count_leaves(c) for c in node["children"])


def _prune_to_user_mags(node, max_ref_context=3):
    """Keep user MAGs and a small number of reference neighbours for context."""
    if not node.get("has_user") and not node.get("is_user"):
        return None
    if not node["children"]:
        return node if node.get("is_user") else None

    user_kids    = [c for c in node["children"] if c.get("has_user") or c.get("is_user")]
    ref_leaves   = [c for c in node["children"]
                    if not c.get("has_user") and not c.get("is_user") and not c["children"]]
    ref_subtrees = [c for c in node["children"]
                    if not c.get("has_user") and not c.get("is_user") and c["children"]]

    pruned = [_prune_to_user_mags(c, max_ref_context) for c in user_kids]
    pruned = [c for c in pruned if c is not None]

    context = ref_leaves[:max_ref_context]
    for rs in ref_subtrees[:1]:          # one collapsed reference subtree for context
        n = _count_leaves(rs)
        context.append({
            "name": f"[{n} ref. genomes]",
            "length": rs["length"],
            "support": None,
            "children": [],
            "is_user": False,
            "has_user": False,
            "is_collapsed": True,
            "collapsed_count": n,
        })

    new_node = dict(node)
    new_node["children"] = pruned + context
    return new_node


def _collapse_single_children(node):
    """Collapse internal nodes that have exactly one child (sum branch lengths)."""
    if not node["children"]:
        return node
    node = dict(node)
    node["children"] = [_collapse_single_children(c) for c in node["children"]]
    if len(node["children"]) == 1:
        child = dict(node["children"][0])
        child["length"] = round(child["length"] + node["length"], 6)
        if not child.get("name") and node.get("name"):
            child["name"] = node["name"]
        if child["support"] is None and node.get("support") is not None:
            child["support"] = node["support"]
        return child
    return node


def _build_phylo_newick(output_dir):
    """
    Parse GTDB-Tk Newick tree for one aviary output directory, prune to user
    MAGs + reference context, and annotate with taxonomy and SingleM abundance.
    Returns JSON-serialisable dict, or None if no tree exists.
    """
    out          = Path(output_dir)
    tree_path    = out / GTDBTK_NEWICK_SUBPATH
    summary_path = out / GTDBTK_SUMMARY_SUBPATH

    if not tree_path.exists():
        return None

    # ── Taxonomy from GTDB-Tk summary TSV ──────────────────────────────
    taxonomy_map = {}
    if summary_path.exists():
        try:
            with open(summary_path) as f:
                hdr   = f.readline().rstrip("\n").split("\t")
                ug_i  = next((i for i, h in enumerate(hdr) if h == "user_genome"),    0)
                cl_i  = next((i for i, h in enumerate(hdr) if h == "classification"), 1)
                for line in f:
                    p = line.rstrip("\n").split("\t")
                    if len(p) > max(ug_i, cl_i):
                        taxonomy_map[p[ug_i]] = p[cl_i]
        except (OSError, ValueError):
            pass

    # ── Abundance from singlem condensed profile ────────────────────────
    singlem_by_deepest = {}   # deepest taxonomy token → coverage_pct
    profile_path = out / SINGLEM_PROFILE_FILENAME
    if profile_path.exists():
        rows = _parse_condensed_profile(str(profile_path))
        if rows:
            total = sum(r["coverage"] for r in rows) or 1.0
            for r in rows:
                deepest = r["taxonomy"].split("; ")[-1].strip()
                singlem_by_deepest[deepest] = round(
                    singlem_by_deepest.get(deepest, 0.0) + r["coverage"] / total * 100, 3)

    # ── Parse Newick ────────────────────────────────────────────────────
    try:
        tree = _parse_newick(tree_path.read_text())
    except Exception as e:
        logger.error(f"[phylo_newick] Parse failed: {e}")
        return None

    _annotate_has_user(tree)
    if not tree.get("has_user"):
        return {"error": "No user MAGs found in tree", "tree": None}

    # ── Prune + collapse ────────────────────────────────────────────────
    pruned = _prune_to_user_mags(tree, max_ref_context=3)
    if not pruned:
        return None
    pruned = _collapse_single_children(pruned)

    # ── Annotate leaves ──────────────────────────────────────────────────
    def _annotate(node):
        if not node["children"]:
            name             = node["name"]
            node["taxonomy"] = taxonomy_map.get(name, "")
            node["is_user"]  = not _is_reference_genome(name) and not node.get("is_collapsed")
            if node["is_user"] and node["taxonomy"]:
                deepest = node["taxonomy"].split(";")[-1].strip()
                node["coverage_pct"] = singlem_by_deepest.get(deepest, 0.0)
            else:
                node["coverage_pct"] = 0.0
        else:
            for c in node["children"]:
                _annotate(c)
    _annotate(pruned)

    # ── Serialise ────────────────────────────────────────────────────────
    def _to_json(node):
        r = {
            "name":         node["name"],
            "length":       node["length"],
            "support":      node.get("support"),
            "is_user":      node.get("is_user", False),
            "taxonomy":     node.get("taxonomy", ""),
            "coverage_pct": node.get("coverage_pct", 0.0),
        }
        if node.get("is_collapsed"):
            r["is_collapsed"]    = True
            r["collapsed_count"] = node.get("collapsed_count", 0)
        if node["children"]:
            r["children"] = [_to_json(c) for c in node["children"]]
        return r

    user_mags = []
    def _collect(node):
        if not node["children"] and node.get("is_user"):
            user_mags.append({"name": node["name"], "taxonomy": node.get("taxonomy", ""),
                               "coverage_pct": node.get("coverage_pct", 0.0)})
        for c in node.get("children", []):
            _collect(c)
    _collect(pruned)

    return {
        "tree":            _to_json(pruned),
        "user_mags":       user_mags,
        "total_user_mags": len(user_mags),
        "has_abundance":   bool(singlem_by_deepest),
    }


# ---------------------------------------------------------------------------
# Flask application
# ---------------------------------------------------------------------------

app = Flask(__name__, template_folder='templates')

def _p(key, default=None):
    """Get a query param from the current Flask request."""
    return request.args.get(key, default)

def _json_r(data, status=200):
    return Response(json.dumps(data), status=status, mimetype='application/json')

def _html_r(html_str, status=200):
    return Response(html_str, status=status, mimetype='text/html; charset=utf-8')


# ---------------------------------------------------------------------------
# Page routes
# ---------------------------------------------------------------------------

@app.route('/')
def landing():
    try:
        return _html_r(LANDING_TEMPLATE_PATH.read_text())
    except OSError as e:
        return _html_r(f"<pre>Template error: {e}</pre>", 500)

@app.route('/dashboard')
def index():
    pinned_dir = _p('output_dir') or DEFAULT_OUTPUT_DIR or ""
    try:
        tmpl = TEMPLATE_PATH.read_text()
        tmpl = tmpl.replace(
            "{{ ('\"' + output_dir + '\"') if output_dir else 'null' }}",
            f'"{pinned_dir}"' if pinned_dir else "null"
        )
        return _html_r(tmpl)
    except OSError as e:
        return _html_r(f"<pre>Template error: {e}</pre>", 500)

@app.route('/docs')
def docs():
    try:
        return _html_r(DOCS_TEMPLATE_PATH.read_text())
    except OSError as e:
        return _html_r(f"<pre>Template error: {e}</pre>", 500)

@app.route('/graph')
def graph():
    try:
        return _html_r(GRAPH_TEMPLATE_PATH.read_text())
    except OSError as e:
        return _html_r(f"<pre>Template error: {e}</pre>", 500)

@app.route('/assembly')
def assembly():
    try:
        return _html_r(ASSEMBLY_TEMPLATE_PATH.read_text())
    except OSError as e:
        return _html_r(f"<pre>Template error: {e}</pre>", 500)

@app.route('/guide')
def guide():
    try:
        return _html_r(GUIDE_TEMPLATE_PATH.read_text())
    except OSError as e:
        return _html_r(f"<pre>Template error: {e}</pre>", 500)

@app.route('/export')
def export():
    try:
        return _html_r(EXPORT_TEMPLATE_PATH.read_text())
    except OSError as e:
        return _html_r(f"<pre>Template error: {e}</pre>", 500)

@app.route('/view')
def phylo():
    try:
        return _html_r(PHYLO_TEMPLATE_PATH.read_text())
    except OSError as e:
        return _html_r(f"<pre>Template error: {e}</pre>", 500)

@app.route('/logo')
def logo():
    if LOGO_PATH.exists():
        return Response(LOGO_PATH.read_bytes(), mimetype='image/png')
    return Response('Logo not found', status=404, mimetype='text/plain')


# ---------------------------------------------------------------------------
# API routes
# ---------------------------------------------------------------------------

@app.route('/api/structure')
def api_structure():
    root = _p('root') or DEFAULT_OUTPUT_DIR or PROJECT_ROOT
    if not root:
        return _json_r({"samples": {}, "root": None})
    dirs = find_all_output_dirs(root)
    root_p = Path(root)

    run_groups = {}
    ungrouped  = []
    for d in dirs:
        try:
            rel = Path(d).relative_to(root_p)
        except ValueError:
            continue
        parts = rel.parts
        if len(parts) >= 3:
            run_groups.setdefault(parts[0], []).append((d, parts))
        else:
            ungrouped.append((d, parts))

    def _run_latest_log_date(entries):
        best_mtime, best_date = 0, None
        for d, _ in entries:
            for lf in list_logs(d):
                mt = lf.stat().st_mtime
                if mt > best_mtime:
                    best_mtime = mt
                    best_date  = lf.name[:10]
        return best_date or "unknown"

    run_dates   = {rid: _run_latest_log_date(entries) for rid, entries in run_groups.items()}
    most_recent = max(run_groups, key=lambda r: run_dates[r]) if run_groups else None

    all_entries = list(ungrouped)
    for run_id, entries in run_groups.items():
        suffix = "" if run_id == most_recent else f" ({run_dates[run_id]})"
        all_entries.extend((d, parts, suffix) for d, parts in entries)

    COMMIT_HASH_RE = re.compile(r'^[0-9a-f]{40}(_\S+)?$')
    samples = {}
    for item in all_entries:
        suffix = item[2] if len(item) == 3 else ""
        d, parts = item[0], item[1]
        if len(parts) == 2:
            assembler, sample = parts[0], parts[1]
            # Skip entries where the "assembler" is actually a commit hash (old 2-level layout)
            if COMMIT_HASH_RE.match(assembler):
                continue
        elif len(parts) == 1:
            assembler, sample = "default", parts[0]
        elif len(parts) >= 3:
            assembler, sample = parts[-2], parts[-1]
        else:
            continue
        sample_key = f"{sample}{suffix}"
        all_logs = list_logs(d)
        log_path = all_logs[0] if all_logs else None
        if log_path:
            data = parse_snakemake_log(log_path)
            data["log_name"] = Path(log_path).name
            data.pop("rule_to_jobid", None)
            data.pop("jobs", None)
        else:
            data = {
                "host": None, "log_file": None, "log_modified": None,
                "total_steps": 0, "completed_steps": 0, "progress_pct": 0,
                "job_stats": {}, "jobs_list": [], "errors": [],
                "is_complete": False, "workflow_failed": False,
                "started_at": None, "counts": {},
                "log_name": None,
            }
        data["available_logs"] = [p.name for p in all_logs]
        data["output_dir"] = d
        data["is_latest_run"] = (suffix == "")
        if sample_key not in samples:
            samples[sample_key] = {}
        samples[sample_key][assembler] = data
    return _json_r({"samples": samples, "root": root})


@app.route('/api/output_dirs')
def api_output_dirs():
    root = _p('root') or DEFAULT_OUTPUT_DIR or PROJECT_ROOT
    if not root:
        return _json_r([])
    dirs = find_all_output_dirs(root)
    root_p = Path(root)
    display = []
    for d in dirs:
        try:
            display.append({"path": d, "label": str(Path(d).relative_to(root_p))})
        except ValueError:
            display.append({"path": d, "label": d})
    return _json_r(display)


@app.route('/api/logs')
def api_logs():
    output_dir = _p('output_dir') or DEFAULT_OUTPUT_DIR
    if not output_dir:
        return _json_r([])
    return _json_r([p.name for p in list_logs(output_dir)])


@app.route('/api/status')
def api_status():
    output_dir = _p('output_dir') or DEFAULT_OUTPUT_DIR
    log_name   = _p('log')
    if not output_dir:
        return _json_r({"error": "No output_dir specified"}, 400)
    if log_name:
        log_path = Path(output_dir) / ".snakemake" / "log" / log_name
    else:
        log_path = find_latest_log(output_dir)
    if not log_path or not Path(log_path).exists():
        return _json_r({"error": f"No snakemake log found in {output_dir}/.snakemake/log/"}, 404)
    data = parse_snakemake_log(log_path)
    data["log_name"] = Path(log_path).name
    data["available_logs"] = [p.name for p in list_logs(output_dir)]
    data.pop("rule_to_jobid", None)
    data.pop("jobs", None)
    return _json_r(data)


@app.route('/api/job_log')
def api_job_log():
    output_dir = _p('output_dir') or DEFAULT_OUTPUT_DIR
    log_path   = _p('log_path')
    if not output_dir or not log_path:
        return _json_r({"error": "Missing parameters"}, 400)
    content, err = read_job_log(output_dir, log_path)
    if err:
        return _json_r({"error": err}, 404)
    available_attempts = []
    try:
        p_full = Path(output_dir) / log_path
        attempt_dir = p_full.parent
        if attempt_dir.is_dir():
            available_attempts = sorted(
                str(sib.relative_to(Path(output_dir)))
                for sib in attempt_dir.glob("attempt*.log")
                if sib.is_file()
            )
    except Exception:
        pass
    return _json_r({"content": content, "log_path": log_path, "available_attempts": available_attempts})


@app.route('/api/benchmark')
def api_benchmark():
    output_dir = _p('output_dir') or DEFAULT_OUTPUT_DIR
    bench_path = _p('path')
    if not output_dir or not bench_path:
        return _json_r({"error": "Missing parameters"}, 400)
    data = read_benchmark(output_dir, bench_path)
    if data is None:
        return _json_r({"error": "Benchmark not found or not yet written"}, 404)
    return _json_r(data)


@app.route('/api/summary')
def api_summary():
    root = _p('root') or _p('output_dir') or DEFAULT_OUTPUT_DIR
    if not root:
        return _json_r({"error": "No output_dir specified"}, 400)
    root_p = Path(root)
    output_dirs = find_all_output_dirs(root)
    batches = {}
    for d in output_dirs:
        sample_dir = Path(d)
        try:
            rel = sample_dir.relative_to(root_p)
        except ValueError:
            rel = sample_dir
        parts = rel.parts
        if len(parts) == 2:
            assembler, sample = parts[0], parts[1]
        elif len(parts) == 1:
            assembler, sample = "default", parts[0]
        else:
            assembler, sample = (parts[-2] if len(parts) >= 2 else "default"), parts[-1]
        logs = list_logs(d)
        if logs:
            log_name = logs[0].name.replace(".snakemake.log", "")
            try:
                dt = datetime.strptime(log_name[:15], "%Y-%m-%dT%H%M%S")
                date_key = dt.strftime("%Y-%m-%d")
            except ValueError:
                date_key = log_name[:10]
        else:
            try:
                date_key = datetime.fromtimestamp(sample_dir.stat().st_mtime).strftime("%Y-%m-%d")
            except OSError:
                date_key = "Unknown"
        bins_dir = sample_dir / "recover" / "bins"
        if not bins_dir.exists():
            bins_dir = sample_dir / "bins"
        bins = _parse_bin_info(bins_dir / "bin_info.tsv")
        final_bins_dir = bins_dir / "final_bins"
        if final_bins_dir.exists():
            mag_count = len(list(final_bins_dir.glob("*.fna")) + list(final_bins_dir.glob("*.fasta")))
        else:
            mag_count = None
        hq = mq = lq = 0
        if bins:
            for b in bins:
                t = _quality_tier(b.get("Completeness"), b.get("Contamination"))
                if t == "high":     hq += 1
                elif t == "medium": mq += 1
                else:               lq += 1
        assembly = _parse_assembly_stats(sample_dir / "www" / "assembly_stats.txt")
        entry = {
            "sample":    sample,
            "assembler": assembler,
            "mag_count": mag_count,
            "hq": hq, "mq": mq, "lq": lq,
            "assembly":  assembly,
            "bins":      bins or [],
            "has_circular": bool(bins and any("Circular contigs" in b for b in bins)),
            "has_taxonomy": bool(bins and any(b.get("classification") for b in bins)),
        }
        batches.setdefault(date_key, []).append(entry)
    return _json_r({"batches": dict(sorted(batches.items(), reverse=True))})


@app.route('/api/gfa_stats')
def api_gfa_stats():
    output_dir = _p('output_dir')
    if not output_dir:
        return _json_r({"error": "output_dir required"}, 400)
    gfa_path = Path(output_dir) / "assembly" / "assembly_graph.gfa"
    if not gfa_path.exists():
        return _json_r({"error": "assembly_graph.gfa not found"}, 404)
    try:
        mtime     = gfa_path.stat().st_mtime
        cache_key = (str(gfa_path), mtime, _GFA_CACHE_VERSION)
        with _gfa_cache_lock:
            if cache_key in _gfa_cache:
                return _json_r(_gfa_cache[cache_key])
        result = _parse_gfa(str(gfa_path))
        with _gfa_cache_lock:
            _gfa_cache[cache_key] = result
        return _json_r(result)
    except Exception as e:
        return _json_r({"error": str(e)}, 500)


@app.route('/api/gfa_available')
def api_gfa_available():
    root = _p('root') or DEFAULT_OUTPUT_DIR or PROJECT_ROOT
    if not root:
        return _json_r({})
    dirs = find_all_output_dirs(root)
    result = {d: (Path(d) / "assembly" / "assembly_graph.gfa").exists() for d in dirs}
    return _json_r(result)


@app.route('/api/taxonomy_tree')
def api_taxonomy_tree():
    root = _p('root') or DEFAULT_OUTPUT_DIR or PROJECT_ROOT
    if not root:
        return _json_r({"error": "No root specified"}, 400)
    try:
        # Prefer singlem profiles when available — better taxonomy + real abundance
        singlem = _build_singlem_tree(root)
        if singlem and singlem.get("total", 0) > 0:
            return _json_r(singlem)
        # Fall back to GTDB-Tk classification strings
        data = _build_taxonomy_tree(root)
        data["data_source"]  = "gtdbtk"
        data["scan_status"]  = {
            "total": len(find_all_output_dirs(root)),
            "done":  0,
        }
        return _json_r(data)
    except Exception as e:
        return _json_r({"error": str(e)}, 500)


@app.route('/api/singlem_status')
def api_singlem_status():
    root = _p('root') or DEFAULT_OUTPUT_DIR or PROJECT_ROOT
    if not root:
        return _json_r({"error": "No root specified"}, 400)
    dirs  = find_all_output_dirs(root)
    total = len(dirs)
    done  = sum(1 for d in dirs if (Path(d) / SINGLEM_PROFILE_FILENAME).exists())
    return _json_r({"total": total, "done": done, "scanning": False})


@app.route('/api/phylo_newick')
def api_phylo_newick():
    root       = _p('root')       or DEFAULT_OUTPUT_DIR or PROJECT_ROOT
    output_dir = _p('output_dir') or DEFAULT_OUTPUT_DIR
    if not output_dir and root:
        for d in find_all_output_dirs(root):
            if (Path(d) / GTDBTK_NEWICK_SUBPATH).exists():
                output_dir = d
                break
    if not output_dir:
        return _json_r({"error": "No output directory found"}, 400)
    try:
        result = _build_phylo_newick(output_dir)
        if not result:
            return _json_r({"error": "No GTDB-Tk Newick tree found in output directory"}, 404)
        return _json_r(result)
    except Exception as e:
        return _json_r({"error": str(e)}, 500)


# ---------------------------------------------------------------------------
# GFA cache warm-up
# ---------------------------------------------------------------------------

def _warmup_gfa_cache(root):
    """Background thread: pre-parse all GFA files into cache at low OS priority."""
    try:
        os.nice(10)  # Lower priority so user requests always take precedence
    except (AttributeError, OSError):
        pass  # os.nice not available on all platforms

    dirs = find_all_output_dirs(root)
    gfa_paths = [
        Path(d) / "assembly" / "assembly_graph.gfa"
        for d in dirs
        if (Path(d) / "assembly" / "assembly_graph.gfa").exists()
    ]

    if not gfa_paths:
        return

    print(f"[warmup] Pre-caching {len(gfa_paths)} GFA file(s) in background...")
    warmed = 0
    for gfa_path in gfa_paths:
        try:
            mtime     = gfa_path.stat().st_mtime
            cache_key = (str(gfa_path), mtime, _GFA_CACHE_VERSION)
            with _gfa_cache_lock:
                already_cached = cache_key in _gfa_cache
            if not already_cached:
                result = _parse_gfa(str(gfa_path))
                with _gfa_cache_lock:
                    _gfa_cache[cache_key] = result
            warmed += 1
        except Exception:
            pass  # Don't let a bad file block the rest
    print(f"[warmup] Done — {warmed}/{len(gfa_paths)} GFA file(s) cached.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    global DEFAULT_OUTPUT_DIR, PROJECT_ROOT

    parser = argparse.ArgumentParser(description="Aviary Web Interface")
    parser.add_argument("--output-dir", "-o", help="Pin to a specific aviary output directory")
    parser.add_argument("--port", "-p", type=int, default=8090)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--reload", action="store_true",
                        help="Auto-reload server when source files change (development)")
    args = parser.parse_args()

    PROJECT_ROOT = os.getcwd()
    if args.output_dir:
        DEFAULT_OUTPUT_DIR = os.path.abspath(args.output_dir)

    # Find the first free port starting from the requested one
    port = args.port
    for candidate in range(args.port, args.port + 20):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind((args.host, candidate))
                port = candidate
                break
            except OSError:
                continue
    else:
        print(f"Error: no free port found in range {args.port}–{args.port + 19}")
        return

    if port != args.port:
        print(f"Port {args.port} is in use — using {port} instead.")

    print(f"Aviary Web Interface →  http://localhost:{port}")
    print(f"Project root:           {PROJECT_ROOT}")
    if args.reload:
        print("Auto-reload:            enabled")
    print("Press Ctrl+C to stop.\n")

    warmup_root = DEFAULT_OUTPUT_DIR or PROJECT_ROOT
    if warmup_root:
        t = threading.Thread(target=_warmup_gfa_cache, args=(warmup_root,), daemon=True)
        t.start()

    app.run(host=args.host, port=port, debug=args.reload, use_reloader=args.reload)


if __name__ == "__main__":
    main()
