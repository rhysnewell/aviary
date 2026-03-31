#!/usr/bin/env python3
"""
Generate a self-contained HTML results report for an Aviary output directory.

Usage:
    python generate_report.py --output-dir /path/to/aviary/output
    python generate_report.py --output-dir /path/to/aviary/output --report ~/my_report.html
"""

import os
import re
import json
import argparse
from pathlib import Path
from datetime import datetime


# ---------------------------------------------------------------------------
# Data parsers
# ---------------------------------------------------------------------------

def parse_assembly_stats(stats_path):
    """Parse www/assembly_stats.txt into a dict."""
    result = {}
    p = Path(stats_path)
    if not p.exists():
        return None
    try:
        text = p.read_text()
    except OSError:
        return None

    patterns = {
        "total_length_mb": r"Main genome scaffold sequence total:\s+([\d.]+) Mb",
        "contig_count":    r"Main genome scaffold total:\s+([\d,]+)",
        "n50":             r"Main genome scaffold N/L50:\s+[\d,]+/([\d,]+)",
        "l50":             r"Main genome scaffold N/L50:\s+([\d,]+)/",
        "max_length":      r"Max scaffold length:\s+([\S]+)",
        "gc_content":      r"^GC\s+([\d.]+)",
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

    # GC from the first line of the file (A C G T N IUPAC Other GC GC_stdev)
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


def parse_bin_info(bin_info_path):
    """Parse bins/bin_info.tsv into list of dicts."""
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
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        row = {}
        for i, h in enumerate(headers):
            row[h] = parts[i] if i < len(parts) else ""
        rows.append(row)
    return rows


def parse_coverage(coverage_path):
    """Parse bins/coverm_abundances.tsv if it exists."""
    p = Path(coverage_path)
    if not p.exists():
        return None
    try:
        lines = p.read_text().splitlines()
    except OSError:
        return None
    if len(lines) < 2:
        return None
    headers = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        rows.append(dict(zip(headers, parts)))
    return {"headers": headers, "rows": rows}


def count_final_bins(bins_dir):
    """Count .fna and .fasta files in final_bins/."""
    d = Path(bins_dir)
    if not d.exists():
        return None
    return len(list(d.glob("*.fna")) + list(d.glob("*.fasta")))


def get_snakemake_logs(output_dir):
    """Return sorted list of (datetime_str, Path) for all snakemake logs."""
    log_dir = Path(output_dir) / ".snakemake" / "log"
    if not log_dir.exists():
        return []
    logs = sorted(log_dir.glob("*.snakemake.log"), key=lambda p: p.stat().st_mtime, reverse=True)
    result = []
    for p in logs:
        name = p.name.replace(".snakemake.log", "")
        try:
            dt = datetime.strptime(name[:26], "%Y-%m-%dT%H%M%S.%f")
            label = dt.strftime("%Y-%m-%d %H:%M")
        except ValueError:
            label = name.replace("T", " ")
        result.append((label, p))
    return result


# ---------------------------------------------------------------------------
# Quality helpers
# ---------------------------------------------------------------------------

def quality_tier(completeness, contamination):
    try:
        c = float(completeness)
        x = float(contamination)
    except (ValueError, TypeError):
        return "low"
    if c >= 90 and x <= 5:
        return "high"
    if c >= 50 and x <= 10:
        return "medium"
    return "low"


def parse_taxonomy(classification_str):
    """Parse a GTDB-Tk classification string into a dict of rank→name."""
    if not classification_str or classification_str.strip() in ("", "N/A", "Unclassified"):
        return {}
    ranks = {}
    for part in classification_str.split(";"):
        part = part.strip()
        if "__" in part:
            prefix, name = part.split("__", 1)
            rank_map = {"d": "Domain", "p": "Phylum", "c": "Class",
                        "o": "Order", "f": "Family", "g": "Genus", "s": "Species"}
            rank = rank_map.get(prefix.strip())
            if rank and name.strip():
                ranks[rank] = name.strip()
    return ranks


TAXONOMY_COLORS = {
    "Domain":  ("#dbeafe", "#1d4ed8"),
    "Phylum":  ("#dcfce7", "#15803d"),
    "Class":   ("#fef9c3", "#854d0e"),
    "Order":   ("#fce7f3", "#9d174d"),
    "Family":  ("#ede9fe", "#5b21b6"),
    "Genus":   ("#ffedd5", "#c2410c"),
    "Species": ("#f1f5f9", "#475569"),
}


# ---------------------------------------------------------------------------
# HTML generation
# ---------------------------------------------------------------------------

def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def render_progress_bar(value, max_val, color, label=None):
    try:
        pct = min(100, max(0, float(value) / float(max_val) * 100))
        display = f"{float(value):.1f}%"
    except (ValueError, TypeError, ZeroDivisionError):
        pct = 0
        display = "N/A"
    lbl = label or display
    return f"""<div class="progress-wrap">
      <div class="progress-bar" style="width:{pct:.1f}%;background:{color}"></div>
      <span class="progress-label">{_esc(lbl)}</span>
    </div>"""


def render_taxonomy_tags(classification_str):
    taxa = parse_taxonomy(classification_str)
    if not taxa:
        return '<span class="no-data">No taxonomy</span>'
    tags = []
    for rank in ["Phylum", "Class", "Order"]:
        name = taxa.get(rank)
        if not name:
            continue
        bg, fg = TAXONOMY_COLORS.get(rank, ("#f1f5f9", "#475569"))
        tags.append(f'<span class="tax-tag" style="background:{bg};color:{fg}">{_esc(rank[0])}: {_esc(name)}</span>')
    return " ".join(tags) if tags else '<span class="no-data">No taxonomy</span>'


def render_quality_badge(tier):
    classes = {"high": "badge-high", "medium": "badge-medium", "low": "badge-low"}
    labels  = {"high": "High Quality", "medium": "Medium Quality", "low": "Low Quality"}
    return f'<span class="badge {classes.get(tier,"badge-low")}">{labels.get(tier,"Low Quality")}</span>'


def render_metric_card(value, label, sub=None, color="#16a34a"):
    sub_html = f'<div class="card-sub">{_esc(sub)}</div>' if sub else ""
    return f"""<div class="metric-card">
  <div class="card-value" style="color:{color}">{_esc(str(value))}</div>
  <div class="card-label">{_esc(label)}</div>
  {sub_html}
</div>"""


def render_assembly_box(stats):
    if not stats:
        return '<div class="box"><div class="box-title">Assembly Statistics</div><div class="no-data-box">No assembly stats available</div></div>'

    rows = [
        ("Total Length",    f"{stats.get('total_length_mb', 'N/A')} Mb"),
        ("Contig Count",    f"{stats.get('contig_count', 'N/A'):,}" if isinstance(stats.get('contig_count'), int) else str(stats.get('contig_count', 'N/A'))),
        ("N50",             f"{stats.get('n50', 'N/A'):,} bp" if isinstance(stats.get('n50'), int) else str(stats.get('n50', 'N/A'))),
        ("L50",             str(stats.get('l50', 'N/A'))),
        ("Max Contig",      str(stats.get('max_length', 'N/A'))),
        ("GC Content",      f"{stats.get('gc_content', 'N/A')}%"),
        ("Scaffolds >50kb", str(stats.get('scaffolds_50kb', 'N/A'))),
    ]
    rows_html = "".join(
        f'<tr><td class="stat-key">{_esc(k)}</td><td class="stat-val">{_esc(v)}</td></tr>'
        for k, v in rows
    )
    return f"""<div class="box">
  <div class="box-title">Assembly Statistics</div>
  <table class="stat-table">{rows_html}</table>
</div>"""


def render_bins_table(bins):
    if not bins:
        return '<div class="box"><div class="box-title">Bins</div><div class="no-data-box">No bin data available</div></div>'

    has_circular  = any("Circular contigs" in b for b in bins)
    has_taxonomy  = any(b.get("classification") for b in bins)

    # Sort by completeness desc
    def sort_key(b):
        try:
            return -float(b.get("Completeness", 0))
        except ValueError:
            return 0
    bins_sorted = sorted(bins, key=sort_key)

    header_cells = ["Bin", "Completeness", "Contamination", "Genome Size", "GC%", "Contigs", "N50"]
    if has_circular:
        header_cells.append("Circular")
    if has_taxonomy:
        header_cells.append("Taxonomy")
    header_cells.append("Quality")

    headers_html = "".join(f"<th>{h}</th>" for h in header_cells)

    rows_html = ""
    for b in bins_sorted:
        name = b.get("Name", "")
        # Shorten ugly bin names like metabat2_refined_bins.tsv_binned_contigs.16 → bin.16
        short = re.sub(r".*[._](\d+)$", r"bin.\1", name) if name else name

        comp = b.get("Completeness", "")
        cont = b.get("Contamination", "")
        tier = quality_tier(comp, cont)

        try:
            size_mb = f"{int(b.get('Genome_Size', 0)) / 1_000_000:.2f} Mb"
        except (ValueError, TypeError):
            size_mb = "N/A"

        try:
            gc = f"{float(b.get('GC_Content', 0)) * 100:.1f}%" if float(b.get('GC_Content', 0)) <= 1 else f"{float(b.get('GC_Content', 0)):.1f}%"
        except (ValueError, TypeError):
            gc = "N/A"

        contig_count = b.get("Total_Contigs", "N/A")
        n50 = b.get("Contig_N50", "N/A")
        try:
            n50 = f"{int(n50):,}"
        except (ValueError, TypeError):
            pass

        row = f"""<tr class="tier-{tier}">
          <td title="{_esc(name)}" class="bin-name">{_esc(short)}</td>
          <td>{render_progress_bar(comp, 100, "#16a34a", f"{float(comp):.1f}%" if comp else "N/A")}</td>
          <td>{render_progress_bar(cont, 20,  "#dc2626", f"{float(cont):.1f}%" if cont else "N/A")}</td>
          <td>{_esc(size_mb)}</td>
          <td>{_esc(gc)}</td>
          <td>{_esc(contig_count)}</td>
          <td>{_esc(n50)}</td>"""

        if has_circular:
            circ = b.get("Circular contigs", "N/A")
            row += f"<td>{_esc(str(circ))}</td>"

        if has_taxonomy:
            row += f"<td>{render_taxonomy_tags(b.get('classification',''))}</td>"

        row += f"<td>{render_quality_badge(tier)}</td></tr>"
        rows_html += row

    return f"""<div class="box full-width">
  <div class="box-title">Bins</div>
  <div class="table-scroll">
    <table class="bins-table">
      <thead><tr>{headers_html}</tr></thead>
      <tbody>{rows_html}</tbody>
    </table>
  </div>
</div>"""


def render_coverage_box(coverage):
    if not coverage:
        return '<div class="box"><div class="box-title">Coverage / Abundance</div><div class="no-data-box">No coverage data available</div></div>'

    headers = coverage["headers"]
    rows    = coverage["rows"]
    headers_html = "".join(f"<th>{_esc(h)}</th>" for h in headers)
    rows_html = ""
    for r in rows[:30]:  # cap at 30 rows to keep report size reasonable
        rows_html += "<tr>" + "".join(f"<td>{_esc(r.get(h,''))}</td>" for h in headers) + "</tr>"

    return f"""<div class="box full-width">
  <div class="box-title">Coverage / Abundance</div>
  <div class="table-scroll">
    <table class="bins-table">
      <thead><tr>{headers_html}</tr></thead>
      <tbody>{rows_html}</tbody>
    </table>
  </div>
</div>"""


def render_quality_summary_box(bins):
    if not bins:
        return '<div class="box"><div class="box-title">Quality Summary</div><div class="no-data-box">No bin data</div></div>'
    hq = sum(1 for b in bins if quality_tier(b.get("Completeness"), b.get("Contamination")) == "high")
    mq = sum(1 for b in bins if quality_tier(b.get("Completeness"), b.get("Contamination")) == "medium")
    lq = sum(1 for b in bins if quality_tier(b.get("Completeness"), b.get("Contamination")) == "low")
    total = len(bins)

    def bar(count, color, label):
        pct = count / total * 100 if total else 0
        return f"""<div class="qs-row">
          <span class="qs-label">{label}</span>
          <div class="qs-bar-wrap">
            <div class="qs-bar" style="width:{pct:.1f}%;background:{color}"></div>
          </div>
          <span class="qs-count">{count}</span>
        </div>"""

    return f"""<div class="box">
  <div class="box-title">Quality Summary</div>
  <div class="qs-total">Total: <strong>{total}</strong> bins</div>
  {bar(hq, "#16a34a", "High Quality")}
  {bar(mq, "#d97706", "Medium Quality")}
  {bar(lq, "#dc2626", "Low Quality")}
</div>"""


def render_batch_section(label, output_dir, batch_index):
    """Render a full batch/run section."""
    output_dir = Path(output_dir)

    # Locate key files (handle bins/ and recover/bins/)
    bins_base = output_dir / "recover" / "bins"
    if not bins_base.exists():
        bins_base = output_dir / "bins"

    bin_info_path   = bins_base / "bin_info.tsv"
    final_bins_dir  = bins_base / "final_bins"
    coverage_path   = bins_base / "coverm_abundances.tsv"
    stats_path      = output_dir / "www" / "assembly_stats.txt"

    bins     = parse_bin_info(bin_info_path)
    stats    = parse_assembly_stats(stats_path)
    coverage = parse_coverage(coverage_path)
    mag_count = count_final_bins(final_bins_dir)

    hq = sum(1 for b in bins if quality_tier(b.get("Completeness"), b.get("Contamination")) == "high") if bins else 0
    mq = sum(1 for b in bins if quality_tier(b.get("Completeness"), b.get("Contamination")) == "medium") if bins else 0
    total_mb = stats.get("total_length_mb", "N/A") if stats else "N/A"

    cards = "".join([
        render_metric_card(mag_count if mag_count is not None else "N/A", "MAGs Recovered", "final bins", "#7c3aed"),
        render_metric_card(hq,      "High Quality",    "≥90% · ≤5% cont.", "#16a34a"),
        render_metric_card(mq,      "Medium Quality",  "≥50% · ≤10% cont.", "#d97706"),
        render_metric_card(f"{total_mb} Mb" if total_mb != "N/A" else "N/A", "Assembly Size", "total scaffold length", "#0284c7"),
    ])

    checked = "checked" if batch_index == 0 else ""
    return f"""
<div class="batch-section">
  <input type="checkbox" id="batch-{batch_index}" class="batch-toggle" {checked}>
  <label for="batch-{batch_index}" class="batch-header">
    <span class="batch-arrow">▶</span>
    <span class="batch-date">{_esc(label)}</span>
    <span class="batch-pills">
      <span class="mini-pill purple">{mag_count if mag_count is not None else "?"} MAGs</span>
      <span class="mini-pill green">{hq} HQ</span>
      <span class="mini-pill amber">{mq} MQ</span>
    </span>
  </label>
  <div class="batch-content">
    <div class="metric-cards">{cards}</div>
    <div class="two-col">
      {render_assembly_box(stats)}
      {render_quality_summary_box(bins)}
    </div>
    {render_bins_table(bins)}
    {render_coverage_box(coverage)}
  </div>
</div>"""


def generate_html(output_dir, report_path):
    output_dir = Path(output_dir)
    logs = get_snakemake_logs(output_dir)

    # Derive sample/assembler from directory structure
    parts = output_dir.resolve().parts
    sample   = parts[-1] if len(parts) >= 1 else "Unknown"
    assembler = parts[-2] if len(parts) >= 2 else "Unknown"

    if logs:
        run_date = logs[0][0]
        status_html = '<span class="status-badge complete">Complete</span>'
    else:
        run_date = datetime.fromtimestamp(output_dir.stat().st_mtime).strftime("%Y-%m-%d %H:%M") if output_dir.exists() else "Unknown"
        status_html = '<span class="status-badge unknown">Unknown</span>'

    # Build batch sections
    if logs:
        batches_html = "\n".join(
            render_batch_section(label, output_dir, i)
            for i, (label, _) in enumerate(logs)
        )
    else:
        batches_html = render_batch_section(run_date, output_dir, 0)

    css = """
    *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
           background: #f1f5f9; color: #1e293b; font-size: 14px; line-height: 1.6; }

    /* Header */
    .site-header { background: linear-gradient(135deg, #0f172a 0%, #1e3a5f 100%);
                   color: white; padding: 32px 40px; }
    .site-header h1 { font-size: 1.8rem; font-weight: 700; letter-spacing: -0.02em; }
    .header-meta { margin-top: 8px; display: flex; align-items: center; gap: 16px;
                   flex-wrap: wrap; font-size: 0.85rem; opacity: 0.85; }
    .header-chip { background: rgba(255,255,255,0.12); border-radius: 20px;
                   padding: 3px 12px; font-size: 0.78rem; }
    .status-badge { padding: 4px 14px; border-radius: 20px; font-size: 0.75rem;
                    font-weight: 700; letter-spacing: 0.05em; text-transform: uppercase; }
    .status-badge.complete { background: #bbf7d0; color: #14532d; }
    .status-badge.unknown  { background: #e2e8f0; color: #475569; }

    /* Main content */
    .main { max-width: 1400px; margin: 0 auto; padding: 32px 24px; }

    /* Batch accordion */
    .batch-toggle { display: none; }
    .batch-section { margin-bottom: 16px; }
    .batch-header { display: flex; align-items: center; gap: 12px; cursor: pointer;
                    background: white; border-radius: 10px; padding: 14px 20px;
                    box-shadow: 0 1px 3px rgba(0,0,0,0.08);
                    border: 1px solid #e2e8f0; user-select: none;
                    transition: box-shadow 0.15s; }
    .batch-header:hover { box-shadow: 0 3px 10px rgba(0,0,0,0.1); }
    .batch-arrow { font-size: 0.7rem; color: #94a3b8; transition: transform 0.2s; }
    .batch-toggle:checked ~ .batch-header .batch-arrow { transform: rotate(90deg); }
    .batch-date { font-weight: 600; font-size: 0.95rem; color: #0f172a; }
    .batch-pills { display: flex; gap: 6px; margin-left: auto; }
    .mini-pill { padding: 2px 10px; border-radius: 20px; font-size: 0.72rem; font-weight: 600; }
    .mini-pill.purple { background: #ede9fe; color: #5b21b6; }
    .mini-pill.green  { background: #dcfce7; color: #15803d; }
    .mini-pill.amber  { background: #fef3c7; color: #92400e; }
    .batch-content { display: none; padding-top: 16px; }
    .batch-toggle:checked ~ .batch-content { display: block; }

    /* Metric cards */
    .metric-cards { display: grid; grid-template-columns: repeat(4, 1fr); gap: 16px; margin-bottom: 16px; }
    .metric-card { background: white; border-radius: 10px; padding: 20px 24px;
                   box-shadow: 0 1px 3px rgba(0,0,0,0.08); border: 1px solid #e2e8f0; }
    .card-value { font-size: 2rem; font-weight: 700; line-height: 1.1; }
    .card-label { font-size: 0.8rem; font-weight: 600; color: #64748b;
                  text-transform: uppercase; letter-spacing: 0.06em; margin-top: 4px; }
    .card-sub   { font-size: 0.72rem; color: #94a3b8; margin-top: 2px; }

    /* Two-column layout */
    .two-col { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; margin-bottom: 16px; }

    /* Boxes */
    .box { background: white; border-radius: 10px; padding: 20px 24px;
           box-shadow: 0 1px 3px rgba(0,0,0,0.08); border: 1px solid #e2e8f0; }
    .box.full-width { margin-bottom: 16px; }
    .box-title { font-weight: 700; font-size: 0.85rem; text-transform: uppercase;
                 letter-spacing: 0.07em; color: #475569; margin-bottom: 14px;
                 padding-bottom: 10px; border-bottom: 1px solid #f1f5f9; }
    .no-data-box { color: #94a3b8; font-size: 0.85rem; font-style: italic; padding: 12px 0; }

    /* Stat table */
    .stat-table { width: 100%; border-collapse: collapse; }
    .stat-table tr { border-bottom: 1px solid #f8fafc; }
    .stat-table tr:last-child { border-bottom: none; }
    .stat-key { color: #64748b; font-size: 0.82rem; padding: 6px 0; width: 50%; }
    .stat-val { font-weight: 600; font-size: 0.82rem; text-align: right; color: #0f172a; }

    /* Quality summary */
    .qs-total { font-size: 0.85rem; color: #64748b; margin-bottom: 12px; }
    .qs-row { display: flex; align-items: center; gap: 10px; margin-bottom: 8px; }
    .qs-label { width: 110px; font-size: 0.78rem; color: #475569; flex-shrink: 0; }
    .qs-bar-wrap { flex: 1; height: 10px; background: #f1f5f9; border-radius: 99px; overflow: hidden; }
    .qs-bar { height: 100%; border-radius: 99px; transition: width 0.4s; }
    .qs-count { width: 28px; text-align: right; font-weight: 700; font-size: 0.82rem; color: #0f172a; }

    /* Bins table */
    .table-scroll { overflow-x: auto; }
    .bins-table { width: 100%; border-collapse: collapse; font-size: 0.8rem; }
    .bins-table th { text-align: left; padding: 8px 10px; background: #f8fafc;
                     color: #64748b; font-size: 0.7rem; font-weight: 700;
                     letter-spacing: 0.07em; text-transform: uppercase;
                     border-bottom: 2px solid #e2e8f0; white-space: nowrap; }
    .bins-table td { padding: 7px 10px; border-bottom: 1px solid #f1f5f9; vertical-align: middle; }
    .bins-table tbody tr:hover td { background: #f8fafc; }
    .bins-table tbody tr:last-child td { border-bottom: none; }
    .bin-name { font-family: "SF Mono", Consolas, monospace; font-size: 0.72rem;
                max-width: 120px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }

    /* Progress bars */
    .progress-wrap { position: relative; height: 18px; background: #f1f5f9;
                     border-radius: 99px; overflow: hidden; min-width: 80px; }
    .progress-bar  { height: 100%; border-radius: 99px; opacity: 0.8; }
    .progress-label { position: absolute; inset: 0; display: flex; align-items: center;
                      justify-content: center; font-size: 0.68rem; font-weight: 700;
                      color: #0f172a; pointer-events: none; }

    /* Quality badges */
    .badge { display: inline-block; padding: 2px 10px; border-radius: 20px;
             font-size: 0.68rem; font-weight: 700; letter-spacing: 0.04em;
             text-transform: uppercase; white-space: nowrap; }
    .badge-high   { background: #dcfce7; color: #14532d; }
    .badge-medium { background: #fef3c7; color: #92400e; }
    .badge-low    { background: #fee2e2; color: #991b1b; }

    /* Row tier shading */
    tr.tier-high   td { background: rgba(220,252,231,0.3); }
    tr.tier-medium td { background: rgba(254,243,199,0.3); }
    tr.tier-low    td { background: rgba(254,226,226,0.2); }
    tr.tier-high:hover   td, tr.tier-medium:hover td, tr.tier-low:hover td { background: #f8fafc !important; }

    /* Taxonomy tags */
    .tax-tag { display: inline-block; padding: 2px 8px; border-radius: 20px;
               font-size: 0.68rem; font-weight: 600; margin-right: 3px; white-space: nowrap; }
    .no-data { color: #cbd5e1; font-size: 0.75rem; font-style: italic; }

    /* Footer */
    .footer { text-align: center; color: #94a3b8; font-size: 0.75rem; padding: 32px 0; }

    @media (max-width: 900px) {
      .metric-cards { grid-template-columns: repeat(2, 1fr); }
      .two-col { grid-template-columns: 1fr; }
    }
    """

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Aviary Results — {_esc(sample)}</title>
<style>{css}</style>
</head>
<body>

<header class="site-header">
  <div style="display:flex;align-items:center;gap:16px;flex-wrap:wrap;">
    <h1>Aviary Run Results</h1>
    {status_html}
  </div>
  <div class="header-meta">
    <span class="header-chip">Sample: {_esc(sample)}</span>
    <span class="header-chip">Assembler: {_esc(assembler)}</span>
    <span class="header-chip">Latest run: {_esc(run_date)}</span>
    <span class="header-chip">{len(logs)} batch{"es" if len(logs) != 1 else ""}</span>
  </div>
</header>

<div class="main">
  {batches_html}
</div>

<div class="footer">Generated {datetime.now().strftime("%Y-%m-%d %H:%M")} · Aviary Results Report</div>

</body>
</html>"""

    Path(report_path).write_text(html)
    print(f"Report written to {report_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate Aviary HTML results report")
    parser.add_argument("--output-dir", "-o", required=True,
                        help="Aviary output directory (e.g. megahit/ERR599108/)")
    parser.add_argument("--report", "-r", default=str(Path.home() / "aviary_results.html"),
                        help="Output HTML file path (default: ~/aviary_results.html)")
    args = parser.parse_args()

    if not Path(args.output_dir).exists():
        print(f"Error: output directory not found: {args.output_dir}")
        return 1

    generate_html(args.output_dir, args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
