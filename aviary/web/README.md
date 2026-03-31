# Aviary Web Interface

A live run monitor and results explorer for Aviary/Snakemake workflows. Powered by Flask and served via a self-contained pixi environment in `aviary/web/`.

---

## Pages

| Page | URL | Description |
|---|---|---|
| **Landing Page** | `/` | Animated splash screen with navigation nodes linking to all pages. |
| **Pipeline Monitor** | `/dashboard` | Live Snakemake job tree grouped by sample and date. View log output, retry attempts, and error-highlighted lines per job. |
| **Bin Quality Report** | `/dashboard` → Bin Quality Report tab | Per-sample bin quality table (HQ / MQ / LQ) and assembly stats. Bins sorted by quality tier then completeness. |
| **Results Visualisation** | `/graph` | Generate interactive charts from results: bar, scatter, histogram, phylum stacked bar, heatmap, and donut. Full colour palette selection, font controls, and axis label customisation. SVG and PNG download. |
| **Phylogenetic Identification** | `/phylo` | Interactive phylogenetic tree built from GTDB-Tk classification. Radial and horizontal layouts. Configurable initial view depth, label visibility, font size/style, and per-phylum colour coding. |
| **Contig Assembly Graph** | `/assembly` | Per-sample assembly graph viewer. Browse the GFA-format assembly graphs produced by megahit and metaSPAdes. |
| **Export Results** | `/export` | Download bin tables as CSV or TSV with column selection and quality filters. Live preview before download. |
| **Documentation** | `/docs` | Inline Aviary documentation, citations, and links to external resources. |
| **User Guide** | `/guide` | How to use the web interface — this guide served in-browser. |

---

## Starting the Server

### First time — install the environment

```bash
cd /path/to/aviary/aviary/web
~/.pixi/bin/pixi install
```

This creates a self-contained environment with Python and Flask. Only needed once.

### Start the server

> **On an HPC?** Open the SSH tunnel first (see [Setup on a Shared HPC](#setup-on-a-shared-hpc) below), then come back here to start the server.

```bash
pixi run -e web server --output-dir /path/to/aviary_output
```

Then open `http://localhost:8090` in your browser.

### Arguments

| Argument | Default | Description |
|---|---|---|
| `--output-dir` / `-o` | current working directory | Root directory to scan for Aviary outputs. The server traverses any subdirectory depth automatically (including `<commit-hash>/assembler/sample` layouts). |
| `--port` / `-p` | `8090` | Port to listen on |
| `--host` | `127.0.0.1` | Bind address. Keep as `127.0.0.1` on shared HPCs to avoid port conflicts with other users. |
| `--reload` | off | Auto-restart when `server.py` or templates are saved. Useful when editing the interface. |

---

## Setup on a Shared HPC

Two terminal windows are required — both opened locally (not already inside an HPC session).

**Window 1 — SSH tunnel (keep open for the whole session):**
```powershell
ssh -L 8090:localhost:8090 <username>@<hpc-address>
```

**Window 2 — start the server on the HPC:**
```powershell
ssh <username>@<hpc-address>
cd /path/to/aviary/aviary/web
pixi run -e web server --output-dir /path/to/aviary_output
```

Then open your browser at `http://localhost:8090`.

> Keep Window 1 open for as long as you want to use the interface. Closing it drops the tunnel and the browser will lose connection.

---

## Pipeline Monitor

The main view groups sample cards by **run date** (derived from Snakemake log timestamps). Each date group shows an overall progress bar for that batch.

Within each date group, every sample has a card showing:
- An **assembler pill** for each assembler (megahit / metaspades) with a live progress bar
  - Teal = in progress · Green = completed · Red = failed
- Expanding an assembler pill reveals every **individual job** (rule) with its status, duration, threads, memory, and a link to its log

The sidebar overview shows **sample count, running, completed, and failed** counts — scoped to the most recent run so historical runs don't inflate the numbers.

### Dismissing Failed Runs

When a run fails, a red banner appears on the sample card listing the failed rule(s). If you have already acknowledged the failure and don't want to see the indicator, click **Dismiss** on the banner. This hides the banner and turns the assembler pill dot grey. The dismissal is saved in `localStorage` and tied to the specific log file — if the sample is re-run and fails again, the red indicator reappears automatically.

### Log Viewer

Click any job to open its log in the right panel. The log viewer:
- **Auto-scrolls to the first error** (`Error`, `Failed`, `Exception`, `Traceback`, etc.) and highlights matching lines in red. If no errors are found it scrolls to the bottom.
- **Attempt selector** — if a job was retried, pills appear above the log (`Attempt 1`, `Attempt 2`, …) to switch between log files without leaving the page.
- **Search** — type to filter lines; matching text is highlighted. Use the × button or clear the field to restore full log view (error highlights are reapplied).

### Auto-Refresh

The page refreshes every 60 seconds. A live countdown (`in Xs`) appears next to the pulse dot in the sidebar so you always know when the next update is coming. You can also hit the manual refresh button at any time.

### Output Directory Switcher

The **directory switcher** in the top bar lets you hot-swap between multiple Aviary project roots without restarting the server. Type a path and press Go, or select from recently-used roots. The last 5 roots are remembered across sessions via `localStorage`.

---

## Bin Quality Report

Navigate to `/dashboard` and toggle from **Pipeline Monitor** to **Bin Quality Report** in the top bar to view final results once runs are complete.

- **Bin table** — per-bin completeness, contamination, strain heterogeneity, N50, size, GC, and GTDB taxonomy. Sorted HQ → MQ → LQ, then by completeness descending within each tier.
  - HQ: ≥ 90% completeness, ≤ 5% contamination
  - MQ: ≥ 50% completeness, ≤ 10% contamination
  - LQ: everything else
- **Assembly stats** — total length, N50, contig count, GC content per sample/assembler.

---

## Results Visualisation (`/graph`)

A dedicated page for generating charts from your results. Six chart types are available:

| Chart | What it shows |
|---|---|
| **Bar** | HQ/MQ/LQ counts, total MAGs, N50, or assembly size per sample — grouped or stacked |
| **Scatter** | Any two numeric bin fields against each other; colour by tier, assembler, sample, or domain; optional size encoding |
| **Histogram** | Distribution of any bin quality field with adjustable bin count; threshold lines for completeness/contamination tiers |
| **Phylum Stacked Bar** | GTDB taxonomy breakdown at domain, phylum, or class level; top-N taxa + Other |
| **Heatmap** | Sample × assembler grid coloured by HQ MAG count, mean completeness, or total wall time |
| **Donut** | Proportional breakdown by quality tier, domain, or assembler |

Charts can be downloaded as **SVG** (vector) or **PNG** (2× scale raster) using the buttons in the top bar.

### Customise Panel

The **Customise** section lets you adjust the appearance of the current chart without regenerating it:

- **Colour palette** — 11 palettes available. The first three colours in each palette map to HQ / MQ / LQ respectively, so quality-tier charts respond immediately to palette changes. Palettes include Classic (green/orange/red), Aviary, Bold, Vivid, Ocean, Sunset, Tropical, Accessible, Forest, Mono, and Amethyst.
- **Font family** — override the default system font with Serif, Monospace, Georgia, or others.
- **X / Y axis labels** — type a label to override an existing axis title or inject a new one on charts that don't have one by default. Applies to the current SVG immediately.

---

## Phylogenetic Identification (`/phylo`)

An interactive phylogenetic tree built from GTDB-Tk classifications across all bins.

### Layout

| Layout | Description |
|---|---|
| **Radial (Circular)** | Classic circular cladogram. Phylum-level clade arcs with colour bands. Best for a full overview. |
| **Horizontal Cladogram** | Left-to-right dendrogram. Better for reading taxonomy labels on deep trees. |

### Controls

- **Initial View Depth** — choose how deep the tree expands on first load: Phylum, Class (default), Order, Family, Genus, or fully expanded. Changing this resets the view.
- **Include Names** — tick to show labels at the chosen depth without expanding branches. The tree structure is unchanged; only label visibility toggles.
- **Quality Filter** — show/hide HQ, MQ, or LQ bins to focus on a specific quality tier.
- **Sample Filter** — include/exclude individual samples from the tree.
- **Search Taxonomy** — type a name to highlight matching nodes across the tree.
- **Display** — label size slider (0.5× to 2.5×) and font style selector. Both update the tree in-place without re-rendering.

### Interaction

- **Click** any node to expand or collapse that branch.
- **Expand All / Collapse** buttons in the top bar reset the entire tree.
- **Drag and scroll** to pan and zoom.
- **Hover** any node for a tooltip showing rank, MAG count, and HQ/MQ/LQ breakdown.

### Download

Export the current tree view as **SVG** or **PNG** using the Download button in the top bar.

---

## Export Results (`/export`)

A two-panel export page:

- **Left panel** — configure your export:
  - Select which samples to include (accordion grouped by date)
  - Choose which columns to include via checkboxes (grouped: identity, quality, assembly, taxonomy)
  - Apply quality filters (minimum completeness, maximum contamination, quality tier)
  - Toggle output format: CSV or TSV
- **Right panel** — live preview of the first 50 rows matching your selection
- **Download** button in the top bar exports the full filtered table

---

## API Reference

The server exposes the following JSON endpoints. All endpoints accept an optional `?root=<path>` query parameter to override the default output directory.

| Endpoint | Method | Description |
|---|---|---|
| `/api/structure` | GET | Scans for `.snakemake/log/` folders, parses the latest Snakemake log; returns sample/assembler status grouped by run date and the resolved `root` path |
| `/api/status` | GET | Detailed status for a single output directory (`?output_dir=`) |
| `/api/summary` | GET | Final results — bins, assembly stats, SingleM data for all runs |
| `/api/logs` | GET | Lists available Snakemake logs for a directory |
| `/api/job_log` | GET | Content of an individual job log (`?log_path=`) plus `available_attempts` (sibling attempt log paths) |
| `/api/benchmark` | GET | Benchmark timing data for a job (`?benchmark_path=`) |
| `/api/output_dirs` | GET | Lists all discovered output directories under the root |
| `/api/gfa_stats` | GET | Parsed GFA assembly graph statistics for a single output directory (`?output_dir=`). Results are cached by file mtime. |
| `/api/gfa_available` | GET | Map of all discovered output directories to a boolean indicating whether `assembly/assembly_graph.gfa` exists |
| `/api/taxonomy_tree` | GET | Taxonomy breakdown tree — prefers SingleM condensed profile data when available, falls back to GTDB-Tk classification strings |
| `/api/singlem_status` | GET | Count of output directories that have a SingleM condensed profile (`total` and `done` keys) |
| `/api/phylo_newick` | GET | GTDB-Tk Newick tree with MAG annotations for the phylogenetic viewer (`?output_dir=` or auto-detected from `root`) |

---

## How It Works

The server is a **Flask** application with a self-contained pixi environment (`aviary/web/pixi.toml`). It serves static HTML templates and a set of JSON API endpoints. All charts, the phylogeny tree, and the export preview are generated entirely in the browser using inline SVG and D3.js — no bundled charting libraries required.

Output directories are discovered by recursively searching for `.snakemake/log/` under the project root. This means any nesting depth is supported automatically.

---

## Development

```bash
cd aviary/web

# Start with auto-reload (server restarts on any file save)
pixi run -e web server --output-dir /path/to/test/output --reload
```

Templates live in `aviary/web/templates/`. The server has no build step — edit a template and reload the page.
