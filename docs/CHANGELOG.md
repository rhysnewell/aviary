# Changelog

## v0.13.0 - 2026-03-31

Forked from [wwood/aviary](https://github.com/wwood/aviary) at v0.12.0 (`myloasm` branch). All changes below are relative to that base.

---

### Added

#### Web Interface (`aviary/web/`)

A full browser-based monitor and results explorer served via Flask with a self-contained pixi environment. Start with:

```bash
ssh -L 8090:localhost:8090 username@address.com
pixi run -e web server --output-dir /path/to/aviary_output
```

Then open `http://localhost:8090` in your browser.

- **Landing page** (`/`) — animated splash screen with navigation nodes linking to all pages
- **Pipeline Monitor** (`/dashboard`) — live Snakemake job tree grouped by sample and run date; log viewer with error highlighting, attempt selector, and search; auto-refresh every 60 seconds; output directory switcher; dismiss failed run banners
- **Bin Quality Report** (`/dashboard` → tab) — per-bin completeness, contamination, strain heterogeneity, N50, size, GC, and GTDB taxonomy; assembly stats per sample/assembler; bins sorted HQ → MQ → LQ then by completeness
- **Results Visualisation** (`/graph`) — six interactive chart types (bar, scatter, histogram, phylum stacked bar, heatmap, donut); 11 colour palette options; font and axis label customisation; SVG and PNG export
- **Tree View** (`/view`) — interactive taxonomy tree from GTDB-Tk bin classifications with three layout modes:
  - **Sunburst** (default) — concentric ring chart; click any segment to zoom into that clade; click the centre to go back; scroll or use topbar buttons to zoom in/out
  - **Radial (Circular)** — circular cladogram with phylum-level colour bands
  - **Horizontal Cladogram** — left-to-right dendrogram for reading deep taxonomy labels
  - Quality and sample filters; taxonomy search; label size and font controls; SVG and PNG export
- **Contig Assembly Graph** (`/assembly`) — per-sample GFA assembly graph viewer for megahit and metaSPAdes outputs
- **Export Results** (`/export`) — filtered bin table download as CSV or TSV; column selection grouped by category; quality filters; live preview of first 50 rows
- **Documentation** (`/docs`) — inline Aviary pipeline documentation, citations, and external resource links
- **User Guide** (`/guide`) — full in-browser guide covering all pages, controls, HPC/SSH tunnel setup, and API reference

#### Web API

All endpoints accept an optional `?root=` parameter to override the default output directory:

| Endpoint | Description |
|---|---|
| `GET /api/structure` | Sample/assembler job tree grouped by run date |
| `GET /api/status` | Detailed status for a single output directory |
| `GET /api/summary` | Bins and assembly stats for all runs |
| `GET /api/logs` | Available Snakemake log files for a directory |
| `GET /api/job_log` | Log file content and available retry attempts |
| `GET /api/benchmark` | Benchmark timing data for a job |
| `GET /api/output_dirs` | All discovered output directories under the root |
| `GET /api/gfa_stats` | Parsed GFA assembly graph statistics for an output directory |
| `GET /api/gfa_available` | Map of output directories to GFA file availability |
| `GET /api/taxonomy_tree` | Taxonomy tree — SingleM condensed profile when available, falls back to GTDB-Tk |
| `GET /api/singlem_status` | Count of output directories with a SingleM condensed profile |
| `GET /api/phylo_newick` | GTDB-Tk Newick tree with MAG annotations for the tree view |

#### Pipeline

- **myloasm assembler support** — myloasm added as an alternative long-read assembler alongside Flye (`--long-read-assembler myloasm`)
- **GFA graph generation for short and long read assembly** — assembly graphs produced and retained for use in the assembly graph viewer
- **`skip_reads_check` parameter** — added to `template_config.yaml` and config handling to support running subcommands without providing reads
- **SingleM metapackage support in integration tests**

#### Environment

- **Web environment in main `pixi.toml`** — flask added as `[feature.web]` so a single `pixi install` covers both the pipeline and the web interface
- **`server` task** — `pixi run -e web server` starts the web interface directly

---

### Fixed

- Updated symlinks to be working correctly
- Added diamond dependency to das-tool environment
- Added gzip dependencies so packages can read compressed files
- Updated `pbsim.fq.gz` binary file
- Fixed quickbin rule to use relative path for `quickbin.sh` after BBTool update
- Added conditional check for `skip_reads_check` to prevent errors when no reads are provided
- Normalised qnames in PAF output to match original fastq read names
- Handle unexpected long read types by raising a clear exception
- Ensured `bam_cache` directory is created in `get_coverage` before use
- Added memory allocation parameter to quickbin rule to prevent crashes
- Added `--only-id` flag to seqkit command in `filter_contigs_by_size` rule
- Refactored log permission handling in `onsuccess`/`onerror` to check file existence first
- Added temporary directory creation for GFA conversion in `assemble_short_reads`
- Unset CUDA environment variables in taxvamb, semibin, and comebin rules for improved compatibility
- Added file locking to NanoPlot command for improved concurrency handling
- Fixed `OUTPUT_DIR` path to use `PIXI_PROJECT_ROOT` environment variable
- Added pandas import to `prepare_binning_files_gather` rule
- Updated blas and blas-devel package versions in `pixi.lock` for compatibility
- Refactored unrefined binners list formatting
- Removed unnecessary `--no-assign-taxonomy` flag from SingleM commands
- Handle `"none"` input for read lists in `ReadContainer` initialisation
- Added scratch directory to gtdbtk rule for better temporary file and memory handling
- Added sleep delays between GPU test submissions to prevent resource contention
- Increased memory allocation for GPU and expensive tests in mqsub commands

---

### Changed (pipeline only)

- Assembly rules now retain GFA graph files for downstream use
- Test output directory naming split for CPU and GPU tests to prevent log/data overwrites
- Log completion message added after refinery process finishes
