# Changelog

## v0.13.0 - 2026-03-31

Forked from [wwood/aviary](https://github.com/wwood/aviary) at v0.12.0 (`myloasm` branch). All changes below are relative to that base.

---

### Added

#### Web Interface (`aviary/web/`) — experimental

An experimental browser-based monitor and results explorer served via Flask. Start with:

```bash
ssh -L 8090:localhost:8090 username@address.com
pixi run -e web server --output-dir /path/to/aviary_output
```

Then open `http://localhost:8090` in your browser. Includes a pipeline monitor for current and past runs (with log information), bin quality report, results visualisation, taxonomy tree view (sunburst, Radial and Horizontal cladogram), assembly graph viewer (GFA files), and export functionality.

#### Pipeline

- **myloasm assembler support** — myloasm added as an alternative long-read assembler alongside Flye (`--long-read-assembler myloasm`)
- **GFA graph generation for short and long read assembly** — assembly graphs produced and retained for use in the assembly graph viewer
- **`skip_reads_check` parameter** — added to `template_config.yaml` and config handling to support running subcommands without providing reads
- **SingleM metapackage support in integration tests**
- **FastQC replaced with RastQC** — RastQC is implemented as a drop-in replacement for FastQC, providing equivalent short-read QC reporting. All files are processed in a single command invocation for clean, readable log output.

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
- Fixed isolate functionality, with updated medaka

---

### Changed (pipeline only)

- Assembly rules now retain GFA graph files for downstream use
- Test output directory naming split for CPU and GPU tests to prevent log/data overwrites
- Log completion message added after refinery process finishes
