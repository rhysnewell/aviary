# Changelog

## Unreleased

### Added

- **`aviary recover`/`aviary complete`: SemiBin2 multi-sample binning** — new `--semibin-mode {single,multi}` flag (default `single`, preserving prior behaviour). In `multi` mode, `-a/--assembly` accepts multiple FASTA files; they are concatenated via `SemiBin2 concatenate_fasta`, reads from every sample are mapped back onto the concatenation, and `SemiBin2 multi_easy_bin` co-bins across samples using cross-sample abundance correlation. Output bins are written per-sample-prefixed into `data/semibin_bins/output_bins/`. Passing multiple assemblies without `--semibin-mode multi` is now a hard error; passing `--semibin-mode multi` with only one assembly warns but proceeds.

## v0.13.1 - 2026-07-08

Patch release focused on repairing database downloads (`aviary configure --download`).

---

### Fixed

- **`aviary configure --download` no longer requires read inputs** — `download_databases` added to `SUBCOMMANDS_WITHOUT_READS`; previously failed with "both long_reads and short_reads_1 are set to none"
- **eggNOG database download** — `eggnogdb.embl.de` was decommissioned; files are now fetched from `eggnog5.embl.de`
- **Metabuli GTDB database download** — upstream relocated the tarball to an `archive/` path. The old command 404'd but exited 0, silently leaving an empty database; it now downloads the archived index directly and fails loudly on error
- **CheckM2 database download** — unsets `CHECKM2DB` and runs under `bash -e -o pipefail` so download failures are no longer swallowed

### Changed

- **pixi 0.71+ compatibility** — `pixi.toml` migrated to rich platforms (CUDA on platform entries); minimum `pixi` bumped to `>=0.71`; lockfile regenerated

---

## v0.13.0 - 2026-03-31

Forked from [wwood/aviary](https://github.com/wwood/aviary) at v0.12.0 (`myloasm` branch). All changes below are relative to that base.

---

### Changed

- **SingleM updated to v0.21.3** — minimum version bumped from 0.20.3 to 0.21.3; `singlem-appraise` environment unpinned from 0.19.0 now that the v0.20 performance regression for `--genome-fasta-files` input is fixed in v0.21
- **GTDB-Tk updated to v2.7.2** — minimum version bumped from 2.6.1 to 2.7.2. v2.7+ removes the `--skip_ani_screen` flag; replaced by `--place_species`
- **Database updated to GTDB R232** — SingleM metapackage updated to `S6.5.0.GTDB_r232` and GTDB-Tk database updated to `release232`; download URL in `aviary configure --download` updated accordingly
- **Benchmark added to `singlem_appraise` rule** — runtime now recorded in `benchmarks/singlem_appraise.benchmark.txt`

### Added

#### Web Interface (`aviary/web/`) — experimental

An experimental browser-based monitor and results explorer served via Flask. Start with:

```bash
ssh -L 8090:localhost:8090 username@address.com
pixi install -e web
pixi run -e web server --output-dir /path/to/aviary_output
```

Then open `http://localhost:8090` in your browser. Includes a pipeline monitor for current and past runs (with log information), bin quality report, results visualisation, taxonomy tree view (sunburst, Radial and Horizontal cladogram), assembly graph viewer (GFA files), and export functionality.

#### Pipeline

- **myloasm assembler support** — myloasm added as an alternative long-read assembler alongside Flye (Myloasm default)
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
- Removed `--skip_ani_screen` flag in gtdbtk rule (removed in GTDB-Tk v2.7.0)
- Added sleep delays between GPU test submissions to prevent resource contention
- Increased memory allocation for GPU and expensive tests in mqsub commands
- Fixed isolate functionality, with medaka updated to `>=2.2.1` (previously restricted to `<2.1.0` due to [nanoporetech/medaka#566](https://github.com/nanoporetech/medaka/issues/566), fixed in 2.2.x) and dnaapler updated to `>=1.0.0`

---

### Changed (pipeline only)

- Assembly rules now retain GFA graph files for downstream use
- Test output directory naming split for CPU and GPU tests to prevent log/data overwrites
- Log completion message added after refinery process finishes
