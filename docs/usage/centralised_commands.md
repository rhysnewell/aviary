---
title: Centralised commands
---

# Centralised commands

These options are not specific to any one subcommand. They're defined once in `aviary.py` as a shared parent parser (`base_group`, plus the `-w`/`--workflow` argument added separately per subcommand) and attached to every subcommand via argparse's `parents=` mechanism — that's why `--clean`, `--dry-run`, `--build`, etc. behave identically whether you're running `aviary recover`, `aviary annotate`, or any other subcommand.

**Note on `bird_tool_utils`:** the flags below are Aviary's own, not part of `bird_tool_utils`. What *does* come from `bird_tool_utils` is the surrounding CLI machinery: the `BirdArgparser` wrapper class, the short-help/`--full-help` split (see the gotcha below), and the overall `--help` formatting. Don't confuse "shared across subcommands" (this page) with "provided by bird_tool_utils" (the CLI plumbing) — they're two different things that happen to overlap here.

## Gotchas

**`--help` vs `--full-help`.** For `assemble`, `recover`, and `complete`, plain `-h`/`--help` prints only a short description and examples — none of the flags on this page (or the subcommand's own flags) are shown. You need `--full-help`/`--full_help` to get the complete flag listing. The other subcommands (`annotate`, `cluster`, `build`, `isolate`, `configure`) use plain argparse help, so `--help` already shows everything for them.

**Boolean flags need an explicit value to turn off.** `--clean`, `--dry-run`, `--strict`, `--request-gpu`, `--build`, and `--build-gpu` all use `nargs='?', const=True` with a custom `str2bool` parser. That means:
- Bare `--clean` (no value) → `True`.
- To disable, you must pass an explicit falsy value: `--clean False`, `--clean no`, `--clean 0`, etc. (accepted values: `yes/true/t/y/1` and `no/false/f/n/0`, case-insensitive.) There's no `--no-clean` form.
- `--clean` in particular defaults to `True`, so if you want to keep intermediate files (e.g. for debugging or to resume from a partial run), you must explicitly pass `--clean False` — omitting the flag does not do this.

**`--build`/`--build-gpu` short-circuit the pipeline.** Passing either causes aviary to build the dependency Conda/Snakemake environments for that subcommand and then exit — it will not run the actual workflow in the same invocation.

**`--workflow` is shared in mechanism, not in default.** Every subcommand exposes `-w`/`--workflow` (the snakemake target rule to run), but each has its own default target and, for `configure`, the help text is even suppressed:

| Subcommand | Default workflow target |
| --- | --- |
| `assemble` | `complete_assembly_with_qc` |
| `recover` | `recover_mags` |
| `annotate` | `annotate` |
| `cluster` | `complete_cluster` |
| `build` | `build` |
| `complete` | `get_bam_indices recover_mags annotate` |
| `isolate` | `dnaapler` |
| `configure` | `download_databases` (help hidden) |

## Performance options

**`-t`**, **`--max-threads`** INT

  Maximum number of threads given to any particular process. If `max_threads` > `n_cores`, `n_cores` is bumped up to match. [default: 8]

**`-p`**, **`--pplacer-threads`** INT

  Threads given to pplacer. Values above `--max-threads` are scaled down to equal it. [default: 8]

**`-n`**, **`--n-cores`** INT

  Maximum number of cores available for use. Set to a multiple of `max_threads` to allow multiple processes in parallel. [default: 16]

**`-m`**, **`--max-memory`** INT

  Maximum memory available, in gigabytes. [default: 250]

**`--local-cores`** INT

  Maximum cores available locally. Only relevant when submitting to a cluster (see `--snakemake-profile`), in which case `--n-cores` restricts cores requested per submitted job. [default: 16]

## Output options

**`-o`**, **`--output`** DIR

  Output directory. [default: `./`]

**`--tmpdir`** DIR

  Directory used for temporary files. Aliases: `--tempdir`, `--tmp-dir`, `--tmp`, `--temp`, `--temp-dir`. If not specified, the `TMPDIR` environment variable is used. Can also be configured via the `configure` subcommand. [default: none]

## Misc options

**`--request-gpu`**

  Request a GPU for the pipeline. Only takes effect when run on a cluster (see `--snakemake-profile`). [default: false]

**`--strict`**

  Ensure each binner completes successfully, rather than skipping ones that fail. [default: false — skip failing binners]

**`--default-resources`** STRING

  Snakemake resources passed through as-is — see the [Snakemake docs](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html?highlight=resources#standard-resources). Note: `tmpdir` is handled by aviary's own `--tmpdir` flag, not through this option.

**`--snakemake-profile`** PROFILE

  Snakemake profile for cluster submission (see the [Snakemake docs](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)). Create the profile at `~/.config/snakemake/[CLUSTER_PROFILE]/config.yaml`. Requires `cluster`, `cluster-status`, `jobs`, `cluster-cancel` to be set. See the Guides section for HPC usage.

**`--cluster-retries`** INT

  Number of times to retry a failed job when using cluster submission (see `--snakemake-profile`). [default: 0]

**`--dry-run`**

  Perform a snakemake dry run: tests workflow order and conda environments without actually running anything.

**`--clean`**

  Clean up all temporary files — most BAM files and any FASTQ files generated from read filtering. Setting to `False` is the equivalent of snakemake's `--notemp`. Useful when running only part of a workflow, since it avoids deleting files needed by later parts. [default: true]

  > Not cleaning makes reruns faster, but will incur the wrath of your sysadmin.

**`--build`**

  Build Aviary's dependency environments, then exit — does not run the actual subcommand pipeline. [default: no]

**`--build-gpu`**

  Same as `--build`, but also builds GPU-enabled environments. Requires a GPU to be present. [default: no]

**`--download`** [DATABASE ...]

  Download the requested databases: `gtdb`, `eggnog`, `singlem`, `checkm2`, `metabuli`. If no arguments are given, all databases are downloaded. [default: none downloaded]

**`--rerun-triggers`** [TRIGGER ...]

  Which kinds of modifications should trigger a rule to rerun: `mtime`, `params`, `input`, `software-env`, `code`. [default: `mtime`]

**`--snakemake-cmds`** STRING

  Additional commands passed through to snakemake as a single string, e.g. `--snakemake-cmds "--print-compilation True"`. Most `snakemake -h` commands are valid, but some may clash with commands aviary supplies directly — check for conflicts before using.

**`-w`**, **`--workflow`** TARGET [TARGET ...]

  The snakemake target rule(s) to run. Default varies per subcommand — see the gotcha above.
