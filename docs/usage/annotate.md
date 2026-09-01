---
title: aviary annotate
---

# aviary annotate

Annotate a given set of MAGs using EggNOG and GTDB-Tk.

```bash
aviary annotate --genome-fasta-directory input_bins/
```

> This subcommand also accepts the common workflow, resource, output and
> execution options described under [Shared options](centralised_commands.md).

## Examples

### Full annotation (GTDB-tk + EggNOG)

```bash
aviary annotate --genome-fasta-directory input_bins/
```

The default `annotate` target runs GTDB-tk taxonomy and EggNOG functional annotation together
(CheckM2 is run separately, as part of the `recover`/`complete` binning pipeline rather than
here). Point at specific database locations with `--gtdb-path`/`--eggnog-db-path` if they
aren't already set via `aviary configure`:
```bash
aviary annotate --genome-fasta-directory input_bins/ --gtdb-path /path/to/gtdb/
```

### Run a single annotator

Use `-w`/`--workflow` to target one annotation step instead of the full set — useful for
re-running just the step that failed, or when you only need one kind of annotation:

Taxonomy only (GTDB-tk):
```bash
aviary annotate --genome-fasta-directory input_bins/ -w gtdbtk
```

Functional annotation only (EggNOG):
```bash
aviary annotate --genome-fasta-directory input_bins/ -w eggnog
```

CheckM2 quality assessment is not run through `aviary annotate` — it runs as part of the
`recover`/`complete` binning pipeline, where it has access to the intermediate binning outputs
it depends on.

### A different FASTA extension

`--fasta-extension` defaults to `fna`; set it to match your files if they use something else
(e.g. bins produced outside aviary):
```bash
aviary annotate --genome-fasta-directory input_bins/ --fasta-extension fa
```

## Input options

**`-d`**, **`--genome-fasta-directory`** DIR

  Directory containing MAGs to annotate.

**`-x`**, **`--fasta-extension`** EXT

  File extension of FASTA files in `--genome-fasta-directory`. [default: fna]

**`-a`**, **`--assembly`** FILE [FILE ...]

  FASTA file(s) containing scaffolded contigs to pass to QUAST for QC.

## Annotation / bin processing options

**`--gtdb-path`** PATH

  Path to local GTDB database files.

**`--eggnog-db-path`** PATH

  Path to local EggNOG database files.

**`--singlem-metapackage-path`** PATH

  Path to local SingleM metapackage.

**`--checkm2-db-path`** PATH

  Path to CheckM2 database.

**`--metabuli-db-path`** PATH

  Path to local Metabuli database.

## Performance options

**`-t`**, **`--max-threads`** INT

  Maximum threads per process. [default: 8]

**`-n`**, **`--n-cores`** INT

  Maximum cores available. [default: 16]

**`-m`**, **`--max-memory`** INT

  Maximum memory in gigabytes. [default: 250]

**`-p`**, **`--pplacer-threads`** INT

  Threads for pplacer during GTDB-tk classification. Values above
  `--max-threads` are capped to that limit. [default: 8]

## Output options

**`-o`**, **`--output`** DIR

  Output directory. [default: ./]

**`--tmpdir`** DIR

  Temporary files directory.

## Inherited QC options

The parser also accepts `--gold-standard-assembly`, `--gsa-mappings`,
`--host-filter`, `--min-read-size`, `--min-mean-q`, `--keep-percent`,
`--min-short-read-length`, `--max-short-read-length`,
`--disable-adapter-trimming`, `--unqualified-percent-limit`,
`--quality-cutoff`, `--extra-fastp-params` and `--skip-qc` because `annotate`
inherits Aviary's QC option group. They do not affect the default `annotate`
target. Use `-a`/`--assembly` to provide assemblies to a separately targeted
QC rule, and see [`aviary assemble` → QC options](assemble.md#qc-options) for
the accepted values.
