---
title: aviary annotate
---

# aviary annotate

Annotate a given set of MAGs using EggNOG, GTDB-tk, and CheckM2.

```
aviary annotate --genome-fasta-directory input_bins/
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

## Output options

**`-o`**, **`--output`** DIR

  Output directory. [default: ./]

**`--tmpdir`** DIR

  Temporary files directory.

## Examples

Annotate a directory of MAGs:
```
aviary annotate --genome-fasta-directory input_bins/
```

Annotate with a specific GTDB path:
```
aviary annotate --genome-fasta-directory input_bins/ --gtdb-path /path/to/gtdb/
```
