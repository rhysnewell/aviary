---
title: aviary cluster
---

# aviary cluster

Dereplicate and choose representative genomes from multiple aviary runs using Galah.

```
aviary cluster --input-runs aviary_output_folder_1/ aviary_output_folder_2/
```

> This subcommand also accepts `--dry-run`, `--clean`, `--strict`, `--request-gpu`, `--build`, `--build-gpu`, `--download`, `--rerun-triggers`, `--default-resources`, `--snakemake-profile`, `--snakemake-cmds`, `--cluster-retries`, `--local-cores`, and `--workflow`, which are shared across every aviary subcommand — see [Centralised commands](centralised_commands.md).

## Input options

**`-i`**, **`--input-runs`** DIR [DIR ...]

  Paths to previous finished aviary runs. Each must contain `bins/checkm.out` and `bins/final_bins`. **(required)**

## Clustering options

**`--ani`** FLOAT

  Overall ANI level to dereplicate at with FastANI. [default: 97]

**`--precluster-ani`** FLOAT

  Minimum dashing-derived ANI for preclustering. [default: 95]

**`--precluster-method`** METHOD

  Method for rough ANI in preclustering: `dashing` (HyperLogLog) or `finch` (MinHash). [default: dashing]

**`--min-completeness`** FLOAT

  Ignore genomes below this completeness percentage. [default: none]

**`--max-contamination`** FLOAT

  Ignore genomes above this contamination percentage. [default: none]

**`--use-checkm2-scores`**

  Use CheckM2 completeness and contamination scores for Galah dereplication.

**`--pggb-params`** STRING

  Parameters for pggb. [default: `-k 79 -G 7919,8069`]

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

## Examples

Cluster two aviary runs:
```
aviary cluster --input-runs run1/ run2/
```

Cluster with custom ANI threshold:
```
aviary cluster --input-runs run1/ run2/ --ani 95
```
