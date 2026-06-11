---
title: Configuration
---

# Configuration

## Environment variables

On first run, Aviary will prompt for database locations if they haven't been set. Use `aviary configure` to set or update them at any time:

```
aviary configure --gtdb-path /path/to/gtdb/ --checkm2-db-path /path/to/checkm2/
```

You can also set them manually in your `.bashrc`:

```bash
export GTDBTK_DATA_PATH=/path/to/gtdb/gtdb_release232/db/
export EGGNOG_DATA_DIR=/path/to/eggnog-mapper/2.1.3/
export SINGLEM_METAPACKAGE_PATH=/path/to/singlem/S6.5.0.GTDB_r232.metapackage_20260319.smpkg.zb
export CHECKM2DB=/path/to/checkm2db/uniref100.KO.1.dmnd
```

See `aviary configure --help` for all available options.

## Thread control

Aviary has three thread control parameters:

**`-t`**, **`--max-threads`**

Controls how many threads any single program can use. If set to 24, each tool with threading support will use up to 24 threads.

**`-n`**, **`--n-cores`**

Controls total cores given to Snakemake. Setting this higher than `--max-threads` allows multiple rules to run concurrently. Defaults to match `--max-threads` if not set.

**`-p`**, **`--pplacer-threads`**

Pplacer gets its own parameter because it can deadlock with too many threads and is memory-intensive at high thread counts.

Example — run up to 4 rules in parallel:
```
aviary recover -1 reads_1.fq.gz -2 reads_2.fq.gz -t 8 -n 32
```

## RAM control

Set maximum memory with `-m/--max-memory` (in gigabytes):

```
aviary recover -1 reads_1.fq.gz -2 reads_2.fq.gz -m 500
```

When using HPC cluster submission, requested job memory increases with each retry and is capped at `--max-memory`.

> **Note:** GTDB-Tk v2.7.0+ requires at least 140 GB of RAM for R232.

## Temporary directory

By default Aviary uses `/tmp`. To change it for a single run:

```
aviary recover -1 reads_1.fq.gz -2 reads_2.fq.gz --tmpdir /scratch/tmp/
```

To set it permanently:

```
aviary configure --tmpdir /scratch/tmp/
```
