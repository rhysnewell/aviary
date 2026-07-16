---
title: HPC & cluster submission
---

# HPC & cluster submission

Aviary and Snakemake are compatible with HPC cluster schedulers. Jobs can be submitted as either a single pipeline job or as individual rule-level jobs spread across nodes.

## Snakemake profile (recommended)

The preferred approach is to use a Snakemake profile via `--snakemake-profile`. Create a profile at `~/.config/snakemake/<profile-name>/config.yaml`:

```yaml
cluster: qsub
cluster-status: qstat
jobs: 10000
cluster-cancel: qdel
```

Then run Aviary with:

```
aviary recover -1 reads_1.fq.gz -2 reads_2.fq.gz --snakemake-profile cluster
```

Use `--cluster-retries` to automatically retry failed jobs with increasing time and memory:

```
aviary recover -1 reads_1.fq.gz -2 reads_2.fq.gz --snakemake-profile cluster --cluster-retries 3
```

CPU and memory bounds set via `--max-threads` and `--max-memory` are used as hard caps for submitted jobs.

Job resources were set based on [empirical data from 1,000 Aviary runs](https://github.com/rhysnewell/aviary/pull/270). The number of CPUs was set based on average mean load (to the nearest multiple of 8). Maximum memory was set based on the nearest power of 2, rounding up from maximum RSS.

## Using --snakemake-cmds

For simpler cluster setups, pass cluster options directly:

```
aviary assemble -1 reads_1.fq.gz -2 reads_2.fq.gz --longreads reads.fastq.gz \
    --long_read_type ont -t 24 -p 24 -n 24 --snakemake-cmds '--cluster qsub '
```

> **Note:** The trailing space after `qsub` is required due to a quirk in Python's argparse module.

## Running the coordinator job

When using cluster submission, Snakemake acts as a lightweight coordinator that submits individual rules as jobs. Submit the coordinator itself with minimal resources:

```
mqsub -m 8 -t 1 -w 48:00:00 --name aviary_recover -- \
    aviary recover -1 reads_1.fq.gz -2 reads_2.fq.gz --snakemake-profile aqua
```

## Local parallel execution

To run multiple rules in parallel locally, set `--n-cores` to a multiple of `--max-threads`:

```
aviary recover -1 reads_1.fq.gz -2 reads_2.fq.gz -t 8 -n 32
```

This allows up to 4 rules to run concurrently (32 cores / 8 threads each).
