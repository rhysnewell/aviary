---
title: Quickstart
---

# Quickstart

This example runs the complete short-read metagenome workflow. It assumes
Aviary is installed and the required databases are configured.

## 1. Check the installation

```bash
aviary --version
aviary complete --full-help
```

## 2. Prepare paired reads

Use one forward and one reverse FASTQ file for the same sample:

```text
reads/
├── sample_R1.fastq.gz
└── sample_R2.fastq.gz
```

Compressed FASTQ input is accepted. Forward and reverse files must be supplied
in matching order when more than one pair is given.

## 3. Preview the workflow

```bash
aviary complete \
  -1 reads/sample_R1.fastq.gz \
  -2 reads/sample_R2.fastq.gz \
  --output sample_aviary \
  --max-threads 8 \
  --n-cores 8 \
  --dry-run
```

The dry run resolves the workflow without executing analysis tools. Remove
`--dry-run` after checking the planned jobs and configured database paths.

## 4. Run

```bash
aviary complete \
  -1 reads/sample_R1.fastq.gz \
  -2 reads/sample_R2.fastq.gz \
  --output sample_aviary \
  --max-threads 8 \
  --n-cores 8
```

`--max-threads` limits an individual tool. `--n-cores` is the total CPU capacity
available to Snakemake, which may schedule several compatible jobs at once.

## 5. Inspect the results

Start with:

```text
sample_aviary/
├── assembly/final_contigs.fasta
├── bins/bin_info.tsv
├── bins/final_bins/
├── benchmarks/
└── logs/
```

Read [`bin_info.tsv` and the rest of the output layout](../guides/output.md)
before using recovered genomes downstream. If a rule fails, find its log in
`sample_aviary/logs/` and see [troubleshooting](../faqs.md).

## Next step

The [first complete analysis](first-analysis.md) explains how the
stages behave and which controls are most important.
