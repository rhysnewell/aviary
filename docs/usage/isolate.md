---
title: aviary isolate
---

# aviary isolate

Step-down hybrid assembly for isolated pure culture sequencing results. For use with isolate (not metagenomic) sequencing data.

```
aviary isolate -1 reads_1.fq.gz -2 reads_2.fq.gz --longreads reads.fastq.gz --long_read_type ont
```

> This subcommand also accepts `--dry-run`, `--clean`, `--strict`, `--request-gpu`, `--build`, `--build-gpu`, `--download`, `--rerun-triggers`, `--default-resources`, `--snakemake-profile`, `--snakemake-cmds`, `--cluster-retries`, `--local-cores`, and `--workflow`, which are shared across every aviary subcommand — see [Centralised commands](centralised_commands.md).

## Input options (short reads)

**`-1`**, **`--pe-1`** FILE [FILE ...]

  Forward short read files.

**`-2`**, **`--pe-2`** FILE [FILE ...]

  Reverse short read files.

**`-i`**, **`--interleaved`** FILE [FILE ...]

  Interleaved read files.

**`-c`**, **`--coupled`** FILE [FILE ...]

  Forward and reverse read files in a coupled space-separated list.

## Input options (long reads)

**`-l`**, **`--longreads`** FILE [FILE ...]

  Long-read files.

**`-z`**, **`--longread-type`** TYPE

  Sequencing platform: `rs`, `sq`, `ccs`, `hifi`, `ont`, `ont_hq`. [default: ont]

**`--medaka-model`** MODEL

  Medaka model for polishing. [default: r941_min_hac_g507]

## Isolate options

**`--genome-size`** INT

  Approximate size of the isolate genome in base pairs. [default: 5000000]

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

Hybrid isolate assembly:
```
aviary isolate -1 reads_1.fq.gz -2 reads_2.fq.gz --longreads reads.fastq.gz --long_read_type ont
```

Short-read-only isolate assembly:
```
aviary isolate -1 reads_1.fq.gz -2 reads_2.fq.gz
```
