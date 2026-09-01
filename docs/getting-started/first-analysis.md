---
title: First complete analysis
---

# First complete analysis

`aviary complete` joins assembly, genome recovery and annotation into one
Snakemake run. Use it when you want Aviary to carry raw reads through to a
curated set of metagenome-assembled genomes (MAGs).

## Data flow

```text
paired short reads ─┐
                    ├─ quality control ─ assembly ─ read mapping
optional long reads ┘                         │
                                              ▼
                         binning ─ refinement ─ quality assessment
                                              │
                                              ▼
                               taxonomy and functional annotation
```

Aviary orchestrates these stages; specialised dependencies perform the
underlying assembly, binning, quality assessment and annotation algorithms.

## Choose the input mode

| Data available | Required input | Typical use |
| --- | --- | --- |
| Paired short reads | `-1`, `-2` | Short-read metagenome assembly and recovery |
| Long reads | `--longreads`, `--long-read-type` | Long-read assembly and recovery |
| Both | short- and long-read options | Hybrid assembly and recovery |
| Existing assembly | `--assembly` plus reads for coverage | Skip de novo assembly and recover MAGs |

See the [input reference](../reference/inputs.md) for pairing rules and accepted
long-read type identifiers.

## Run in a dedicated output directory

```bash
aviary complete \
  -1 reads/sample_R1.fastq.gz \
  -2 reads/sample_R2.fastq.gz \
  --longreads reads/sample_ont.fastq.gz \
  --long-read-type ont \
  --output sample_aviary \
  --max-threads 16 \
  --n-cores 32 \
  --max-memory 250
```

The default maximum memory value in the workflow configuration is 250 GB.
Set the option to the actual hard limit available to the run; it is not a
prediction of typical memory consumption.

## Monitor and resume

Rule-specific messages are written beneath `logs/`, while Snakemake benchmark
records are written beneath `benchmarks/`. Re-running the same command and
output directory normally resumes from existing valid outputs. Do not delete
`.snakemake/` or intermediate files while diagnosing an interrupted run.

For scheduler submission, resource caps and retries, see
[HPC and cluster submission](../guides/hpc.md). For precise command options, see
the [`complete` CLI reference](../usage/complete.md).
