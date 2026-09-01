---
title: Metagenome assembly
---

# Metagenome assembly

Assembly transforms quality-controlled reads into contigs for genome recovery
or other downstream analysis. Aviary supports short-read, long-read and hybrid
inputs and chooses applicable workflow rules from the data supplied.

## Select an input strategy

- Paired short reads require matching `-1` and `-2` lists.
- Long reads require `--longreads` and a `--long-read-type` value.
- Hybrid assembly supplies both short and long reads.
- Multiple samples can be co-assembled where the selected options support it.

Input quality control happens before assembly unless it is explicitly skipped.
Host references supplied through `--host-filter` are used to remove matching
reads before the assembly stage.

## Basic usage

```bash
aviary assemble \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  --output assembly_run \
  --max-threads 16 \
  --n-cores 16
```

For a hybrid dataset, add long reads and their platform type:

```bash
aviary assemble \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  --longreads sample_ont.fastq.gz \
  --long-read-type ont \
  --output hybrid_run
```

## Important controls

`--long-read-assembler` selects the supported long-read assembly implementation.
Quality-control thresholds affect which reads reach the assembler and should be
chosen for the sequencing data, not copied mechanically between projects.
`--min-contig-size` is accepted by `assemble` through its shared parser, but it
is a downstream binning filter and does not alter the default assembly target;
set it on `recover` or `complete` when filtering contigs for binning.

## Result

The canonical assembly is stored in `data/final_contigs.fasta` and exposed as
`assembly/final_contigs.fasta`. Reports are written under `www/`; rule logs and
benchmarks are kept in `logs/` and `benchmarks/`.

See the [`assemble` reference](../usage/assemble.md) for every option and the
[output reference](../guides/output.md) for the run layout.
