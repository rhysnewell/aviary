---
title: Worked analyses
---

# Worked analyses

These examples show how Aviary stages connect. Replace the illustrative paths
with real input files; no example dataset is bundled with the repository.

## Hybrid metagenome from reads to MAGs

### Input

```text
reads/
├── sample_R1.fastq.gz
├── sample_R2.fastq.gz
└── sample_ont.fastq.gz
```

### Preview

```bash
aviary complete \
  -1 reads/sample_R1.fastq.gz \
  -2 reads/sample_R2.fastq.gz \
  --longreads reads/sample_ont.fastq.gz \
  --long-read-type ont \
  --output hybrid_complete \
  --max-threads 16 \
  --n-cores 32 \
  --dry-run
```

### Run

Remove `--dry-run` after checking the planned targets and database paths:

```bash
aviary complete \
  -1 reads/sample_R1.fastq.gz \
  -2 reads/sample_R2.fastq.gz \
  --longreads reads/sample_ont.fastq.gz \
  --long-read-type ont \
  --output hybrid_complete \
  --max-threads 16 \
  --n-cores 32 \
  --max-memory 250
```

### Interpret

Inspect `assembly/final_contigs.fasta` for the final assembly,
`bins/bin_info.tsv` for the MAG summary and `bins/final_bins/` for genome FASTA
files. Check `logs/` for rule execution details and `benchmarks/` before
adjusting scheduler requests.

## Separate assembly and recovery

Splitting stages is useful when you want to inspect the assembly before
committing to genome recovery.

```bash
aviary assemble \
  -1 reads/sample_R1.fastq.gz \
  -2 reads/sample_R2.fastq.gz \
  --output 01_assembly \
  --max-threads 16 \
  --n-cores 16
```

```bash
aviary recover \
  --assembly 01_assembly/assembly/final_contigs.fasta \
  -1 reads/sample_R1.fastq.gz \
  -2 reads/sample_R2.fastq.gz \
  --output 02_recovery \
  --max-threads 16 \
  --n-cores 32
```

The reads are supplied again because coverage across the assembly contributes
to genome recovery and abundance estimation.

## Annotate an existing genome collection

```bash
aviary annotate \
  --genome-fasta-directory 02_recovery/bins/final_bins \
  --fasta-extension fna \
  --output 03_annotation \
  --max-threads 16 \
  --n-cores 16
```

Keep the annotation run with the reference database releases used. See the
[annotation guide](guide/annotation.md) for interpretation and the
[CLI reference](usage/README.md) for advanced controls.
