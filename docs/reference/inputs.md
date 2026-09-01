---
title: Input reference
---

# Input reference

## Short reads

Aviary accepts paired reads as separate forward and reverse lists (`-1` and
`-2`), interleaved reads (`--interleaved`), or a coupled list (`--coupled`).
Choose one of these layouts; `-1`, `--interleaved` and `--coupled` are
mutually exclusive.

```text
reads/
├── sample_A_R1.fastq.gz
├── sample_A_R2.fastq.gz
├── sample_B_R1.fastq.gz
└── sample_B_R2.fastq.gz
```

When supplying lists to `-1` and `-2`, keep samples in the same order. FASTQ
files may be gzip-compressed. With multiple files, assembly behaviour depends
on `--coassemble`; all supplied reads may still contribute differential
coverage during genome recovery.

=== "Paired lists"

    ```bash
    aviary complete \
      -1 sample_A_R1.fastq.gz sample_B_R1.fastq.gz \
      -2 sample_A_R2.fastq.gz sample_B_R2.fastq.gz
    ```

=== "Interleaved"

    ```bash
    aviary complete --interleaved sample_A.interleaved.fastq.gz
    ```

=== "Coupled list"

    ```bash
    aviary complete --coupled \
      sample_A_R1.fastq.gz sample_A_R2.fastq.gz \
      sample_B_R1.fastq.gz sample_B_R2.fastq.gz
    ```

!!! warning "Short-read identity option limitation"

    In the current parser, `--min-percent-read-identity-short` is in the same
    mutually exclusive group as `-1`, `--interleaved` and `--coupled`.
    Supplying it together with one of those inputs is rejected before the
    workflow starts. Its default remains `95`; an explicit override is only
    accepted with input shapes that do not use those mutually exclusive flags.

## Long reads

Supply one or more files with `--longreads`. `--long-read-type` accepts:

| Value | Sequencing data |
| --- | --- |
| `ont` | Oxford Nanopore reads |
| `ont_hq` | high-quality/Q20 Oxford Nanopore reads |
| `rs` | PacBio RS II reads |
| `sq` | PacBio Sequel reads |
| `ccs` | PacBio CCS reads |
| `hifi` | PacBio HiFi reads |

The default is `ont`. This selection affects parameters passed to upstream
tools, so set it to the actual sequencing technology.

```bash
aviary assemble --longreads sample.nanopore.fastq.gz --long-read-type ont
```

## Assemblies

`--assembly` accepts FASTA files containing scaffolded metagenome contigs.
`recover` uses the supplied reads to calculate coverage. SemiBin2 multi-sample
mode accepts multiple assemblies and requires `--semibin-mode multi`.

```bash
aviary recover --assembly week1.fasta week2.fasta week3.fasta \
  -1 week1_R1.fq.gz week2_R1.fq.gz week3_R1.fq.gz \
  -2 week1_R2.fq.gz week2_R2.fq.gz week3_R2.fq.gz \
  --semibin-mode multi
```

## Genome collections

`annotate` reads FASTA files from `--genome-fasta-directory`. Use
`--fasta-extension` when files do not use the default `fna` extension.

```bash
aviary annotate --genome-fasta-directory external_bins/ --fasta-extension fa
```

## Previous Aviary runs

`cluster --input-runs` expects completed Aviary output directories containing
the final-bin collection and bin summary needed for dereplication.

```bash
aviary cluster --input-runs runs/sample_A runs/sample_B runs/sample_C
```
