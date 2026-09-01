---
title: Genome recovery
---

# Genome recovery

Genome recovery groups assembled contigs into candidate genomes, refines the
candidate sets and evaluates the resulting MAGs. Reads are still important:
mapping them back to the assembly provides coverage information used by
binning and abundance calculations.

## Basic usage

```bash
aviary recover \
  --assembly assembly/final_contigs.fasta \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  --output recovery_run \
  --max-threads 16 \
  --n-cores 32
```

## How the stage behaves

Aviary prepares coverage inputs, invokes the enabled binners, refines candidate
bins and combines evidence into a final set. Quality assessment and taxonomy
then add evidence needed to interpret each MAG. A binning algorithm producing
a FASTA file does not by itself establish that the genome is complete,
uncontaminated or correctly classified.

## Important controls

`--min-bin-size` excludes candidate genomes smaller than the configured number
of bases from later processing. Binner selection options change which upstream
algorithms contribute candidates. GPU variants require a compatible GPU and
the corresponding Aviary environments; enabling a GPU flag does not make every
stage GPU accelerated.

SemiBin2 multi-sample mode accepts multiple assemblies and should be used only
with inputs organised for that mode. See the exact validation and accepted
values in the [`recover` reference](../usage/recover.md).

## Interpreting results

Use `bins/final_bins/` as the recovered FASTA collection and
`bins/bin_info.tsv` as its main summary. Retain the run's database versions and
software environment because taxonomy and quality estimates can change with
reference data and dependency versions.

Continue with [annotation](annotation.md), the
[output reference](../guides/output.md), or
[reproducibility guidance](../advanced/reproducibility.md).
