---
title: Genome annotation
---

# Genome annotation

Annotation adds taxonomic and functional context to recovered genomes. Aviary
prepares the MAG inputs and coordinates external annotation tools; the
reference databases and upstream software perform the classifications.

## Basic usage

```bash
aviary annotate \
  --genome-fasta-directory recovery_run/bins/final_bins \
  --output annotation_run \
  --max-threads 16 \
  --n-cores 16
```

Use the literal option names shown by `aviary annotate --help` for the installed
version; the [CLI reference](../usage/annotate.md) documents aliases and
directory/extension controls.

## Reference data matters

Taxonomic assignments depend on the configured GTDB data. Functional
annotations depend on the configured EggNOG data. Record database releases
alongside the Aviary version, and do not compare annotations produced from
different releases as if the reference context were identical.

## Outputs

Annotation products are written under `annotation/`, while raw taxonomic output
may also appear under `taxonomy/` depending on the selected targets. Consult
the [output reference](../guides/output.md) and cite the upstream tools used by
the run using the [citation guide](../citations.md).
