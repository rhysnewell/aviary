---
title: Reproducibility
---

# Reproducibility

A reproducible Aviary analysis requires more than retaining the final MAGs.
Software environments, reference databases, the command line and workflow
metadata all affect the result.

## Record before running

```bash
aviary --version
conda list --explicit > aviary-environment.txt
```

Also record:

- full input paths and file checksums;
- the complete Aviary command;
- GTDB, EggNOG, CheckM2, SingleM and Metabuli releases used by enabled stages;
- the Snakemake profile and scheduler configuration;
- CPU, memory, GPU and temporary-storage settings;
- whether quality control, taxonomy, abundance or individual binners were
  skipped.

When installed from this repository, the checked-in `aviary/pixi.lock` pins the
dependency solution used by Pixi. Preserve the lockfile revision alongside the
Aviary source revision.

## Retain the run state

Keep the generated run configuration, `.snakemake/` metadata, `logs/`,
`benchmarks/`, final output links and their canonical targets in `data/`.
Copying only `bins/final_bins/` loses the evidence needed to explain how those
genomes were produced.

## Database sensitivity

Taxonomy, quality assessment and functional annotation depend on external
reference data. A rerun with newer databases is a new analysis and can produce
different classifications even when sequence inputs are unchanged.

## Rerun triggers

Aviary defaults `--rerun-triggers` to `mtime`. Other accepted triggers are
`params`, `input`, `software-env` and `code`. Changing trigger policy can change
which rules Snakemake considers out of date; record it whenever it differs from
the default.

## Publication archive

Archive the command, environment description, database releases, key logs,
summaries and checksums with the outputs used in analysis. Cite Aviary and each
upstream tool used by the selected workflow; see [citations](../citations.md).
