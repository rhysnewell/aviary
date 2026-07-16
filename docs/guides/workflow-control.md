---
title: Workflow control
---

# Workflow control

## Targeting specific rules

Aviary uses Snakemake under the hood, which means you can target specific rules rather than running a full module. Use `-w/--workflow` to specify the rule to run up to:

```
aviary recover -w rosella --assembly scaffolds.fasta \
    -1 reads_1.fq.gz -2 reads_2.fq.gz --output output_dir/ --max-threads 12
```

All steps leading up to the targeted rule will still run if they haven't been completed yet.

The available rules for each module are defined in the module's Snakefile under `aviary/modules/`.

## Resuming interrupted runs

If an Aviary workflow is interrupted, re-running the same command will resume from the last completed step — Snakemake tracks completed outputs automatically.

## Dry run

Test the workflow order and verify conda environments without executing any rules:

```
aviary recover -1 reads_1.fq.gz -2 reads_2.fq.gz --dry-run
```

## Skipping steps

Several steps can be skipped to speed up runs when not needed:

| Flag | Effect |
| --- | --- |
| `--skip-qc` | Skip read quality control |
| `--skip-taxonomy` | Skip GTDB-tk taxonomy assignment |
| `--skip-abundances` | Skip CoverM abundance calculations |
| `--skip-singlem` | Skip SingleM recovery assessment |
| `--binning-only` | Stop after binning, skip all downstream steps |

## Temporary file cleanup

By default Aviary removes temporary files (BAM files, intermediate FASTQs) on completion. To keep them for partial reruns:

```
aviary recover -1 reads_1.fq.gz -2 reads_2.fq.gz --clean false
```

> **Note:** Keeping temporary files avoids re-generating them on reruns but will use significantly more disk space.
