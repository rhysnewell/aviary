---
title: Usage
---

# Usage

Aviary provides several subcommands for different stages of the metagenomics workflow:

| Subcommand | Description |
| --- | --- |
| `assemble` | Step-down hybrid assembly using long and short reads, or assembly using only short or long reads |
| `recover` | Recover MAGs from an assembly using multiple binning algorithms |
| `annotate` | Annotate a set of MAGs using EggNOG, GTDB-tk, and CheckM2 |
| `complete` | Run the full pipeline: assembly → binning → refinement → annotation |
| `cluster` | Dereplicate/choose representative genomes from multiple aviary runs |
| `isolate` | Hybrid isolate assembly for pure sequencing results |
| `configure` | Set or reset environment variables used by aviary |
