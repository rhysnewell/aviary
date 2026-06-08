---
title: Output files
---

# Output files

All Aviary subcommands write to the directory specified by `-o/--output` (default: `./`).

Rule benchmarks are in `benchmarks/`. Logs and error messages are in `logs/`. The `data/` folder contains intermediate and final outputs, with final results symlinked into their respective top-level folders.

## assemble

| Path | Description |
| --- | --- |
| `assembly/final_contigs.fasta` | Assembled contigs in FASTA format |
| `www/assembly_stats.txt` | Summary assembly statistics |
| `www/fastqc/` | FastQC output |
| `www/nanoplot/` | NanoPlot output (long reads) |

## recover

If assembly is performed during `recover`, the assemble outputs will also be present.

| Path | Description |
| --- | --- |
| `bins/bin_info.tsv` | Per-MAG summary: taxonomy, QC, size, N50, abundance |
| `bins/checkm_minimal.tsv` | Minimal CheckM2 summary |
| `bins/coverm_abundances.tsv` | MAG abundance per sample |
| `bins/final_bins/` | Final MAGs in FASTA format |
| `diversity/` | SingleM output |
| `taxonomy/` | GTDB-Tk output |

## annotate

| Path | Description |
| --- | --- |
| `annotation/` | EggNOG-mapper output |

## complete

Produces all outputs from assemble, recover, and annotate.
