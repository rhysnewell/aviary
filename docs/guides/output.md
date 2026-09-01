---
title: Output reference
---

# Output reference

Every command writes beneath `--output` (default `./`). Aviary keeps working
products in `data/` and exposes important final products through stable
top-level paths, often as symbolic links. The exact tree depends on the command,
input types and skipped stages.

## Common run directories

```text
output_dir/
├── assembly/       final assembly link
├── bins/           final MAGs and genome summaries
├── annotation/     functional annotation products
├── taxonomy/       raw taxonomic classification products
├── diversity/      SingleM products when enabled
├── www/            quality-control reports and summaries
├── benchmarks/     per-rule Snakemake benchmark records
├── logs/           per-rule standard output and error logs
└── data/           working data and canonical generated files
```

Do not assume every directory is present. For example, `assemble` does not
produce `bins/`, and `--skip-taxonomy` suppresses the corresponding taxonomy
stage.

## Assembly outputs

| Path | Format | Status | Meaning |
| --- | --- | --- | --- |
| `assembly/final_contigs.fasta` | FASTA | final link | Assembly intended for downstream use |
| `data/final_contigs.fasta` | FASTA | final, canonical | Generated assembly targeted by the stable link |
| `www/assembly_stats.txt` | text | final report | Assembly size statistics |
| `www/fastp.html` | HTML | report, when short-read QC runs | Interactive short-read filtering report |
| `www/rastqc/` | HTML/support files | report, when enabled | Short-read contamination/QC report |
| `www/rastqc_long/` | HTML/support files | report, when enabled | Long-read contamination/QC report |
| `www/nanoplot/` | HTML/support files | report, with long reads | Long-read quality summary |

## Genome recovery outputs

| Path | Format | Status | Meaning |
| --- | --- | --- | --- |
| `bins/final_bins/` | directory of FASTA | final | Recovered MAG sequences |
| `bins/bin_info.tsv` | TSV | final | Main per-MAG quality, taxonomy and assembly summary |
| `bins/checkm_minimal.tsv` | TSV | final | Compact CheckM2 quality summary |
| `bins/coverm_abundances.tsv` | TSV | final when abundance runs | CoverM abundance estimates by sample |
| `taxonomy/` | directory | final/supporting | Raw GTDB-Tk results |
| `diversity/metagenome.combined_otu_table.csv` | CSV | final when SingleM runs | Read-based community profile |
| `diversity/singlem_appraisal.tsv` | TSV | final when SingleM appraisal runs | Comparison of read and genome recovery evidence |
| `diversity/singlem_appraise.svg` | SVG | final figure | Visual SingleM appraisal summary |

Treat `bin_info.tsv` as a summary, not a substitute for retaining the raw
quality and taxonomy outputs. Column sets can depend on enabled stages; inspect
the header from the installed version before writing downstream parsers.

## Annotation outputs

Functional and taxonomic products are written below `annotation/` and
`taxonomy/` according to the selected workflow targets. Preserve their logs and
database release information when interpreting or publishing results.

## Logs and benchmarks

Files below `logs/` are the first place to investigate a failed rule. Files
below `benchmarks/` are Snakemake benchmark records for completed rules and can
support local resource planning. They are observations from the current input
and environment, not general performance guarantees.

## Intermediate files and cleanup

The `data/` directory contains both canonical generated products and working
files. Do not remove it wholesale. With `--clean` enabled (the default),
Snakemake removes files declared temporary after their consumers finish. Use
`--clean false` only when retaining those intermediates is worth the storage
cost or they are needed for a planned partial workflow.
