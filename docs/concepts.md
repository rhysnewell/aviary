---
title: Core concepts
---

# Core concepts

Aviary is a workflow orchestrator for metagenomic assembly and genome recovery.
Understanding four distinctions makes its commands and outputs easier to use.

## Workflow commands describe stopping points

`assemble`, `recover` and `annotate` are composable stages. `complete` runs the
full path permitted by the supplied inputs, while `isolate` uses an assembly
path intended for a cultured isolate rather than a mixed community.

| Command | Begins with | Principal result |
| --- | --- | --- |
| `assemble` | sequencing reads | assembled contigs |
| `recover` | assembly plus reads used for coverage | recovered MAGs |
| `annotate` | MAG FASTA files | taxonomy and functional annotations |
| `complete` | reads or an existing assembly | results through annotation |
| `cluster` | completed Aviary runs | dereplicated representative genomes |
| `isolate` | reads from a cultured isolate | isolate assembly |

## Assemblies, bins and MAGs

An assembly joins overlapping sequence evidence into contigs. Metagenomic
binning then groups contigs that appear to originate from the same population.
A final bin is treated as a metagenome-assembled genome (MAG), but its
completeness, contamination and taxonomy remain estimates. Use
`bins/bin_info.tsv` to assess each recovered genome rather than treating all
FASTA files as equivalent.

## Aviary and upstream tools

Aviary owns workflow decisions: it validates command-line input, writes the run
configuration, selects Snakemake targets, passes resource limits and links
final products into stable output locations. Upstream programs perform the
underlying read filtering, assembly, mapping, binning, quality estimation,
taxonomic classification and functional annotation. The exact tools invoked
depend on command options and available data.

## Requested resources have different scopes

`--max-threads` is the maximum made available to one tool. `--n-cores` is the total
capacity Snakemake may schedule concurrently. `--local-cores` limits work kept
on the coordinator node when using a cluster profile. `--max-memory` is a hard
workflow cap in gigabytes, not a guarantee that every tool uses that amount.

## Workflow state and resumability

Snakemake determines whether an output is current from its inputs, rules and
metadata. A repeated Aviary command against the same output directory normally
continues incomplete work. Options such as `--rerun-triggers`, `--clean` and
`--unlock` deliberately change this behaviour; use them only after reading the
[workflow-control guide](guides/workflow-control.md).

## Inputs and outputs

Aviary accepts FASTQ sequencing reads and FASTA assemblies or genomes. BAM
files are produced internally when reads are mapped to assemblies. Final
results are presented through stable directories such as `assembly/`, `bins/`,
`taxonomy/` and `annotation/`; working files remain in `data/`.

Continue with the [assembly guide](guide/assembly.md),
[genome-recovery guide](guide/genome-recovery.md), or
[annotation guide](guide/annotation.md).
