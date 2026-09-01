---
title: CLI reference
---

# CLI reference

Aviary exposes eight subcommands. The installed command's `--help` output is
the final authority for the installed version.

Start with [command syntax](syntax.md) for global-option placement, boolean
values, list inputs, aliases and quoting rules.

| Command | Purpose |
| --- | --- |
| [`assemble`](assemble.md) | Quality-control reads and assemble contigs |
| [`recover`](recover.md) | Recover and assess MAGs from an assembly |
| [`annotate`](annotate.md) | Annotate a directory of genome FASTA files |
| [`complete`](complete.md) | Run the applicable metagenome stages end to end |
| [`cluster`](cluster.md) | Dereplicate genomes across completed runs |
| [`isolate`](isolate.md) | Assemble reads from a cultured isolate |
| [`configure`](configure.md) | Store reference-data paths and request downloads |
| [`build`](build.md) | Build dependency environments; add `--gpu` for GPU environments |

Shared workflow and resource options are explained in
[shared options](centralised_commands.md).
