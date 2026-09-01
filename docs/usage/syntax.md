---
title: Command syntax
---

# Command syntax

Every Aviary invocation has a command followed by that command's options:

```text
aviary [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]
```

```bash
aviary --verbosity 5 recover --assembly assembly.fasta \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz --output results/
```

Global options must come **before** the command. Command options come after it.

## Global options

These three options belong to the top-level `aviary` parser and are not valid
after a subcommand.

| Option | Value | Meaning |
| --- | --- | --- |
| `--version` | none | Print the installed Aviary version and exit. |
| `--verbosity` | `1`–`5` | Set logging to critical, error, warning, info or debug. [default: `4`] |
| `--log` | `FILE` | Write logging to a new file. Aviary refuses to overwrite an existing file. |

```bash
aviary --version
aviary --verbosity 5 assemble -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz
aviary --log aviary.log complete -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz
```

## Help

```bash
aviary --help
aviary recover --help
aviary recover --full-help
```

For `assemble`, `recover` and `complete`, `--help` shows a concise command
summary and examples; use `--full-help` (or `--full_help`) for the complete
option list. The other commands show their complete option list with
`--help`.

## Files and directories

Options shown as `FILE [FILE ...]` or `DIR [DIR ...]` accept a space-separated
list. Shell globs are also valid when they expand in the intended sample order:

```bash
aviary recover --assembly assembly.fasta \
  -1 sample_A_R1.fastq.gz sample_B_R1.fastq.gz \
  -2 sample_A_R2.fastq.gz sample_B_R2.fastq.gz

aviary cluster --input-runs runs/sample_A runs/sample_B runs/sample_C
```

Keep paired forward and reverse lists in matching order. Quote paths that
contain spaces. A value beginning with `-` must usually be quoted as part of a
single string, as with `--snakemake-cmds` and `--pggb-params`.

## Boolean values

Boolean options accept either a bare flag or an explicit value:

```bash
aviary recover --assembly assembly.fasta -1 R1.fq.gz -2 R2.fq.gz --strict
aviary recover --assembly assembly.fasta -1 R1.fq.gz -2 R2.fq.gz --strict true
aviary recover --assembly assembly.fasta -1 R1.fq.gz -2 R2.fq.gz --clean false
```

A bare flag means `true`. Accepted true values are `yes`, `true`, `t`, `y` and
`1`; false values are `no`, `false`, `f`, `n` and `0`, case-insensitively.
There are no `--no-clean`-style negated flags. This syntax applies to
`--request-gpu`, `--strict`, `--dry-run`, `--clean`, `--build`, `--build-gpu`,
`--disable-adapter-trimming`, `--skip-qc`, `--use-unicycler`, `--use-megahit`,
`--coassemble`, `--binning-only`, `--skip-abundances`, `--skip-taxonomy`,
`--skip-singlem`, `--use-checkm2-scores` and `build --gpu`.

## Aliases

The reference uses the canonical hyphenated spelling. Aviary also accepts
underscore aliases for most multiword options, for example
`--max-threads`/`--max_threads` and `--long-read-type`/`--long_read_type`.
Prefer hyphens in scripts and publications because that is the spelling shown
in this manual and in normal help output.

Common non-mechanical aliases are:

| Canonical form | Also accepted |
| --- | --- |
| `--pe-1`, `--pe-2` | `-1`, `-2`, `--paired-reads-1`, `--paired-reads-2`, `--pe1`, `--pe2` |
| `--longreads` | `-l`, `--long-reads`, `--long_reads` |
| `--longread-type` | `-z`, `--long-read-type`, `--long_read_type`, `--longread_type` |
| `--tmpdir` | `--tempdir`, `--tmp-dir`, `--tmp`, `--temp`, `--temp-dir` and underscore forms |
| `--dry-run` | `--dryrun`, `--dry_run` |
| `--coassemble` | `--co-assemble`, `--co_assemble` |
| `--extra-binners` | singular `--extra-binner` and underscore forms |
| `--skip-binners` | singular `--skip-binner` and underscore forms |

## Passing lists and option strings

List-valued options consume values until the next recognised Aviary option.
Put them before the next flag, and do not comma-separate them:

```bash
aviary recover --assembly assembly.fasta \
  --extra-binners concoct comebin \
  --skip-binners vamb metabat
```

Options that forward a command fragment take one quoted string:

```bash
aviary recover --assembly assembly.fasta -1 R1.fq.gz -2 R2.fq.gz \
  --snakemake-cmds "--printshellcmds --keep-going"
```

See [Input files](../reference/inputs.md) for all supported read layouts and
[Shared options](centralised_commands.md) for workflow, resource and output
controls.
