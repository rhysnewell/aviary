---
title: aviary isolate
---

# aviary isolate

Step-down hybrid assembly for isolated pure culture sequencing results. For use with isolate (not metagenomic) sequencing data.

```bash
aviary isolate -1 reads_1.fq.gz -2 reads_2.fq.gz --longreads reads.fastq.gz --long-read-type ont
```

> This subcommand also accepts the common workflow, resource, output and
> execution options described under [Shared options](centralised_commands.md).

## Examples

### Long-read only

The minimum viable input — Flye assembly, racon and medaka polishing, then dnaapler
reorientation. No illumina polishing round runs, since no short reads are given:
```bash
aviary isolate --longreads reads.fastq.gz --long-read-type ont
```

### Hybrid: long reads + short reads for polishing

The typical case for a bacterial isolate closed with ONT and cleaned up with Illumina — adds one
extra Pilon/racon polishing round using the short reads on top of the long-read assembly:
```bash
aviary isolate -1 reads_1.fq.gz -2 reads_2.fq.gz --longreads reads.fastq.gz --long-read-type ont
```

Short reads accept the same input shapes as the other subcommands:
```bash
aviary isolate -i sample_interleaved.fq.gz --longreads reads.fastq.gz --long-read-type ont
aviary isolate -c sample_1.fq.gz sample_2.fq.gz --longreads reads.fastq.gz --long-read-type ont
```

`aviary isolate` has no `--long-read-assembler`/`--use-megahit`/`--use-unicycler` flags — unlike
`assemble`/`recover`/`complete`, the isolate assembly path is fixed to Flye rather than
user-selectable, since it targets a single pure-culture genome rather than a mixed community.

### Sizing the assembly

`--genome-size` is present in the CLI but is not currently consumed by the
workflow, so changing it does not alter the assembly. The following is accepted
for compatibility but behaves like the default invocation:
```bash
aviary isolate --longreads reads.fastq.gz --long-read-type ont --genome-size 4500000
```

## Input options (short reads)

Long reads (below) are required — isolate assembly is Flye-based and has no short-read-only
path. Short reads are optional and, when given, are used for one additional Pilon/racon
polishing round on top of the long-read assembly.

**`-1`**, **`--pe-1`** FILE [FILE ...]

  Forward short read files.

**`-2`**, **`--pe-2`** FILE [FILE ...]

  Reverse short read files.

**`-i`**, **`--interleaved`** FILE [FILE ...]

  Interleaved read files.

**`-c`**, **`--coupled`** FILE [FILE ...]

  Forward and reverse read files in a coupled space-separated list.

## Input options (long reads)

**`-l`**, **`--longreads`** FILE [FILE ...]

  Long-read files.

**`-z`**, **`--longread-type`** TYPE

  Sequencing platform: `rs`, `sq`, `ccs`, `hifi`, `ont`, `ont_hq`. [default: ont]

**`--medaka-model`** MODEL

  Model passed to the shared polishing helper. [default: r941_min_hac_g507]

> Reads supplied above are mapped back onto the assembly for coverage and
> polishing. See [Read mappers](mappers.md) for the
> `--short-read-mapper`/`--long-read-mapper` family of flags that control
> which aligner is used.

## Isolate options

**`--guppy-model`** MODEL

  Medaka model used by the isolate-specific `polish_isolate_medaka` rule.
  Despite the historical option name, this is a Medaka model identifier.
  The default is passed as a scalar; the current parser stores an explicit
  override as a one-item list, so custom values should be treated as a known
  compatibility limitation. [default: r941_min_hac_g507]

**`--genome-size`** INT

  Accepted by the current CLI for compatibility, but not read by the current
  isolate workflow and therefore does not tune Flye. [default: 5000000]

## QC options

Short-read and long-read filtering uses the same controls as `assemble`:
`--host-filter`, `--gold-standard-assembly`, `--gsa-mappings`,
`--min-read-size`, `--min-mean-q`, `--keep-percent`,
`--min-short-read-length`, `--max-short-read-length`,
`--disable-adapter-trimming`, `--unqualified-percent-limit`,
`--quality-cutoff`, `--extra-fastp-params` and `--skip-qc`. See
[`aviary assemble` → QC options](assemble.md#qc-options) for types, defaults
and examples.

## Performance options

**`-t`**, **`--max-threads`** INT

  Maximum threads per process. [default: 8]

**`-n`**, **`--n-cores`** INT

  Maximum cores available. [default: 16]

**`-m`**, **`--max-memory`** INT

  Maximum memory in gigabytes. [default: 250]

**`-p`**, **`--pplacer-threads`** INT

  Accepted through the shared parser but not used by the default isolate
  workflow. [default: 8]

## Output options

**`-o`**, **`--output`** DIR

  Output directory. [default: ./]

**`--tmpdir`** DIR

  Temporary files directory.

## Inherited compatibility options

`isolate` also inherits the binning and database-path groups used by the
metagenome workflows. The default `dnaapler` isolate target does not use them:

- `--min-contig-size`, `--min-bin-size`, `--coverage-job-strategy`,
  `--coverage-samples-per-job`, `--semibin-model`, `--semibin-mode`,
  `--refinery-max-iterations`, `--refinery-max-retries`, `--extra-binners`,
  `--skip-binners`, `--binning-only`, `--skip-abundances`, `--skip-taxonomy`,
  `--skip-singlem`, `--min-completeness` and `--max-contamination`
- `--min-percent-read-identity-short`, `--min-percent-read-identity-long`
- `--gtdb-path`, `--eggnog-db-path`, `--singlem-metapackage-path`,
  `--checkm2-db-path`, `--metabuli-db-path`

See [`aviary recover`](recover.md#binning-options) for the binning values and
[`aviary annotate`](annotate.md#annotation-bin-processing-options) for database
paths.
