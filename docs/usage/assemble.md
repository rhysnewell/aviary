---
title: aviary assemble
---

# aviary assemble

Step-down hybrid assembly using long and short reads, or assembly using only short or long reads.

```
aviary assemble -1 reads_1.fq.gz -2 reads_2.fq.gz --longreads reads.fastq.gz --long_read_type ont
```

## Input options (short reads)

**`-1`**, **`--pe-1`** FILE [FILE ...]

  Forward read files. If multiple files are provided and no longreads are given, all samples will be co-assembled with megahit or metaspades.

**`-2`**, **`--pe-2`** FILE [FILE ...]

  Reverse read files.

**`-i`**, **`--interleaved`** FILE [FILE ...]

  Interleaved read files.

**`-c`**, **`--coupled`** FILE [FILE ...]

  Forward and reverse read files in a coupled space-separated list.

## Input options (long reads)

**`-l`**, **`--longreads`** FILE [FILE ...]

  Long-read files. The first file will be used for assembly unless `--coassemble` is set.

**`-z`**, **`--longread-type`** TYPE

  Sequencing platform: `rs` (PacBio RSII), `sq` (PacBio Sequel), `ccs` (PacBio CCS), `hifi` (PacBio HiFi), `ont` (Oxford Nanopore), `ont_hq` (ONT high quality, Guppy5+ or Q20). [default: ont]

**`--long-read-assembler`** ASSEMBLER

  Long-read assembler to use. `myloasm` (default) or `flye`. [default: myloasm]

**`--medaka-model`** MODEL

  Medaka model for polishing long reads. [default: r941_min_hac_g507]

## Assembly options

**`--use-unicycler`**

  Use Unicycler to re-assemble the metaSPAdes hybrid assembly. Not recommended for complex metagenomes.

## QC options

**`-r`**, **`--host-filter`** FILE [FILE ...]

  Host reference FASTA files for removal of contaminant reads prior to assembly.

**`--skip-qc`**

  Skip quality control steps.

**`--min-read-size`** INT

  Minimum long read size when filtering using Filtlong. [default: 100]

**`--min-mean-q`** INT

  Minimum long read mean quality threshold. [default: 10]

**`--min-short-read-length`** INT

  Minimum length of short reads to keep. [default: 15]

**`--max-short-read-length`** INT

  Maximum length of short reads to keep, 0 = no maximum. [default: 0]

**`--disable-adapter-trimming`**

  Disable adapter trimming of short reads.

**`--quality-cutoff`** INT

  Phred quality value threshold for short reads. [default: 15]

**`--unqualified-percent-limit`** INT

  Percentage of bases allowed to be unqualified. [default: 40]

**`--extra-fastp-params`** STRING

  Extra parameters to pass to fastp, e.g. `--extra-fastp-params "-V -e 10"`.

## Performance options

**`-t`**, **`--max-threads`** INT

  Maximum threads given to any particular process. [default: 8]

**`-n`**, **`--n-cores`** INT

  Maximum cores available. Setting to multiples of `--max-threads` allows parallel processes. [default: 16]

**`-m`**, **`--max-memory`** INT

  Maximum memory in gigabytes. [default: 250]

**`-p`**, **`--pplacer-threads`** INT

  Threads given to pplacer. [default: 8]

**`--local-cores`** INT

  Maximum local cores when submitting to a cluster. [default: 16]

## Output options

**`-o`**, **`--output`** DIR

  Output directory. [default: ./]

**`--tmpdir`** DIR

  Temporary files directory. Uses `TMPDIR` environment variable if not set.

## Misc options

**`--snakemake-profile`** PROFILE

  Snakemake profile for cluster submission. See the Guides section for HPC usage.

**`--cluster-retries`** INT

  Number of retries for failed cluster jobs. [default: 0]

**`--dry-run`**

  Perform a snakemake dry run.

**`--clean`**

  Clean up temporary files after completion. [default: True]

**`--snakemake-cmds`** STRING

  Additional snakemake commands as a single string.

## Examples

Hybrid assembly from paired short reads and long reads:
```
aviary assemble -1 reads_1.fq.gz -2 reads_2.fq.gz --longreads reads.fastq.gz --long_read_type ont
```

Short-read-only assembly:
```
aviary assemble -1 reads_1.fq.gz -2 reads_2.fq.gz
```

Long-read-only assembly:
```
aviary assemble --longreads reads.fastq.gz --long_read_type ont
```
