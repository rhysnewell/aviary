---
title: aviary assemble
---

# aviary assemble

Step-down hybrid assembly using long and short reads, or assembly using only short or long reads.

```
aviary assemble -1 reads_1.fq.gz -2 reads_2.fq.gz --longreads reads.fastq.gz --long_read_type ont
```

> This subcommand also accepts `--build`, `--build-gpu`, `--download`, `--rerun-triggers`, `--default-resources`, and `--workflow`, which are shared across every aviary subcommand — see [Centralised commands](centralised_commands.md).

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

**`--use-megahit`**

  Use MEGAHIT instead of metaSPAdes for short-read-only assembly. [default: false]

**`--coassemble`**, **`--co-assemble`**

  When multiple read sets are given, coassemble them together. If false, aviary uses only the first short-read and first long-read set for assembly (all read sets are still used for differential-coverage binning). [default: false]

**`-k`**, **`--kmer-sizes`** INT [INT ...]

  Manually specify the k-mer sizes used by SPAdes during assembly. Space-separated odd integers less than 128, or `auto`. [default: auto]

**`--min-cov-long`** INT

  Automatically include Flye contigs with long-read coverage ≥ this value. High long-read coverage indicates the overlap-layout-consensus assembly is more likely to be correct. [default: 5]

**`--min-cov-short`** INT

  Automatically include Flye contigs with short-read coverage ≤ this value. Low short-read coverage indicates metaSPAdes would not assemble this contig better. [default: 5]

**`--exclude-contig-cov`** INT

  Automatically exclude Flye contigs with long-read coverage ≤ this value, provided their length is also ≤ `--exclude-contig-size`. [default: 10]

**`--exclude-contig-size`** INT

  Automatically exclude Flye contigs with length ≤ this value, provided their long-read coverage is also ≤ `--exclude-contig-cov`. [default: 2500]

**`--include-contig-size`** INT

  Automatically include Flye contigs with length ≥ this value. [default: 10000]

## QC options

**`-r`**, **`--host-filter`** FILE [FILE ...]

  Host reference FASTA files for removal of contaminant reads prior to assembly.

**`-g`**, **`--gold-standard-assembly`** FILE [FILE ...]

  A gold-standard assembly to compare the resulting (or a given input) assembly against.

**`--gsa-mappings`** FILE

  CAMI I & II gold-standard-assembly mappings, used alongside `--gold-standard-assembly`.

**`--keep-percent`** INT

  **Deprecated.** Percentage of reads passing quality thresholds kept by Filtlong. [default: 100]

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
