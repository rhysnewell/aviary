---
title: aviary complete
---

# aviary complete

Performs all steps in the Aviary pipeline: Assembly → Binning → Refinement → Annotation.

```
aviary complete -1 reads_1.fq.gz -2 reads_2.fq.gz --longreads reads.fastq.gz --long_read_type ont
```

> This subcommand also accepts `--build`, `--build-gpu`, `--download`, `--rerun-triggers`, `--default-resources`, and `--workflow`, which are shared across every aviary subcommand — see [Centralised commands](centralised_commands.md).

## Input options

**`-a`**, **`--assembly`** FILE [FILE ...]

  One or more FASTA files containing scaffolded contigs of metagenome assemblies. If not provided, aviary will assemble first. Provide multiple assemblies for SemiBin2 multi-sample binning (requires `--semibin-mode multi`).

**`-1`**, **`--pe-1`** FILE [FILE ...]

  Forward short read files.

**`-2`**, **`--pe-2`** FILE [FILE ...]

  Reverse short read files.

**`-i`**, **`--interleaved`** FILE [FILE ...]

  Interleaved read files.

**`-c`**, **`--coupled`** FILE [FILE ...]

  Forward and reverse read files in a coupled space-separated list.

**`-l`**, **`--longreads`** FILE [FILE ...]

  Long-read files.

**`-z`**, **`--longread-type`** TYPE

  Sequencing platform: `rs`, `sq`, `ccs`, `hifi`, `ont`, `ont_hq`. [default: ont]

## Assembly options

**`--long-read-assembler`** ASSEMBLER

  Long-read assembler: `myloasm` or `flye`. [default: myloasm]

**`--medaka-model`** MODEL

  Medaka model for long read polishing. [default: r941_min_hac_g507]

**`--use-unicycler`**

  Use Unicycler to re-assemble the metaSPAdes hybrid assembly. Not recommended for complex metagenomes.

**`--use-megahit`**

  Use MEGAHIT instead of metaSPAdes for short-read-only assembly. [default: false]

**`--coassemble`**, **`--co-assemble`**

  When multiple read sets are given, coassemble them together. If false, aviary uses only the first short-read and first long-read set for assembly (all read sets are still used for differential-coverage binning). [default: false]

**`-k`**, **`--kmer-sizes`** INT [INT ...]

  Manually specify the k-mer sizes used by SPAdes during assembly. Space-separated odd integers less than 128, or `auto`. [default: auto]

**`--min-cov-long`** INT

  Automatically include Flye contigs with long-read coverage ≥ this value. [default: 5]

**`--min-cov-short`** INT

  Automatically include Flye contigs with short-read coverage ≤ this value. [default: 5]

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

## Binning options

**`-s`**, **`--min-contig-size`** INT

  Minimum contig size in base pairs for binning. [default: 1500]

**`-b`**, **`--min-bin-size`** INT

  Minimum bin size in base pairs for a MAG. [default: 200000]

**`--extra-binners`** BINNER [BINNER ...]

  Extra binning algorithms to run: `maxbin2`, `concoct`, `comebin`, `taxvamb`, `quickbin`. These are skipped by default due to long runtimes.

**`--skip-binners`** BINNER [BINNER ...]

  Binning algorithms to skip: `rosella`, `semibin`, `metabat1`, `metabat2`, `metabat`, `vamb`, `quickbin`.

**`--semibin-model`** MODEL

  SemiBin environment model: `human_gut`, `dog_gut`, `ocean`, `soil`, `cat_gut`, `human_oral`, `mouse_gut`, `pig_gut`, `built_environment`, `wastewater`, `global`. [default: global]

**`--semibin-mode`** `{single,multi}`

  SemiBin2 mode to use. `single` runs `single_easy_bin` on one assembly at a time. `multi` runs `multi_easy_bin`, co-binning multiple assemblies together — pass two or more files to `--assembly` to use it. Multi mode ignores `--semibin-model`, as pre-trained environments aren't supported for multi-sample binning. [default: single]

  ```
  aviary complete --assembly sample1.fasta sample2.fasta sample3.fasta \
    -1 sample1_1.fq.gz sample2_1.fq.gz sample3_1.fq.gz \
    -2 sample1_2.fq.gz sample2_2.fq.gz sample3_2.fq.gz \
    --semibin-mode multi
  ```

  Assemblies are concatenated (with unique `sample:contig` headers) into one SemiBin2 input, and reads from every sample are mapped back onto that concatenation, so co-abundance across samples improves binning. Output bins are written per-sample-prefixed into `data/semibin_bins/output_bins/`.

**`--refinery-max-iterations`** INT

  Maximum Rosella refinery iterations. Set to 0 to skip. [default: 5]

**`--refinery-max-retries`** INT

  Maximum Rosella refinery retries per iteration. [default: 3]

**`--binning-only`**

  Stop after binning. Skip SingleM, GTDB-tk, and CoverM.

**`--skip-abundances`**

  Skip CoverM post-binning abundance calculations.

**`--skip-taxonomy`**

  Skip GTDB-tk post-binning taxonomy assignment.

**`--skip-singlem`**

  Skip SingleM post-binning recovery assessment. [default: True]

**`--min-completeness`** FLOAT

  Minimum CheckM2 completeness percentage for annotation. [default: 50.0]

**`--max-contamination`** FLOAT

  Maximum CheckM2 contamination percentage for annotation. [default: 5.0]

**`--coverage-job-strategy`** STRATEGY

  Strategy for coverage calculation across many samples: `default`, `never`, `always`. [default: default]

**`--coverage-samples-per-job`** INT

  Number of samples per coverage job when splitting. [default: 5]

**`--min-percent-read-identity-short`** FLOAT

  Minimum percent read identity used by CoverM for short reads when calculating genome abundances. [default: 95]

**`--min-percent-read-identity-long`** FLOAT

  Minimum percent read identity used by CoverM for long reads when calculating genome abundances. [default: 85]

## Annotation / bin processing options

**`--gtdb-path`** PATH

  Path to local GTDB database files.

**`--eggnog-db-path`** PATH

  Path to local EggNOG database files.

**`--singlem-metapackage-path`** PATH

  Path to local SingleM metapackage.

**`--checkm2-db-path`** PATH

  Path to CheckM2 database.

**`--metabuli-db-path`** PATH

  Path to local Metabuli database.

## Performance options

**`-t`**, **`--max-threads`** INT

  Maximum threads per process. [default: 8]

**`-n`**, **`--n-cores`** INT

  Maximum cores available. [default: 16]

**`-m`**, **`--max-memory`** INT

  Maximum memory in gigabytes. [default: 250]

**`-p`**, **`--pplacer-threads`** INT

  Threads for pplacer. [default: 8]

**`--local-cores`** INT

  Maximum cores available locally. Only relevant when submitting to a cluster (see `--snakemake-profile`), in which case `--n-cores` restricts cores requested per submitted job. [default: 16]

## Output options

**`-o`**, **`--output`** DIR

  Output directory. [default: ./]

**`--tmpdir`** DIR

  Temporary files directory.

## Misc options

**`--snakemake-profile`** PROFILE

  Snakemake profile for cluster submission. See the Guides section for HPC usage.

**`--cluster-retries`** INT

  Retries for failed cluster jobs. [default: 0]

**`--dry-run`**

  Perform a snakemake dry run.

**`--clean`**

  Clean up temporary files. [default: True]

**`--strict`**

  Ensure each binner completes successfully. [default: skip failing binners]

**`--request-gpu`**

  Request a GPU for the pipeline (taxvamb, comebin, semibin). Only takes effect when run on a cluster. [default: false]

**`--snakemake-cmds`** STRING

  Additional commands passed through to snakemake as a single string, e.g. `--snakemake-cmds "--print-compilation True"`. Most `snakemake -h` commands are valid, but some may clash with commands aviary supplies directly — check for conflicts before using.

## Examples

Run the full pipeline with paired reads:
```
aviary complete -1 reads_1.fq.gz -2 reads_2.fq.gz
```

Run the full pipeline with long reads:
```
aviary complete -1 reads_1.fq.gz -2 reads_2.fq.gz --longreads reads.fastq.gz --long_read_type ont
```

Start from an existing assembly:
```
aviary complete --assembly scaffolds.fasta -1 reads_1.fq.gz -2 reads_2.fq.gz
```

Recover from multiple assemblies with SemiBin2 multi-sample binning:
```
aviary complete --assembly sample1.fasta sample2.fasta \
  -1 sample1_1.fq.gz sample2_1.fq.gz -2 sample1_2.fq.gz sample2_2.fq.gz \
  --semibin-mode multi
```
