---
title: aviary recover
---

# aviary recover

Recover metagenome-assembled genomes (MAGs) from an assembly using multiple binning algorithms, followed by quality assessment and taxonomic classification.

```bash
aviary recover --assembly scaffolds.fasta -1 reads_1.fq.gz -2 reads_2.fq.gz
```

If no assembly is provided, aviary will first run the assembly pipeline.

> This subcommand also accepts the common workflow, resource, output and
> execution options described under [Shared options](centralised_commands.md).

| Jump to | |
| --- | --- |
| [Input options](#input-options) | Read files and assemblies |
| [Assembly options](#assembly-options) | Assembler choice, k-mers, Flye contig filtering |
| [QC options](#qc-options) | Read filtering and quality thresholds |
| [Binning options](#binning-options) | Binner selection, SemiBin2 mode/model, refinement, quality cutoffs |
| [Annotation / bin processing options](#annotation-bin-processing-options) | Local database paths |
| [Performance options](#performance-options) | Threads, cores, memory |
| [Output options](#output-options) | Output and temp directories |
| [Misc options](#misc-options) | Cluster profile, dry run, GPU request |
| [Examples](#examples) | Common flag combinations |

## Examples

The basic shape — an existing assembly plus reads for coverage:
```bash
aviary recover --assembly scaffolds.fasta -1 reads_1.fq.gz -2 reads_2.fq.gz
```

If `--assembly` is omitted, aviary assembles first from the same reads before binning; see
[`aviary assemble`](assemble.md) for the read-input shapes (`-1`/`-2`, `-i`, `-c`, `-l`) and
assembler choices (`--use-megahit`, `--long-read-assembler`) that apply equally here.

The rest of these examples focus on what `recover` itself controls: which binners run, how
SemiBin2 is configured, and what happens after binning.

### Choosing which binning algorithms run

`rosella`, `semibin`, `metabat1`/`metabat2` and `vamb` run by default; `maxbin2`, `concoct`,
`comebin`, `taxvamb` and `quickbin` are skipped by default because of their runtime. Turn extras
on, or drop defaults you don't want (e.g. `vamb` is memory-hungry on large assemblies):
```bash
aviary recover --assembly scaffolds.fasta -1 reads_1.fq.gz -2 reads_2.fq.gz \
  --extra-binners concoct comebin

aviary recover --assembly scaffolds.fasta -1 reads_1.fq.gz -2 reads_2.fq.gz \
  --skip-binners vamb metabat
```

### Matching the SemiBin2 model to the sample's environment

`--semibin-model` selects a pre-trained SemiBin2 environment model instead of the generic
`global` default, which usually improves binning when the sample's environment is one SemiBin2
has a model for:
```bash
aviary recover --assembly gut_sample.fasta -1 reads_1.fq.gz -2 reads_2.fq.gz \
  --semibin-model human_gut

aviary recover --assembly seawater_sample.fasta -1 reads_1.fq.gz -2 reads_2.fq.gz \
  --semibin-model ocean
```

For a time series or multiple related samples, co-bin them instead of picking one model —
SemiBin2 learns cross-sample abundance correlation, which `--semibin-model` doesn't need for
this mode and is ignored:
```bash
aviary recover --assembly week1.fasta week2.fasta week3.fasta \
  -1 week1_1.fq.gz week2_1.fq.gz week3_1.fq.gz \
  -2 week1_2.fq.gz week2_2.fq.gz week3_2.fq.gz \
  --semibin-mode multi
```

### Tuning MAG quality thresholds and refinement

Raise the completeness bar and tighten contamination for a stricter final bin set (defaults are
50% / 5%):
```bash
aviary recover --assembly scaffolds.fasta -1 reads_1.fq.gz -2 reads_2.fq.gz \
  --min-completeness 70 --max-contamination 5
```

Rosella's refinement step iterates by default; disable it for a faster, less-refined pass, or
give it more retries on a difficult assembly:
```bash
aviary recover --assembly scaffolds.fasta -1 reads_1.fq.gz -2 reads_2.fq.gz \
  --refinery-max-iterations 0

aviary recover --assembly scaffolds.fasta -1 reads_1.fq.gz -2 reads_2.fq.gz \
  --refinery-max-retries 6
```

### Large sample cohorts

Coverage calculation is normally submitted as one job per rule; with many samples, splitting it
into smaller parallel jobs shortens wall time on a cluster:
```bash
aviary recover --assembly scaffolds.fasta -1 s01_1.fq.gz ... s40_1.fq.gz -2 s01_2.fq.gz ... s40_2.fq.gz \
  --coverage-job-strategy always --coverage-samples-per-job 8
```

### Skip individual downstream steps

Stop after binning (skip SingleM, GTDB-tk and CoverM abundance calculation entirely):
```bash
aviary recover --assembly scaffolds.fasta -1 reads_1.fq.gz -2 reads_2.fq.gz --binning-only
```

Keep binning and abundance, but skip taxonomy assignment specifically (useful when GTDB-tk's
database isn't configured yet, or its runtime isn't needed for a given run):
```bash
aviary recover --assembly scaffolds.fasta -1 reads_1.fq.gz -2 reads_2.fq.gz --skip-taxonomy
```

## Input options

**`-a`**, **`--assembly`** FILE [FILE ...]

  One or more FASTA files containing scaffolded contigs of metagenome assemblies. If not provided, aviary will assemble first. Provide multiple assemblies for SemiBin2 multi-sample binning (requires `--semibin-mode multi`; see [Core concepts](../concepts.md) if binning and assembly terms are new to you).

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

> Reads supplied above are mapped back onto the assembly for coverage,
> abundance and polishing. See [Read mappers](mappers.md) for the
> `--short-read-mapper`/`--long-read-mapper` family of flags that control
> which aligner is used.

## Assembly options

**`--long-read-assembler`** ASSEMBLER

  Long-read assembler: `myloasm` or `flye`. [default: myloasm]

**`--medaka-model`** MODEL

  Medaka model for long-read polishing. [default: r941_min_hac_g507]

**`--use-unicycler`**

  Use Unicycler to re-assemble the metaSPAdes hybrid assembly. Not recommended for complex metagenomes.

**`--use-megahit`**

  Use MEGAHIT instead of metaSPAdes for short-read-only assembly. [default: false]

**`--coassemble`**, **`--co-assemble`**

  When multiple read sets are given, coassemble them together. If false, aviary uses only the first short-read and first long-read set for assembly (all read sets are still used for differential-coverage binning). [default: false]

**`-k`**, **`--kmer-sizes`** INT [INT ...]

  Manually specify the k-mer sizes used by SPAdes during assembly. Space-separated odd integers less than 128, or `auto`. [default: auto]

These five flags all decide which Flye contigs make it into the assembly before binning. They interact, so read them as one rule set rather than five independent options:

| Flag | Compares | Rule | Default |
|---|---|---|---|
| `--min-cov-long` | long-read coverage | keep if coverage ≥ value | 5 |
| `--min-cov-short` | short-read coverage | keep if coverage ≤ value | 5 |
| `--exclude-contig-cov` | long-read coverage | drop **only if** coverage ≤ value **and** length ≤ `--exclude-contig-size` | 10 |
| `--exclude-contig-size` | contig length | drop **only if** length ≤ value **and** coverage ≤ `--exclude-contig-cov` | 2500 |
| `--include-contig-size` | contig length | keep if length ≥ value (checked first, before the exclude/include-coverage rules below even apply) | 10000 |

!!! example "How these combine"
    Aviary checks `--include-contig-size` first: any contig at or above this length (or marked circular by Flye) is kept immediately, regardless of coverage. Only shorter contigs fall through to the coverage rules — a 2,000 bp contig at 8x long-read coverage is then dropped (length ≤ 2500 **and** coverage ≤ 10), but the same contig at 12x coverage is kept (fails the coverage half of the exclude rule). A 15,000 bp contig is always kept, since it clears `--include-contig-size` before the exclude/coverage checks are ever reached.

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

??? note "Show all 13 options"

    **`-r`**, **`--host-filter`** FILE [FILE ...]

      Host reference FASTA files for removal of contaminant reads prior to assembly.

    **`-g`**, **`--gold-standard-assembly`** FILE [FILE ...]

      A gold-standard assembly to compare the resulting (or a given input) assembly against.

    **`--gsa-mappings`** FILE

      CAMI I & II gold-standard-assembly mappings, used alongside `--gold-standard-assembly`.

    **`--keep-percent`** INT

    !!! warning "Deprecated"
        Percentage of reads passing quality thresholds kept by Filtlong. [default: 100]

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

  Extra binning algorithms to run: `maxbin`/`maxbin2` (equivalent), `concoct`, `comebin`, `taxvamb`, `quickbin`. These are skipped by default due to long runtimes.

**`--skip-binners`** BINNER [BINNER ...]

  Binning algorithms to skip: `rosella`, `semibin`, `metabat1`, `metabat2`, `metabat`, `vamb`, `quickbin`.

**`--semibin-model`** MODEL

  SemiBin2 environment model: `human_gut`, `dog_gut`, `ocean`, `soil`, `cat_gut`, `human_oral`, `mouse_gut`, `pig_gut`, `built_environment`, `wastewater`, `global`. [default: global]

**`--semibin-mode`** `{single,multi}`

  SemiBin2 mode to use. `single` runs `single_easy_bin` on one assembly at a time. `multi` runs `multi_easy_bin`, co-binning multiple assemblies together — pass two or more files to `--assembly` to use it. Multi mode ignores `--semibin-model`, as pre-trained environments aren't supported for multi-sample binning. [default: single]

  ```bash
  aviary recover --assembly sample1.fasta sample2.fasta sample3.fasta \
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

  Skip SingleM post-binning recovery assessment. [default: true]

**`--min-completeness`** FLOAT

  Minimum CheckM2 completeness percentage for bins retained in downstream
  processing. [default: 50.0]

**`--max-contamination`** FLOAT

  Maximum CheckM2 contamination percentage for bins retained in downstream
  processing. [default: 5.0]

**`--coverage-job-strategy`** STRATEGY

  Strategy for coverage calculation across many samples: `default`, `never`, `always`. [default: default]

**`--coverage-samples-per-job`** INT

  Number of samples per coverage job when splitting. [default: 5]

**`--min-percent-read-identity-short`** FLOAT

  Minimum percent read identity used by CoverM for short reads when calculating genome abundances. [default: 95]

**`--min-percent-read-identity-long`** FLOAT

  Minimum percent read identity used by CoverM for long reads when calculating genome abundances. [default: 85]

## Annotation / bin processing options

??? note "Show all 5 options"

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

??? note "Show all 5 options"

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

??? note "Show all 2 options"

    **`-o`**, **`--output`** DIR

      Output directory. [default: ./]

    **`--tmpdir`** DIR

      Temporary files directory.

## Misc options

??? note "Show all 7 options"

    **`--snakemake-profile`** PROFILE

      Snakemake profile for cluster submission. See the Guides section for HPC usage.

    **`--cluster-retries`** INT

      Retries for failed cluster jobs. [default: 0]

    **`--dry-run`**

      Perform a snakemake dry run.

    **`--clean`**

      Clean up temporary files. [default: true]

    **`--strict`**

      Ensure each binner completes successfully. [default: skip failing binners]

    **`--request-gpu`**

      Request a GPU for the pipeline (taxvamb, comebin, semibin). Only takes effect when run on a cluster. [default: false]

    **`--snakemake-cmds`** STRING

      Additional commands passed through to snakemake as a single string, e.g. `--snakemake-cmds "--print-compilation True"`. Most `snakemake -h` commands are valid, but some may clash with commands aviary supplies directly — check for conflicts before using.
