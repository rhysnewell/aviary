---
title: aviary recover
---

# aviary recover

Recover metagenome-assembled genomes (MAGs) from an assembly using multiple binning algorithms, followed by quality assessment and taxonomic classification.

```
aviary recover --assembly scaffolds.fasta -1 reads_1.fq.gz -2 reads_2.fq.gz
```

If no assembly is provided, aviary will first run the assembly pipeline.

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

  Skip SingleM post-binning recovery assessment. [default: True]

**`--min-completeness`** FLOAT

  Minimum CheckM2 completeness percentage for annotation. [default: 50.0]

**`--max-contamination`** FLOAT

  Maximum CheckM2 contamination percentage for annotation. [default: 5.0]

**`--coverage-job-strategy`** STRATEGY

  Strategy for coverage calculation across many samples: `default`, `never`, `always`. [default: default]

**`--coverage-samples-per-job`** INT

  Number of samples per coverage job when splitting. [default: 5]

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

## Examples

Recover MAGs from an existing assembly with paired reads:
```
aviary recover --assembly scaffolds.fasta -1 reads_1.fq.gz -2 reads_2.fq.gz
```

Recover MAGs from an assembly with long reads:
```
aviary recover --assembly scaffolds.fasta --longreads reads.fastq.gz --long_read_type ont
```

Recover MAGs and assemble in one step:
```
aviary recover -1 reads_1.fq.gz -2 reads_2.fq.gz --longreads reads.fastq.gz --long_read_type ont
```
