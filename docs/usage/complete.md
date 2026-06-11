---
title: aviary complete
---

# aviary complete

Performs all steps in the Aviary pipeline: Assembly → Binning → Refinement → Annotation.

```
aviary complete -1 reads_1.fq.gz -2 reads_2.fq.gz --longreads reads.fastq.gz --long_read_type ont
```

## Input options

**`-1`**, **`--pe-1`** FILE [FILE ...]

  Forward short read files.

**`-2`**, **`--pe-2`** FILE [FILE ...]

  Reverse short read files.

**`-l`**, **`--longreads`** FILE [FILE ...]

  Long-read files.

**`-z`**, **`--longread-type`** TYPE

  Sequencing platform: `rs`, `sq`, `ccs`, `hifi`, `ont`, `ont_hq`. [default: ont]

**`-a`**, **`--assembly`** FILE

  Optional pre-existing assembly FASTA to skip the assembly step.

## Assembly options

**`--long-read-assembler`** ASSEMBLER

  Long-read assembler: `myloasm` or `flye`. [default: myloasm]

**`--medaka-model`** MODEL

  Medaka model for long read polishing. [default: r941_min_hac_g507]

## Binning options

**`-s`**, **`--min-contig-size`** INT

  Minimum contig size in base pairs for binning. [default: 1500]

**`-b`**, **`--min-bin-size`** INT

  Minimum bin size in base pairs. [default: 200000]

**`--extra-binners`** BINNER [BINNER ...]

  Extra binning algorithms: `maxbin2`, `concoct`, `comebin`, `taxvamb`, `quickbin`.

**`--skip-binners`** BINNER [BINNER ...]

  Binning algorithms to skip: `rosella`, `semibin`, `metabat1`, `metabat2`, `metabat`, `vamb`, `quickbin`.

**`--semibin-model`** MODEL

  SemiBin environment model. [default: global]

**`--min-completeness`** FLOAT

  Minimum CheckM2 completeness for annotation. [default: 50.0]

**`--max-contamination`** FLOAT

  Maximum CheckM2 contamination for annotation. [default: 5.0]

## Annotation / bin processing options

**`--gtdb-path`** PATH

  Path to local GTDB database files.

**`--eggnog-db-path`** PATH

  Path to local EggNOG database files.

**`--checkm2-db-path`** PATH

  Path to CheckM2 database.

## Performance options

**`-t`**, **`--max-threads`** INT

  Maximum threads per process. [default: 8]

**`-n`**, **`--n-cores`** INT

  Maximum cores available. [default: 16]

**`-m`**, **`--max-memory`** INT

  Maximum memory in gigabytes. [default: 250]

## Output options

**`-o`**, **`--output`** DIR

  Output directory. [default: ./]

**`--tmpdir`** DIR

  Temporary files directory.

## Misc options

**`--snakemake-profile`** PROFILE

  Snakemake profile for cluster submission.

**`--dry-run`**

  Perform a snakemake dry run.

**`--clean`**

  Clean up temporary files. [default: True]

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
