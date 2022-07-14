---
title: Usage
---

Usage
========

## Basic Usage

To perform hybrid assembly:
```
aviary assemble -1 *.1.fq.gz -2 *.2.fq.gz --longreads *.nanopore.fastq.gz --long_read_type ont -t 24 -n 48
```
Aviary is compatible with both Nanopore and PacBio long read technologies. 
Note: Aviary can also perform assembly using just short or long reads as well.
```
aviary assemble -1 *.1.fq.gz -2 *.2.fq.gz -t 24 -n 48

OR

aviary assemble --longreads *.nanopore.fastq.gz --long_read_type ont -t 24 -n 48
```


To perform mag recovery:
```
aviary recover --assembly scaffolds.fasta -1 sr1.1.fq sr2.1.fq.gz -2 sr1.2.fq sr2.2.fq.gz --longreads nanopore.fastq.gz -z ont --output output_dir/ --max_threads 12 --n_cores 24 --gtdb_path /path/to/gtdb/release/
```
If no assembly file is provided, then aviary will first perform the assembly pipeline to produce an assembly using the 
input reads.

If at any point the Aviary workflow is interrupted, the pipeline can be restarted and pick up from the last completed
step.

## Advanced Usage

Often users are required to send long running jobs off on to high performance clusters. Aviary and snakemake are
perfectly compatible with clusters and can be sent off as either a single pipeline via PBS script or equivalent.
Alternatively, snakemake can send individual jobs in a pipeline off into a cluster to share the load across nodes. 
You can make use of this feature in Aviary via the `--snakemake-cmds` parameter, E.g.
```
aviary assemble -1 *.1.fq.gz -2 *.2.fq.gz --longreads *.nanopore.fastq.gz --long_read_type ont -t 24 -p 24 -n 24 --snakemake-cmds '--cluster qsub '
```
NOTE: The space after `--cluster qsub ` is required due to a strange quirk in how python's `argparse` module works.

## Helpful parameters and commands

### Environment variables
Upon first running Aviary, you will be prompted to input the location for several database folders if
they haven't already been provided. If at any point the location of these folders change you can
use the the `aviary configure` module to update the environment variables used by aviary.

These environment variables can also be configured manually, just set the following variables in your `.bashrc` file:
```
GTDBTK_DATA_PATH
BUSCO_DB
CONDA_ENV_PATH
```

Make sure to reactivate your conda environment or re-source your `.bashrc` for aviary to be able to access these variables.

### Thread control
Aviary has three thread contol options:

#### `-t, --threads`

- Controls how many threads any individual program can use. If set to 24, then each program that has threading options 
will use 24 threads when they run.

#### `-n, --n-cores, --n_cores`

- Controls how many cores (CPUs) snakemake will be given. If this value is set higher than `--threads`, then potentially
multiple programs will run concurrently providing a great boost in performance. If this value is not set then it defaults 
to being the same value as `--threads`

#### `-p, --pplacer-threads, --pplacer_threads`

- Pplacer is special and gets its own thread parameter. Why? Because it randomly deadlocks when given too many threads and 
can also be kind of memory intensive when given extra threads.

### RAM control

When performing assembly, users are required to estimate how much RAM they will need to use via `-m, --max-memory, --max_memory`

### Workflow control

Often users may not want to run a complete aviary module, as such specific rules can be targeted via the `-w, --workflow`
parameter. For example, if a user wanted to only run a specific binning algorithm then that rule can be specified directly:
```
aviary recover -w rosella --assembly scaffolds.fasta -1 sr1.1.fq sr2.1.fq.gz -2 sr1.2.fq sr2.2.fq.gz --longreads nanopore.fastq.gz --output output_dir/ --max_threads 12
```
NOTE: Every step up to the targeted rule still has to be run if it hasn't been run before. The specific rules that can be 
used can be found within each modules specific snakemake file.