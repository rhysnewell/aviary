---
title: Examples
---

Examples
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

## Batch Processing

Aviary allows users to supply a batch file to the `aviary batch` command. This will cause aviary to run on every line within
the input batch file individually. Example batch files can be found at [here](/examples/example_batch.tsv) and [here](/examples/example_batch.csv).

## Advanced Usage

Often users are required to send long running jobs off on to high performance clusters. Aviary and snakemake are
perfectly compatible with clusters and can be sent off as either a single pipeline via PBS script or equivalent.
Alternatively, snakemake can send individual jobs in a pipeline off into a cluster to share the load across nodes. 
You can make use of this feature in Aviary via the `--snakemake-cmds` parameter, E.g.
```
aviary assemble -1 *.1.fq.gz -2 *.2.fq.gz --longreads *.nanopore.fastq.gz --long_read_type ont -t 24 -p 24 -n 24 --snakemake-cmds '--cluster qsub '
```
NOTE: The space after `--cluster qsub ` is required due to a strange quirk in how python's `argparse` module works.

## HPC cluster submission

Often users are required to send long running jobs off on to high performance clusters. Aviary and snakemake are
perfectly compatible with clusters and can be sent off as either a single pipeline via PBS script or equivalent.
Alternatively, snakemake can send individual jobs in a pipeline off into a cluster to share the load across nodes.
You can make use of this feature in Aviary by providing a `--snakemake-profile cluster`, where `cluster` refers to
a Snakemake profile at `~/.config/snakemake/cluster/config.yaml`. An example `config.yaml` is provided below.
Additionally, `--cluster-retries` can be used to set the # of retries per job, with automatically increasing time
and memory requirements. CPU and memory bounds provided to Aviary are used as a hard-cap for job submission.

```yaml
cluster: qsub
cluster-status: qstat
jobs: 10000
cluster-cancel: qdel
```

Job resources were set based on [empirical data from 1,000 Aviary runs](https://github.com/rhysnewell/aviary/pull/270).
The number of CPUs was set based on the average mean load, to the nearest multiple of 8.
The number of CPUs for binners was standardized to 24, since this was the most common usage pattern.
Maximum memory was set based on the nearest power of 2, rounding up from the maximum RSS memory.
These resources can be restricted by setting `--max_memory` and `--max_threads` in Aviary.

## Expected output

Aviary will produce a variety of different outputs depending on the parameters provided. The following is a list of the expected outputs from the different subcommands.

In general, you will find rule benchmarks in the `benchmarks/` folder. Logs and error messages will be found in the `logs/` folder.

The `data/` folder contains various output from the different programs that have run. Aviary attempts to keep this folder clean, but there can be a lot of superfluous content in here. The `data/` folder is also where the final outputs from the pipeline will be found, they will be symlinked out of this folder into their respective output folders.

#### 1. Assemble

- `assembly/final_contigs.fasta` - The assembled genome in FASTA format
- `www/assembly_stats.txt` - A summary of the assembly statistics
- `www/fastqc` - A folder containing the output from fastqc
- `www/nanoplot` - A folder containing the output from NanoPlot

#### 2. Recover

If assembly is performed during the recover process, then the outputs from the assembly step will also be present.
- `bins/bin_info.tsv` - A file containing all of the information about the recovered MAGs. Taxonomy, QC, size, N50, etc.
- `bins/checkm_minimal.tsv` - A minimal version of the `bin_info.tsv` file.
- `bins/coverm_abundances.tsv` - The abundance of each MAG in each sample.
- `bins/final_bins` - A folder containing the final MAGs in FASTA format.
- `diversity/singlem_out` - A folder containing the output from singlem.
- `taxonomy/` - A folder containing the output from GTDB-tk

#### 3. Annotate

- `annotation/` - A folder containing the output from EggNOG-mapper.



## Helpful parameters and commands

### Environment variables
Upon first running Aviary, you will be prompted to input the location for several database folders if
they haven't already been provided. If at any point the location of these folders change you can
use the the `aviary configure` module to update the environment variables used by aviary.

These environment variables can also be configured manually, just set the following variables in your `.bashrc` file:
```
GTDBTK_DATA_PATH
EGGNOG_DATA_DIR
SINGLEM_METAPACKAGE_PATH
CHECKM2DB
CONDA_ENV_PATH
```

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

When performing assembly, users are required to estimate how much RAM they will need to use via `-m, --max-memory, --max_memory`.
With HPC cluster submission (see above), requested job memory is increased with each rerun and capped at `max_memory`.

### Temporary directory

By default, Aviary will use `/tmp` to store temporary files during many processes throughtout assembly and MAG recovery.
If you would like to use a different directory, you can specify this by using the flexible `--tmp` parameter.
If you would permanently like to change the temporary directory, you can use `aviary configure --tmp /new/tmp` to 
change the `TMPDIR` environment variable within your current conda environment. 

### Workflow control

Often users may not want to run a complete aviary module, as such specific rules can be targeted via the `-w, --workflow`
parameter. For example, if a user wanted to only run a specific binning algorithm then that rule can be specified directly:
```
aviary recover -w rosella --assembly scaffolds.fasta -1 sr1.1.fq sr2.1.fq.gz -2 sr1.2.fq sr2.2.fq.gz --longreads nanopore.fastq.gz --output output_dir/ --max_threads 12
```
NOTE: Every step up to the targeted rule still has to be run if it hasn't been run before. The specific rules that can be 
used can be found within each modules specific snakemake file.
