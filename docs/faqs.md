---
title: FAQs & Troubleshooting
---

FAQs & Troubleshooting
========

This page is just meant for general questions that I notice are asked with some frequency. If you feel like something
is missing from here and you'd like to see it included, feel free to ask it by raising an issue on GitHub.


### I want to perform both MAG recovery and assembly, how do I do that?

If you supply all reads to the `recover` command then Aviary will perform assembly first and then perform MAG recovery. You can perform assembly first by using the `assemble` command and then using the output assembly file as the input to the `recover` command. However, it is much simpler to let aviary handle this process for you.

### An error occurred but I don't know where to look for the error message

All of the error messages are stored in the `logs/` folder. The error messages are stored in the `logs/` folder with the name of the file corresponding the name of the rule that they originate from. If you had an error occur in the `qc_short_reads` rule then consult the `logs/qc_short_reads.log` file for the error message.

### Is Aviary cluster compatible?

Yes! Consult the examples page for more information.

### I have access to a GPU, can I use it?

Yes! Aviary supports the use of GPUs for the assembly process. If the GPU is on a local machine, you must first install the `cuda` package into your conda environment. Then, programs that use GPUs should automatically detect its presence.

If you are using a cluster, you can supply the `--request-gpu` flag and Aviary will attempt to place rules that use GPUs on to a machine that has GPUs available.

### Error in perpare_binning_files

This error is almost always caused by the user running out of storage in their `/tmp` folder when `coverm` performs the mapping process. To fix this, you can either increase the amount of storage available to the `/tmp` folder or you can change the location of the temporary folder by setting the `TMPDIR` environment variable to a folder with more storage. Aviary also allows the user to specify the location of the temporary folder by using the `--tmpdir` parameter.

### I wish to remove host contamination from my reads

Aviary supports the removal of host contamination during the assembly process via the `-r`, `--reference-filter` parameter. This flag can take one or more compressed or non-compressed fasta files. Aviary will then compare the reads to these references and remove any reads that map to them.

### SPAdes error: "Error code: -9" or other errors

The most likely solution to this is that you are running out of memory. SPAdes is a memory intensive program and will exit unexpectedly if it reaches the maximum memory limit of your machine or supplied by aviary.
To increase the amount of memory available to SPAdes, you can either increase the amount of memory available to the entire pipeline by using the `-m` parameter.



### qsub and pysam - ModuleNotFoundError

A known issue with using snakemake + pysam + qsub results in the a break in the pipeline. The issue arises because pysam 
does not activate correctly when using qsub by default. To fix this you just need to add the `-V ` parameter to your qsub
command.

### Which databases do I download?

It is probably best to just let Aviary handle the downloading of your databases via the `--download` parameter. But, if you
would like to set them up yourself, please read ahead

For the GTDB:
* [GTDB](https://gtdb.ecogenomic.org/downloads) Required for taxonomic annotation
Download and point the GTDB environment variable to the `db/` folder inside of that download.

The **optional** databases are as follows:
* [EggNog](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.7#setup).
Download this databse and point to the root folder of the databse.

Aviary will ask for the paths to these database files if they don't exist, otherwise you can place these lines into
the `activate.d/aviary.sh` or `.bashrc` files changing the specific paths:
```
export GTDBTK_DATA_PATH=/path/to/gtdb/gtdb_release207/db/ # https://gtdb.ecogenomic.org/downloads
export EGGNOG_DATA_DIR=/path/to/eggnog-mapper/2.1.7/ # https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.7#setup
export CONDA_ENV_PATH=/path/to/conda/envs/
```

### Why the name "Aviary"? Why the bird names in general?

Put all your birds in one place.

### Where did the logo come from?

I made it (among other bird based + CoverM logos) using [GIMP](https://www.gimp.org/) and based the idea off of this 
[tutorial](https://www.youtube.com/watch?v=fSOR7mPwb4I). They are very easy to make so just follow that video if you 
feel like making something similar.

### Where's the paper?

You sound like my supervisor.