---
title: FAQs & Troubleshooting
---

FAQs & Troubleshooting
========

This page is just meant for general questions that I notice are asked with some frequency. If you feel like something
is missing from here and you'd like to see it included, feel free to ask it by raising an issue on GitHub.

### qsub and pysam - ModuleNotFoundError

A known issue with using snakemake + pysam + qsub results in the a break in the pipeline. The issue arises because pysam 
does not activate correctly when using qsub by default. To fix this you just need to add the `-V ` parameter to your qsub
command.

### Which databases do I download?

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