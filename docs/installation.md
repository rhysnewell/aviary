---
title: Installation
---

Installation
========

## Requirements

Your conda channels should be configured ideally in this:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Your resulting `.condarc` file should look something like:
```
channels:
  - conda-forge
  - bioconda
  - defaults
```

#### Option 1: Install from Bioconda

Conda can handle the creation of the environment for you directly:

```
conda create -n aviary -c bioconda aviary
```

Or install into existing environment:
```
conda install -c bioconda aviary
```

#### Option 2: Install from pip

Create the environment using the `aviary.yml` file then install from pip:
```
conda env create -n aviary -f aviary.yml
conda activate aviary
pip install aviary-genome
```

#### Option 3: Install from source

To install from source, we recommend using [pixi](https://pixi.sh/). First clone
the aviary repository from GitHub:
```
git clone https://github.com/rhysnewell/aviary.git
cd aviary
```

Then install the main environment using pixi:
```
pixi run --manifest-path aviary/pixi.toml postinstall
```

Then aviary can be run using `pixi run` (or via `pixi shell`).
```
pixi run --manifest-path aviary/pixi.toml aviary --help
```

When installed this way, aviary is installed in an "editable" way (similar to `pip install -e .`), meaning that any changes made to aviary source are immediately available via the `aviary` command. This is useful for development and debugging.

When run this way, the databases required for aviary (e.g. `CHECKM2DB`) can be symlinked from a `db/` directory in the aviary repository. An activation hook then ensures that these are available when in the pixi environments. To do this, create a `db/` directory in the aviary repository and symlink the required databases into it. For example, as of writing:
```
$ ls db -l
lrwxrwxrwx - woodcrob 23 Apr 07:56 2.1.3 -> /mnt/hpccs01/work/microbiome/db/eggnog-mapper/2.1.3
lrwxrwxrwx - woodcrob 23 Apr 07:55 2015_01_16_v2 -> /work/microbiome/db/checkm/2015_01_16_v2
lrwxrwxrwx - woodcrob 23 Apr 07:54 CheckM2_database -> /work/microbiome/db/CheckM2_database
lrwxrwxrwx - woodcrob 23 Apr 07:57 gtdb207 -> /work/microbiome/db/metabuli/gtdb207
lrwxrwxrwx - woodcrob 23 Apr 07:56 release220 -> /work/microbiome/db/gtdb/gtdb_release220/auxillary_files/gtdbtk_package/full_package/release220
lrwxrwxrwx - woodcrob 23 Apr 07:55 S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb -> /work/microbiome/db/singlem/S4.3.0.GTDB_r220
```
To check the expected database symlink names, see `admin/set_env_vars.sh` in the
aviary repository. The advantage of this approach is that locations of the
databases are not tracked in the repository, since they are specific to the
computing cluster of the user.

## Databases

Aviary uses programs which require access to locally stored databases. 
These databases can be quite large, as such we recommend setting up one instance of Aviary and these databases per machine or machine cluster.

The **required** databases are as follows:
* [GTDB](https://gtdb.ecogenomic.org/downloads)
* [EggNog](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.8#setup)
* [CheckM2](https://github.com/chklovski/CheckM2)

### Installing databases

Aviary can handle the download and installation of these databases via use of the `--download` flag. Using `--download`
will download and install the databases into the folders corresponding to their associated environment variables. Aviary will
ask you to set these environment variables upon first running and if they are not already available. Otherwise, users can use
the `aviary configure` subcommand to reset the environment variables:

```commandline
aviary configure -o logs/ --eggnog-db-path /shared/db/eggnog/ --gtdb-path /shared/db/gtdb/ --checkm2-db-path /shared/db/checkm2db/ --download
```

This command will check if the databases exist at those given locations, if they don't then aviary will download and change
the conda environment variables to match those paths. 

**N.B.** Again, these databases are VERY large. Please talk to your sysadmin/bioinformatics specialist about setting a shared
location to install these databases to prevent unnecessary storage use. Additionally, the `--download` flag can be used within
any aviary module to check that databases are configured properly.

### Environment variables

Upon first running Aviary, you will be prompted to input the location for several database folders if
they haven't already been provided. If at any point the location of these folders change you can
use the the `aviary configure` module to update the environment variables used by aviary.

These environment variables can also be configured manually, just set the following variables in your `.bashrc` file:
```
export GTDBTK_DATA_PATH=/path/to/gtdb/gtdb_release207/db/ # https://gtdb.ecogenomic.org/downloads
export EGGNOG_DATA_DIR=/path/to/eggnog-mapper/2.1.7/ # https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.7#setup
export SINGLEM_METAPACKAGE_PATH=/path/to/singlem_metapackage.smpkg/
export CHECKM2DB=/path/to/checkm2db/
```
