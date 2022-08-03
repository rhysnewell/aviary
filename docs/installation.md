---
title: Installation
---

Installation
========

## Requirements

Your conda channels should be configured ideally in this order with strict channel priority order
turned on:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Your resulting `.condarc` file should look something like:
```
channels:
  - conda-forge
  - bioconda
  - defaults
channel_priority: strict
```

#### Option 1) Install from Bioconda

Conda can handle the creation of the environment for you directly:

```
conda create -n aviary -c bioconda aviary
```

Or install into existing environment:
```
conda install -c bioconda aviary
```

#### Option 2) Install from pip

Create the environment using the `aviary.yml` file then install from pip:
```
conda env create -n aviary -f aviary.yml
conda activate aviary
pip install aviary-genome
```

#### Option 3) Install from source

Initial requirements for aviary can be downloaded using the `aviary.yml`:
```
git clone https://github.com/rhysnewell/aviary.git
cd aviary
conda env create -n aviary -f aviary.yml
conda activate aviary
pip install -e .
```

Whatever option you choose, running `aviary --help` should return the following
output:
```

                    ......:::::: AVIARY ::::::......

           A comprehensive metagenomics bioinformatics pipeline

Metagenome assembly, binning, and annotation:
        assemble  - Perform hybrid assembly using short and long reads, 
                    or assembly using only short reads
        recover   - Recover MAGs from provided assembly using a variety 
                    of binning algorithms 
        annotate  - Annotate MAGs using EggNOG and GTBD-tk
        genotype  - Perform strain diversity analysis of MAGs using Lorikeet
        complete  - Runs each stage of the pipeline: assemble, recover, 
                    annotate, genotype in that order.
        cluster   - Combines and dereplicates the MAGs from multiple Aviary runs
                    using Galah

Isolate assembly, binning, and annotation:
        isolate   - Perform isolate assembly **PARTIALLY COMPLETED**
        
Utility modules:
        configure - Set or overwrite the environment variables for future runs.

```

Upon first running aviary you will be prompted to input the location for where you would like
your conda environments to be stored, the GTDB release installed on your system, the location of your
EggNog database, and the location of your BUSCO database. These locations will be stored as environment
variables, and automatically sourced by Aviary at runtime.

These environment variables can be reset using `aviary configure`

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
export CONDA_ENV_PATH=/path/to/conda/envs/
```
