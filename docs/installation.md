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
conda config --add channels conda-forge
conda config --add channels bioconda
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

Initial requirements for aviary can be downloaded using the `aviary.yml`:

```
git clone https://github.com/rhysnewell/aviary.git
cd aviary
conda env create -n aviary -f aviary.yml
conda activate aviary
pip install -e .
aviary --help
```

The resulting output should contain a list of the available aviary modules:
```

                    ......:::::: AVIARY ::::::......

A comprehensive metagenomics bioinformatics pipeline

Metagenome assembly, binning, and annotation:
        cluster   - Clusters and dereplicates bins across multiple aviary runs
        assemble  - Perform hybrid assembly using short and long reads,
                    or assembly using only short reads
        recover   - Recover MAGs from provided assembly using a variety
                    of binning algorithms
        annotate  - Annotate MAGs **TBC**
        genotype  - Perform strain level analysis of MAGs **TBC**
        complete  - Runs each stage of the pipeline: assemble, recover,
                    annotate, genotype in that order.

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

Aviary uses programs which require access to locally stored databases. These databases can be quite large, as such we recommend setting up one instance of Aviary and these databases per machine or machine cluster.

The **required** databases are as follows:
* [GTDB](https://gtdb.ecogenomic.org/downloads) Required for taxonomic annotation

The **optional** databases are as follows:
* [EggNog](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.7#setup) Will become required soon.

**If you do not have the optional databases installed, then when aviary asks you to specify these databse passes when configuring just press enter and specify no path.**
