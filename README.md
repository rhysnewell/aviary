# Aviary
A snakemake pipeline for binning metagenomic assemblies

# Installation

```
git clone https://github.com/rhysnewell/aviary.git
cd aviary
conda env create -n aviary -f aviary.yml
conda activate aviary
pip install -e .
aviary recover --help
```

# Requirements

Your conda channels should be configured ideally in this order:
```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

Initial requirements for aviary can be downloaded using the `aviary.yml`:
```
conda env create -n aviary -f aviary.yml
```

# Usage

To perform mag recovery:
```
aviary recover --assembly scaffolds.fasta -1 sr1.1.fq sr2.1.fq.gz -2 sr1.2.fq sr2.2.fq.gz --longreads nanopore.fastq.gz --output output_dir/ --max_threads 12 --n_cores 24 --gtdb_path /path/to/gtdb/release/
```
