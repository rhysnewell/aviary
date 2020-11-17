# BinSnek
A snakemake pipeline for binning metagenomic assemblies

# Installation

```
git clone https://github.com/rhysnewell/BinSnek.git
cd BinSnek
conda env create -n binsnek -f binsnek.yml
conda activate binsnek
pip install --editable .
binsnek recover --help
```

# Requirements

Initial requirements for binsnek can be downloaded using the `binsnek.yml`:
```
conda env create -n binsnek -f binsnek.yml
```

# Usage

To perform mag recovery:
```
binsnek recover --assembly scaffolds.fasta --short_reads_1 sr1.1.fq sr2.1.fq.gz --short_reads_2 sr1.2.fq sr2.2.fq.gz --longreads nanopore.fastq.gz --output output_dir/ --max_threads 24 --n_cores 24 --gtdb_path /path/to/gtdb/release/
```

# Batch Files

Instead of providing BinSnek with an assembly and reads, you can provide it a batch file in the following format:

```
/Absolute/Path/to/Assembly1.fasta    Unique_ID_1    /absolute/path/to/read_set_1.1.fq.gz    /absolute/path/to/read_set_1.2.fq.gz    /absolute/path/to/read_set_2.1.fq.gz    /absolute/path/to/read_set_2.2.fq.gz
/Absolute/Path/to/Assembly2.fasta    Unique_ID_2    /absolute/path/to/read_set_3.1.fq.gz    /absolute/path/to/read_set_3.2.fq.gz    /absolute/path/to/read_set_4.1.fq.gz    /absolute/path/to/read_set_4.2.fq.gz
```

Then specify the path to your batch file in the config.yaml, ignoring inputs for fasta and any reads, and use the following command:
`snakemake --use-conda --cores 24 run_batch`
