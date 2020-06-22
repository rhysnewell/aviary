# BinSnek
A snakemake pipeline for binning metagenomic assemblies

# Requirements

You'll need `snakemake` installed. It is recommended you set up a `conda` environment for this:
```
conda create -n binsnek -c bioconda snakemake conda
```

# Usage

1. Copy the BinSnek pipeline into your working directory and activate your conda environment
`conda activate binsnek`

2. Setup your config.yaml file with the correct paths to your reads and assembly.

3. Run the pipeline:
`snakemake --use-conda --cores 24 recover_mags`

# Batch Files

Instead of providing BinSnek with an assembly and reads, you can provide it a batch file in the following format:

```
/Absolute/Path/to/Assembly1.fasta    Unique_ID_1    /absolute/path/to/read_set_1.1.fq.gz    /absolute/path/to/read_set_1.2.fq.gz    /absolute/path/to/read_set_2.1.fq.gz    /absolute/path/to/read_set_2.2.fq.gz
/Absolute/Path/to/Assembly2.fasta    Unique_ID_2    /absolute/path/to/read_set_3.1.fq.gz    /absolute/path/to/read_set_3.2.fq.gz    /absolute/path/to/read_set_4.1.fq.gz    /absolute/path/to/read_set_4.2.fq.gz
```

Then specify the path to your batch file in the config.yaml, ignoring inputs for fasta and any reads, and use the following command:
`snakemake --use-conda --cores 24 run_batch`
