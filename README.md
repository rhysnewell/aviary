# BinSnek
A snakemake pipeline for binning metagenomic assemblies

# Requirements

You'll need `snakemake` installed. It is recommended you set up a `conda` environment for this:
```
conda create -n binsnek -c bioconda snakemake
```

# Usage

1. Copy the BinSnek pipeline into your working directory and activate your conda environment
`conda activate binsnek`

2. Setup your config.yaml file with the correct paths to your reads and assembly.

3. Run the pipeline:
`snakemake --use-conda --cores 24 recover_mags`
