![](/images/aviary_logo.png)

[![Anaconda-Server Badge](https://anaconda.org/bioconda/aviary/badges/license.svg)](https://anaconda.org/bioconda/aviary)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/aviary/badges/version.svg)](https://anaconda.org/bioconda/aviary)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/aviary/badges/latest_release_relative_date.svg)](https://anaconda.org/bioconda/aviary)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/aviary/badges/platforms.svg)](https://anaconda.org/bioconda/aviary)

# Aviary
An easy to use for wrapper for a robust snakemake pipeline for metagenomic hybrid assembly, binning, and annotation. 
The pipeline currently includes a step-down iterative 
hybrid assembler, an isolate hybrid assembler, a quality control module and a 
comprehensive binning pipeline. Each module can be run independently or as a single pipeline depending on provided input.

# Module details
|__method__ |__description__ |
| --- | --- |
|`cluster`|Dereplicate/choose representative genomes from multiple aviary runs|
|`assemble`|Perform quality control and assembly of provided reads. Will provide hybrid assembly if given long and short reads|
|`recover`|Recover MAGs from provided assembly using a variety of binning algorithms. Also perform quality checks on recovered MAGs and taxonomic classification.|
|`annotate`|Provide taxonomic and functional annotations for a given set of MAGs|
|`complete`|Performs the complete workflow up to last possible rule given the provided inputs|
|`isolate` |Performs hybrid isolate assembly. For use with isolated pure sequencing results.  |
|`configure` |Set or reset environment variables used by aviary  |


## Overview
![](/figures/aviary_workflow.png)

## Citation

If you use Aviary in your research, please cite:

> Newell RJP, Aroney STN, Zaugg J, Sternes P, Tyson GW, Woodcroft BJ.
> **Aviary: Hybrid assembly and genome recovery from metagenomes with Aviary.**
> Zenodo (2024). https://doi.org/10.5281/zenodo.10806928

## License

Code is GPL-3.0 

## GitHub

[Aviary](https://github.com/rhysnewell/aviary)