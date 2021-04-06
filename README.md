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

# Workflow
The current complete workflow for aviary. This is constantly being updated and will eventually include and assembly stage and
post binning analysis of MAGs
![Aviary workflow](figures/aviary_workflow.png)

# Citations
If you use aviary then please be aware that you are using a great number of other programs and aviary wrapping around them.
You should cite all of these tools as well, or whichever tools you know that you are using. To make this easy for you
we have provided the following list of citations for you to use in alphabetical order. This list will be updated as new
modules are added to aviary.

```
1.Kolmogorov, M., Yuan, J., Lin, Y. & Pevzner, P. A. Assembly of long, error-prone reads using repeat graphs. Nature Biotechnology 37, 540–546 (2019).
2.Alneberg, J. et al. Binning metagenomic contigs by coverage and composition. Nat Methods 11, 1144–1146 (2014).
3.Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P. & Tyson, G. W. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Res. 25, 1043–1055 (2015).
4.Hunt, M. et al. Circlator: automated circularization of genome assemblies using long sequencing reads. Genome Biology 16, 294 (2015).
5.Huerta-Cepas, J. et al. eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. Nucleic Acids Research 47, D309–D314 (2019).
6.Vaser, R., Sović, I., Nagarajan, N. & Šikić, M. Fast and accurate de novo genome assembly from long uncorrected reads. Genome Res 27, 737–746 (2017).
7.Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P. & Parks, D. H. GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. Bioinformatics 36, 1925–1927 (2020).
8.Nissen, J. N. et al. Improved metagenome binning and assembly using deep variational autoencoders. Nature Biotechnology 1–6 (2021) doi:10.1038/s41587-020-00777-4.
9.Olm, M. R. et al. inStrain profiles population microdiversity from metagenomic data and sensitively detects shared microbial strains. Nature Biotechnology 1–10 (2021) doi:10.1038/s41587-020-00797-0.
10.Kang, D. D. et al. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. PeerJ 7, (2019).
11.Li, H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34, 3094–3100 (2018).
12.De Coster, W., D’Hert, S., Schultz, D. T., Cruts, M. & Van Broeckhoven, C. NanoPack: visualizing and processing long-read sequencing data. Bioinformatics 34, 2666–2669 (2018).
13.Walker, B. J. et al. Pilon: An Integrated Tool for Comprehensive Microbial Variant Detection and Genome Assembly Improvement. PLOS ONE 9, e112963 (2014).
14.Sieber, C. M. K. et al. Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. Nature Microbiology 3, 836–843 (2018).
15.Köster, J. & Rahmann, S. Snakemake—a scalable bioinformatics workflow engine. Bioinformatics 28, 2520–2522 (2012).
16.Bankevich, A. et al. SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. J Comput Biol 19, 455–477 (2012).
17.Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078–2079 (2009).
18.Wick, R. R., Judd, L. M., Gorrie, C. L. & Holt, K. E. Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads. PLOS Computational Biology 13, e1005595 (2017).
```