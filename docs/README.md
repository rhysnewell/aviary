![](/images/lorikeet_logo.png)

![](https://travis-ci.com/rhysnewell/Lorikeet.svg?branch=master)
![](https://anaconda.org/bioconda/lorikeet-genome/badges/license.svg)
![](https://anaconda.org/bioconda/lorikeet-genome/badges/version.svg)
![](https://anaconda.org/bioconda/lorikeet-genome/badges/latest_release_relative_date.svg)
![](https://anaconda.org/bioconda/lorikeet-genome/badges/platforms.svg)


Lorikeet is a within-species variant analysis pipeline for metagenomic communities that utilizes both long and short reads.
Lorikeet utilizes a re-implementaion of the GATK HaplotypeCaller algorithm, performing local re-assembly of potentially active
regions within candidate genomes. Called variants can be clustered into likely strains using a combination of UMAP and HDBSCAN.
Additional statistics, like consensus ANI, population ANI, and subpopulation ANI will also be calculated for each input
geome providing values for each sample compared to the reference and also compared to all other samples.

Lorikeet has a variety of subcommands with the main being `call` and `genotype`. The `call` pipeline will take any number
of input genomes and samples and perform robust variant calling and ANI calculations. The `genotype` algorithm takes this
a step further and attempts to reconstruct strain haplotypes from the called variants and return complete strain genomes.

## Additional resources

The variant calling algorithm is basically a one-to-one re-implementation of the algorithm used in GATK HaplotypeCaller.
As such, many of the FAQs and documentation for HaplotypeCaller can be useful in understanding how Lorikeet actually
finds variants. An overview of the HaplotypeCaller pipeline can be found here: [HaplotypeCaller Docs](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller).

Lorikeet makes use of a couple new and daunting algorithms. UMAP in particular is an amazing algorithm but might be cause 
for concern since it is difficult to understand how it works and what it is doing. So please look over this amazing article 
by Andy Coenen and Adam Pearce: [Understanding UMAP](https://pair-code.github.io/understanding-umap/)

## Citation

On its way :P

## License

Code is [GPL-3.0](LICENSE)