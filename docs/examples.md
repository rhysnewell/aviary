---
title: Examples
---

Examples
========

## Forcing Variant Calls

Through use of the `-f, --features-vcf` argument it is possible to provide lorikeet with a list of known variants or variants of concern that will be forcibly called by
the algorithm if there is any activity occurring at that genomic location in the provided samples. It is possible that lorikeet
will also call any potential variation events surrounding a given location as well. This can be useful if you are more 
interested in the activity occurring around some given locations rather than specific variation events.

The list of variants must be provided in VCF format with some caveats on how the variant locations are written. E.g.
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##source=lorikeet-v0.6.0
##contig=<ID=random10000~random_sequence_length_10000_1,length=10000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER
random10000~random_sequence_length_10000_1	223	.	G	C	1445.63	.
random10000~random_sequence_length_10000_1	435	.	G	C	1941.63	.
random10000~random_sequence_length_10000_1	949	.	T	A	383.629	.
```

Note that you need to specific both the genome name and contig name in the `CHROM` column separated by the `~` character
 e.g. `random10000~random_sequence_length_10000_1`. Additionally, the standard `FORMAT` and `INFO` columns are optional when
 provding a VCF file.
 
Finally, the provided VCF file must be compressed using `bgzip` and indexed using `bcftools index` e.g.
```
bgzip -c random10000.vcf > random10000.vcf;
bcftools index random10000.vcf.gz;
lorikeet call -r random10000.fna -1 forward_reads.fastq -2 reverse_reads.fastq -l longread.bam -f random10000.vcf.gz
```

## Visualizing ANI results

The ANI tables that Lorikeet produces can feel overwhelming at first, but luckily they can be fairly easily visualized 
into interpretable and digestable heatmaps. The following is going to be a short walkthrough on the processes I have taken 
previously to generate heatmap results for a singular novel MAG.

Imagine we've used [SingleM](https://github.com/wwood/singlem) to search the SRA for hits for our novel organism and then
used [Kingfisher](https://github.com/wwood/kingfisher-download) to download all of the samples where our organism is apparently
appearing in excess of 5 times coverage and contains 10 or more gene hits. We would then pass all of these samples to Lorikeet
with our novel MAG to generate the variants, ANI tables, and optionally the strain genotypes as well. If we just wanted 
the variants and ANI tables we'd use this single command:

```
lorikeet call -r novel_mag.fna -c sra_samples/*_[12].fastq* -t 30 --bam-file-cache-directory bams/ -o lorikeet_out/
```

Since this example involves mapping over 90 samples to our reference genome, I've chosen to use the `--bam-file-cache-directory`
flag so that the read mappins are saved and we can reuse them if we experience a crash or wish to change the settings. If we 
wanted to rerun Lorikeet using the BAM files generate from the mappings then we'd just change our command to:

```
lorikeet call -r novel_mag.fna -b bams/short/*.bam -t 30 -o lorikeet_out/
```

Once Lorikeet has finished up, we should have an output folder filled with ANI tables in `lorikeet_out/novel_mag/`. We then
want to use some R libraries like `data.table`, `ggplot2`, and `plotly` to visualize what the results look like.

A simple example of this would look like:

```R
library(data.table)
library(ggplot2)
library(plotly)

pop_ani <- fread("lorikeet_out/novel_mag/novel_mag_population_ani.tsv", skip="SampleID",) # skip here so we skip all the contig header info
con_ani <- fread("lorikeet_out/novel_mag/novel_mag_consensus_ani.tsv", skip="SampleID")
sub_pop_ani <- fread("lorikeet_out/novel_mag/novel_mag_population_ani.tsv", skip="SampleID",)

## An example using consensus ANI

# sample_labels
con_ani[, SampleID:=factor(SampleID)]

tmp <-heatmaply(con_ani[, -1, with=FALSE],
                  xlab = "Sample Name",
                  ylab = "Sample Name",
                  file = "con_ani_heatmap.html",
                  k_row = 1, k_col = 1, plot_method="ggplot")

```

This will produce a `HTML` file using `plotly` that should look something like [heatmap](/figures/con_ani_heatmap.html) where
dark cells indicate more distant ANI values and lighter cells represent closely related samples. This example only explores
the use of consensus ANI, but we strongly recommend using both population and subpopulation ANI to assess how your communities
are differing.
