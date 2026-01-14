---
title: Concepts
---

Concepts
========

Aviary provides the user with a lot of different outputs for each genome and some of the information present
in those outputs requires.

# File types

## FASTA

A FASTA file is a text based file used to represent either genomic nucleotide sequences or amino acids. It consists of 
of headers (lines starting with `>`) and blocks of sequences immediately following the headers. Fasta files are the format
used for the input reference/MAGs that Aviary uses and creates. The extension for such files is usually `.fasta`, `.fa`, or `.fna`.

For more info refer to the [wikipedia article](https://en.wikipedia.org/wiki/FASTA_format)

## FASTQ

FASTQ or the file format used to store data resulting from sequencing. The sequences present in FASTQ files represent short 
genomic sequences of DNA. FASTQ files are used to build assemblies, MAG binnings, genomic coverage etc. You can provide 
both paired end and unpaired reads to Aviary, as well as short and long reads from a variety of different sequencing platforms.
The file extension for FASTQ files is generally `.fastq`, but often they have been compressed so the extension ends in `.gz`.
Compressed FASTQ files are accepted as input to Aviary so you do not have to uncompress them.

For more info refer to the [wikipedia article](https://en.wikipedia.org/wiki/FASTQ_format)

## BAM/SAM

BAM and SAM (Sequence Alignment/Map) format files are the standard format for indicating the alignment start, end, and quality
of FASTQ files to FASTA files. BAM files are the binary format of SAM files, as such can not be read by conventional means.
When performing read mapping the output from the alignment tool will most likely be in SAM/ BAM files. 

For more info refer to the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
