---
title: myloasm k-mer counting (KMC)
---

# myloasm k-mer counting (`--myloasm-kmc`)

By default myloasm counts k-mers in memory with a bloom filter. With
`--myloasm-kmc`, Aviary counts them on disk first, using
[`myloasm-kmc`](https://github.com/bluenote-1577/myloasm-kmc), and hands the
result to myloasm. It is opt-in, applies to myloasm only (myloasm 0.7.0 or
later), and is off unless you ask for it.

```bash
aviary assemble --longreads reads.fastq.gz --long-read-type ont \
  --myloasm-kmc --myloasm-kmc-ram 32 
```

## When to use it

On the two samples benchmarked below, `--myloasm-kmc` **made assembly faster and
gave an essentially identical assembly, but did not lower peak memory**. Use it
when:

- you want shorter runtimes and have scratch disk to spare, or
- the in-memory counting stage is what limits your job.

Skip it when scratch space is tight: counting writes temporary files comparable
to the size of the long reads.

## How it runs

Counting is its own Snakemake step (`count_kmers_kmc`) that runs before
assembly. This keeps two things separate from the assembly job:

- its scratch space, which follows `--tmpdir`;
- its memory budget (`--myloasm-kmc-ram`, in GB), which does not inflate the
  assembly job's memory request.

myloasm then loads the finished database with its own `--kmc-stranded-db`
option, so the in-memory counting stage is skipped.

`myloasm-kmc` is a small Rust wrapper around a fork of KMC 3, shipped as a
separate program because KMC is GPL-3 while myloasm is MIT/Apache. Its README
describes it as written with AI assistance ("user beware"), so check
assemblies you depend on. If you use it, cite KMC 3: Kokot, Długosz and
Deorowicz, *Bioinformatics* 33(17), 2759-2761 (2017).

## Benchmark

Two long-read ONT metagenomes (`aale`: 34.3 GB of compressed reads, `damh`:
34.7 GB), assembled once each in both modes with myloasm 0.7.0, 16 threads and
180 GB of memory. Counting used a 32 GB `--myloasm-kmc-ram` budget. Each sample
has one normal and one KMC run, with no replicates. The finished k-mer
database was about 6 GB on disk for each sample (5.9 GB for damh, 6.5 GB for
aale), kept alongside the assembly. The scratch space used while counting was
not measured.

### Time

| | aale normal | aale KMC | damh normal | damh KMC |
| --- | ---: | ---: | ---: | ---: |
| k-mer counting | 1129 s | 22 s + 316 s in KMC | 851 s | 21 s + 321 s in KMC |
| Whole assembly | 7086 s | 4669 s | 4585 s | 3859 s |
| Saved | | 34% | | 16% |

The counting stage is where KMC reliably wins: about 790 s on aale and 510 s on
damh once the separate KMC step is paid for. The rest of aale's saving came from
stages KMC does not touch (contained-read removal, overlaps), so it is probably
run-to-run or node variation. Treat the damh figure as the conservative one.

### Memory

Peak memory in GB, as logged by myloasm at each stage:

| Stage | aale normal | aale KMC | damh normal | damh KMC |
| --- | ---: | ---: | ---: | ---: |
| **k-mer counting** | **62.9** | **9.3** | **45.1** | **8.5** |
| Mapping reads to outer reads, round 1 | 78.2 | 78.2 | 79.2 | 79.3 |
| Getting overlaps | 26.8 | 29.6 | 31.7 | 31.8 |
| Progressive coverage filtering | 40.3 | 43.4 | 36.5 | 36.5 |
| **Overall peak** | **78.2** | **78.2** | **79.2** | **79.3** |

KMC cuts the counting stage by 36-54 GB (the counting itself runs in the
separate KMC step, capped by `--myloasm-kmc-ram`). The overall peak is set later,
when reads are mapped to the outer reads, and KMC does not change that stage. So
the whole-process peak was the same in both modes (84-86 GB measured with GNU
`time`).

### Assembly quality

| | aale normal | aale KMC | damh normal | damh KMC |
| --- | ---: | ---: | ---: | ---: |
| Contigs | 61,351 | 61,353 | 36,632 | 36,632 |
| Total size | 1358.10 Mb | 1358.11 Mb | 963.73 Mb | 963.72 Mb |
| N50 | 26,897 | 26,897 | 39,098 | 39,098 |
| Largest contig | 4,906,750 | 4,906,727 | 4,159,873 | 4,159,871 |
| Contigs ≥100 kb | 299 | 299 | 475 | 475 |
| Contigs ≥1 Mb | 2 | 2 | 3 | 3 |

The two assemblies are not byte-identical. No second normal run was made, so
how much of the difference is KMC and how much is ordinary run-to-run variation
is not known. Comparing normal against KMC:

| | aale | damh |
| --- | ---: | ---: |
| Contigs with identical sequence | 84.7% | 88.5% |
| 60 bp lines absent from the other assembly | 15.11% | 13.57% |
| Bases aligned at ≥99.9% identity | 99.52% | 99.69% |
| Mean identity of aligned bases | 99.9937% | 99.9952% |

The 60 bp line figure is a strict upper bound: one insertion or deletion shifts
the line frame, so every later line on that contig counts as different. The
base identity is the better measure of how different the sequences really are
(about 0.006% of bases).

