---
title: Aviary
---

<div class="aviary-hero" markdown>
  <svg class="aviary-hero__motif" viewBox="0 0 900 460" preserveAspectRatio="xMidYMid slice" aria-hidden="true">
    <path d="M -40 380 C 160 380, 220 300, 340 260 C 460 220, 520 140, 700 90" />
    <path d="M -40 420 C 200 430, 280 360, 420 320 C 560 280, 620 200, 820 150" />
    <path d="M 340 260 C 380 200, 440 190, 480 130" />
    <path d="M 420 320 C 470 300, 500 250, 560 220" />
    <circle cx="700" cy="90" r="5"></circle>
    <circle cx="820" cy="150" r="5"></circle>
    <circle cx="480" cy="130" r="4"></circle>
    <circle cx="560" cy="220" r="4"></circle>
    <circle cx="340" cy="260" r="4"></circle>
    <circle cx="420" cy="320" r="4"></circle>
  </svg>

<div class="aviary-hero__inner" markdown>
<div markdown>
<span class="aviary-hero__eyebrow">Metagenomics · Snakemake workflow</span>

# Aviary

Aviary assembles metagenomes, recovers metagenome-assembled genomes (MAGs),
and coordinates their annotation in reproducible local or HPC Snakemake
workflows. It supports short-read, long-read and hybrid analyses.

<div class="aviary-hero__actions" markdown>
[Install Aviary](installation.md){ .md-button .md-button--primary }
[Quickstart](getting-started/quickstart.md){ .md-button }
[Read the PDF manual](https://snh-star.github.io/aviary/pdf/aviary-manual.pdf){ .md-button }
</div>
</div>

<div class="aviary-terminal aviary-reveal">
<div class="aviary-terminal__bar"><span></span><span></span><span></span></div>
<pre class="aviary-terminal__body"><span class="aviary-terminal__prompt">$</span><code class="aviary-terminal__type" data-commands="aviary complete, aviary assemble, aviary recover, aviary annotate, aviary cluster, aviary isolate">aviary complete</code></pre>
</div>
</div>
</div>

## A complete analysis in one command

```bash
aviary complete \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  --output sample_aviary \
  --max-threads 8 \
  --n-cores 8
```

See the [first complete analysis](getting-started/first-analysis.md) walkthrough
for what each stage does and where its outputs land.

## The commands

Six analysis subcommands cover the workflow. `complete` runs the metagenome
stages end to end; the others let you enter or leave at a specific point.
Two utility commands, [`configure`](usage/configure.md) and
[`build`](usage/build.md), manage reference paths and dependency environments.

<div class="aviary-commands" markdown>

<div class="aviary-command" markdown>
<div class="aviary-command__head">
<a class="aviary-command__name" href="usage/complete/"><code>complete</code></a>
<span class="aviary-command__flow">reads <span aria-hidden="true">→</span> annotated MAGs</span>
</div>
<p class="aviary-command__body">Every stage in one command: assembly, binning, refinement, then annotation. Works from short reads, long reads, or both. Already have an assembly? Pass <code>--assembly</code> and Aviary picks up from binning instead of starting over.</p>
</div>

<div class="aviary-command" markdown>
<div class="aviary-command__head">
<a class="aviary-command__name" href="usage/assemble/"><code>assemble</code></a>
<span class="aviary-command__flow">reads <span aria-hidden="true">→</span> contigs</span>
</div>
<p class="aviary-command__body">Step-down hybrid assembly using long and short reads together, playing each to its strengths — long reads for contiguity, short reads for accuracy. Either type alone works too. Several short-read samples with no long reads are co-assembled with MEGAHIT or metaSPAdes.</p>
</div>

<div class="aviary-command" markdown>
<div class="aviary-command__head">
<a class="aviary-command__name" href="usage/recover/"><code>recover</code></a>
<span class="aviary-command__flow">assembly <span aria-hidden="true">→</span> MAGs</span>
</div>
<p class="aviary-command__body">Sorts contigs into metagenome-assembled genomes using several binning algorithms rather than trusting any single one, then refines the result and reports quality and taxonomy per MAG. With no assembly it runs the assembly pipeline first; with several it enables SemiBin2 multi-sample binning.</p>
</div>

<div class="aviary-command" markdown>
<div class="aviary-command__head">
<a class="aviary-command__name" href="usage/annotate/"><code>annotate</code></a>
<span class="aviary-command__flow">MAGs <span aria-hidden="true">→</span> annotations</span>
</div>
<p class="aviary-command__body">Adds functional annotations with EggNOG and taxonomic classifications with GTDB-Tk. Assemblies can be passed alongside when targeting QUAST QC. CheckM2 quality assessment belongs to the <code>recover</code>/<code>complete</code> binning path, not the default annotation target.</p>
</div>

<div class="aviary-command" markdown>
<div class="aviary-command__head">
<a class="aviary-command__name" href="usage/cluster/"><code>cluster</code></a>
<span class="aviary-command__flow">many runs <span aria-hidden="true">→</span> representatives</span>
</div>
<p class="aviary-command__body">Run enough samples and the same organism turns up repeatedly. Galah collapses those near-duplicates across finished Aviary runs and picks one representative genome per cluster — 97% ANI by default, adjustable with <code>--ani</code>.</p>
</div>

<div class="aviary-command" markdown>
<div class="aviary-command__head">
<a class="aviary-command__name" href="usage/isolate/"><code>isolate</code></a>
<span class="aviary-command__flow">pure culture <span aria-hidden="true">→</span> genome</span>
</div>
<p class="aviary-command__body">A long-read-first Flye assembly and polishing path for a single organism from pure culture. Optional short reads add another polishing round. For metagenomic data, use <code>assemble</code> or <code>recover</code> instead.</p>
</div>

</div>

## Find what you need

<div class="grid cards" markdown>

- :material-compass-outline:{ .lg .middle } **[Guide](guide/README.md)**

    ---

    How reads become assemblies, candidate genomes and annotated MAGs.

- :material-book-open-page-variant-outline:{ .lg .middle } **[Worked analyses](examples.md)**

    ---

    End-to-end examples with real commands.

- :material-console:{ .lg .middle } **[CLI reference](usage/README.md)**

    ---

    Every command, option, default and output.

- :material-lifebuoy:{ .lg .middle } **[Troubleshooting](faqs.md)**

    ---

    Validation messages, failed rules, resource limits.

</div>

Using Aviary in research? See the [citation guide](citations.md) for Aviary
and its upstream tools. Aviary is released under GPL-3.0.
