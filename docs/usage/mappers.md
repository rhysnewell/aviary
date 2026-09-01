---
title: Read mappers
---

# Read mappers

`recover`, `assemble`, `complete` and `isolate` all map reads back onto an
assembly for coverage/abundance and (where applicable) polishing. This page
documents the mapper selection flags shared across those subcommands — see
[Shared options](centralised_commands.md) for the resource/output/execution
flags, and each subcommand's own page for flags specific to it.

## Choosing a mapper

| Flag | Options |
|---|---|
| `--short-read-mapper` | `strobealign` (default), `minimap2`, `rammap`, `minibwa`, `bwa-mem`, `bwa-mem2`, `strobealign-aemb` |
| `--long-read-mapper` | `rammap` (default), `minimap2`, `minibwa` |

The `-x`/preset for long reads is chosen from `--longread-type` by default, so
`--long-read-mapper` changes only which aligner runs, not which preset it uses.
Note that strobealign, bwa-mem and bwa-mem2 are short-read only, which is why
the two flags are separate.

The defaults are a reasonable choice for almost everyone. The alternatives
exist for continuity and for cases where a particular aligner is preferred.

### My coverage / abundance numbers differ from an older Aviary version

The default read mappers changed. Aviary now uses **strobealign** for short
reads and **rammap** for long reads, where previous versions used minimap2
throughout. Both are faster, and rammap is a minimap2-compatible implementation,
but a different aligner makes different alignment decisions — so coverage
depths, bin abundances and polished contigs will not be bit-identical to those
from an earlier release.

Nothing needs to change to keep using Aviary: every existing command still runs,
and no new flag is required. To reproduce the previous behaviour exactly, ask
for minimap2 explicitly:

```bash
aviary recover --short-read-mapper minimap2 --long-read-mapper minimap2 ...
```

This matters most if you are partway through an analysis, or comparing against
results produced by an earlier version. For new work, the defaults are the
faster option.

## Choosing a preset explicitly with `--*-mapper-model`

**`--short-read-mapper-model`** `{sr,no-preset}`

**`--long-read-mapper-model`** `{lr-hq,ont,pb,hifi,no-preset}`

`minimap2` and `rammap` are the only families with more than one CoverM
preset. Use these flags to pick one directly instead of relying on the
`--longread-type` default:

```bash
aviary recover --long-read-mapper rammap --long-read-mapper-model hifi ...
```

This also fixes a previously dead code path: `ont_hq` and `hifi` used to be
silently mapped with the `ont`/`pb` presets respectively, because CoverM's
`lr-hq` and `hifi` presets were never actually reachable. They are now used
automatically for `--long-read-type ont_hq` / `hifi`, which changes coverage
numbers slightly for those two read types compared to earlier versions of
this flag. Giving `--*-mapper-model` for `strobealign`, `minibwa`, `bwa-mem`,
`bwa-mem2` or `strobealign-aemb` (which have no selectable model) is an error.

Models are also checked against read length: `sr` and `no-preset` are the only
valid `--short-read-mapper-model` values, and `lr-hq`, `ont`, `pb`, `hifi` and
`no-preset` the only valid `--long-read-mapper-model` values. Crossing them
(e.g. `--short-read-mapper-model ont`) is rejected up front, because CoverM
would otherwise accept it and return a well-formed coverage table of near-zero
depths — a wrong number rather than an error.

### `ont` and `pb` are legacy presets, `lr-hq` and `hifi` are the defaults

`--longread-type ont` (and `ont_hq`) now default to the `lr-hq` model
(minimap2/rammap `-x lr:hq`), not `ont` (`-x map-ont`). `--longread-type ccs`
now defaults to `hifi` (`-x map-hifi`), not `pb` (`-x map-pb`). Per
[minimap2's own docs](https://lh3.github.io/minimap2/minimap2.html):

- `map-ont` "align[s] noisy long reads of ~10% error rate" — accurate for
  older ONT chemistry, far too liberal for current chemistry v14 reads (~99%
  accuracy).
- `map-pb` "is effectively deprecated by HiFi... unless you work on very old
  data, you probably want to use `map-hifi` or `lr:hq`."
- `lr:hq` "was recommended by ONT developers for recent Nanopore reads
  produced with chemistry v14... shown to work better for accurate Nanopore
  reads than `map-hifi`."

CCS reads are, by PacBio's current terminology, HiFi reads (CCS at sufficient
pass count *is* HiFi), so routing them through the CLR-era `map-pb` preset was
the same class of problem.

**`ont` (`map-ont`) and `pb` (`map-pb`) remain fully available** — they were
not removed, only demoted from being the default. Use them explicitly for
genuinely noisy/legacy data:

```bash
aviary recover --long-read-mapper minimap2 --long-read-mapper-model ont ...
```

`rs`/`sq` (PacBio RSII/Sequel) are unaffected — they are genuine older CLR
chemistry, so `map-pb` is still the correct preset for them, not a
deprecated-by-mistake one. `hifi` was already on `map-hifi` and is unchanged.

## Mappers are checked against the reads you supply

Naming a mapper for reads you did not provide is an error rather than a
silently ignored flag:

```bash
# error: --short-read-mapper given, but no short reads to map
aviary recover --assembly contigs.fasta --short-read-mapper bwa-mem ...
```

The same applies to `--long-read-mapper` / `--long-read-mapper-model` without
`-l`/`--longreads`, and to `--minibwa-params` when `minibwa` is not the
selected mapper. Defaults (`strobealign` for short reads, `rammap` for long
reads) are unaffected — they only apply to reads that were actually supplied,
so a short-read-only run never trips the long-read check.

## `bwa-mem` / `bwa-mem2`

These are the reference BWA implementations, usable as `--short-read-mapper`
for coverage, abundance and racon polishing. They are the slowest option here
(no inline indexing, and polishing re-indexes the reference every racon
round), but are included for continuity with pipelines built around them.

## `strobealign-aemb`

This is not a normal aligner choice like the others — it is CoverM's
`-m strobealign-aemb` fast direct abundance estimator, which shells out to
`strobealign --aemb` and skips alignment/pileup entirely. It has narrower
scope than the other `--short-read-mapper` values:

- Only used for the per-contig binning coverage step (`data/coverm.cov`).
  Per-genome relative abundance (`bins/coverm_abundances.tsv`) always falls
  back to plain `strobealign` instead, since CoverM cannot run
  `strobealign-aemb` through `coverm genome` at all.
- Faster than a full alignment, but less precise — treat it as a speed/
  precision tradeoff, not a strictly-better default.
- No selectable model: `--short-read-mapper-model` is an error alongside it,
  the same as for `strobealign`, `minibwa`, `bwa-mem` and `bwa-mem2`.
- Short-read only, like `strobealign`/`bwa-mem`/`bwa-mem2`.

## `minibwa` for long reads

`minibwa` has no long-read preset of its own in CoverM. Use
`--minibwa-params` to pass its native preset flag directly, e.g.:

```bash
aviary recover --long-read-mapper minibwa --minibwa-params "-x lr" ...
```

`minibwa` can be used for racon polishing too: aviary calls `minibwa map -f -x
lr` (`map` is minibwa's own subcommand for read alignment; `-x` is one of its
options, not a top-level flag like minimap2's) to get the PAF racon needs.

## Raw per-aligner passthrough params

CoverM exposes a raw params passthrough for every aligner it wraps, and
aviary exposes all of them, mirroring `--minibwa-params`:

| Flag | Applies to |
|---|---|
| `--bwa-params` | `--short-read-mapper bwa-mem` / `bwa-mem2` |
| `--strobealign-params` | `--short-read-mapper strobealign` (not `strobealign-aemb`, which does not accept it) |
| `--minimap2-params` | `--short-read-mapper` or `--long-read-mapper` `minimap2` |
| `--rammap-params` | `--short-read-mapper` or `--long-read-mapper` `rammap` |
| `--minibwa-params` | `--short-read-mapper` or `--long-read-mapper` `minibwa` |

Each is checked against the reads you supply and the mapper you selected the
same way `--minibwa-params` always was — giving one without the matching
mapper selected is an error rather than a silently ignored flag.
