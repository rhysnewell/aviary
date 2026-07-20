---
title: Aviary Skill File
---

# Aviary

**Comprehensive Snakemake-based metagenomics pipeline: assembly → binning → annotation → dereplication.**
Supports short-read, long-read, and hybrid workflows. Each module runs independently or as a chained pipeline.

> **GitHub:** https://github.com/rhysnewell/aviary
> **Docs:** https://rhysnewell.github.io/aviary
> **Latest version:** 0.13.0


---

## 1. INSTALLATION

```bash
# Set conda channels first (required)
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Option A — Bioconda
conda create -n aviary -c bioconda aviary
# OR into an existing environment:
conda install -c bioconda aviary

# Option B — Pip (create environment from admin/requirements.txt first)
pip install aviary-genome

# Option C — Source with pixi (recommended for development)
git clone https://github.com/rhysnewell/aviary.git
cd aviary
pixi run postinstall     # installs all dependencies

# Editable install for development
pip install -e .
```


---

## 2. QUICK REFERENCE — SUBCOMMANDS

```bash
aviary --version            # show installed version
aviary --help               # top-level help
aviary <subcommand> --help  # help for a specific subcommand
```

| Subcommand   | Description                                                                     |
|--------------|---------------------------------------------------------------------------------|
| `assemble`   | QC + de novo assembly (short, long, or hybrid reads)                           |
| `recover`    | MAG recovery from an assembly using multiple binning algorithms                 |
| `annotate`   | Functional (EggNOG) + taxonomic (GTDB-tk) annotation of MAGs                   |
| `complete`   | Runs assemble → recover → annotate in sequence                                  |
| `cluster`    | Combines and dereplicates the MAGs from multiple Aviary runs (Galah-based)      |
| `isolate`    | Isolate (single-organism) hybrid assembly — separate from the meta pipeline     |
| `configure`  | Set database paths and environment variables                                    |


---

## 3. REQUIRED DATABASES

| Database  | Tool              | Typical size | Required?                            |
|-----------|-------------------|-------------|--------------------------------------|
| GTDB      | GTDB-tk (taxonomy) | ~85 GB     | Yes (for taxonomy in recover/annotate)|
| EggNOG    | EggNOG-mapper      | ~50 GB     | Yes (for functional annotation)       |
| CheckM2   | CheckM2 (quality)  | ~3 GB      | Yes (for MAG QC)                      |
| SingleM   | SingleM (profiling)| ~14 GB     | Yes (for recovery assessment)         |
| Metabuli  | TaxVAMB binner     | varies     | Optional (only for taxvamb binner)    |

Use `aviary configure --download` to auto-install all databases (see Section 4).

### Environment variable paths (set in `~/.bashrc` or conda activate script)

```bash
export GTDBTK_DATA_PATH=/path/to/gtdb/db/
export EGGNOG_DATA_DIR=/path/to/eggnog/
export SINGLEM_METAPACKAGE_PATH=/path/to/singlem.smpkg/
export CHECKM2DB=/path/to/checkm2db/
export CONDA_ENV_PATH=/path/to/conda/envs/
```


---

## 4. CONFIGURE — set database paths and optionally download

```bash
aviary configure \
  -o logs/ \
  --eggnog-db-path /path/to/eggnog/ \
  --gtdb-path /path/to/gtdb/ \
  --checkm2-db-path /path/to/checkm2/ \
  --singlem-metapackage-path /path/to/singlem/ \
  --metabuli-db-path /path/to/metabuli/ \    # optional, for taxvamb
  --download                                  # auto-download missing databases

# Also set a persistent tmpdir via configure
aviary configure --tmpdir /scratch/tmp/
```

| Flag                          | Description                                         |
|-------------------------------|-----------------------------------------------------|
| `-o`                          | Output / log directory                              |
| `--eggnog-db-path`            | EggNOG functional annotation database path          |
| `--gtdb-path`                 | GTDB taxonomy database path                         |
| `--checkm2-db-path`           | CheckM2 quality assessment database path            |
| `--singlem-metapackage-path`  | SingleM marker gene database path                   |
| `--metabuli-db-path`          | Metabuli database (optional; required for taxvamb)  |
| `--download`                  | Auto-download all missing databases                 |
| `--tmpdir`                    | Persist a custom TMPDIR for all future runs         |


---

## 5. ASSEMBLE — de novo assembly from sequencing reads

QC runs automatically before assembly unless `--skip-qc` is set.

```bash
# Short-read only (metaSPAdes by default)
aviary assemble \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o output_dir/ \
  --n-cores 16

# Short-read only with MEGAHIT (faster, lower memory than metaSPAdes)
aviary assemble \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o output_dir/ \
  --n-cores 16 \
  --use-megahit

# Long-read only (myloasm assembler, default)
aviary assemble \
  -l sample.fastq.gz \
  -z ont \
  -o output_dir/ \
  --n-cores 16

# Long-read only with Flye
aviary assemble \
  -l sample.fastq.gz \
  -z ont \
  --long-read-assembler flye \
  -o output_dir/ \
  --n-cores 16

# Hybrid (short + long reads)
aviary assemble \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -l sample.fastq.gz \
  -z ont \
  -o output_dir/ \
  --n-cores 16

# Host contamination removal before assembly
aviary assemble \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -r host_genome.fasta \
  -o output_dir/ \
  --n-cores 16

# Co-assembly of multiple samples
aviary assemble \
  -1 s1_R1.fq.gz s2_R1.fq.gz \
  -2 s1_R2.fq.gz s2_R2.fq.gz \
  --coassemble \
  -o output_dir/ \
  --n-cores 16
```

### Read input flags

| Flag                          | Description                                               |
|-------------------------------|-----------------------------------------------------------|
| `-1 PATH ..`                  | Forward short reads                                       |
| `-2 PATH ..`                  | Reverse short reads                                       |
| `-i / --interleaved PATH ..`  | Interleaved short read files                              |
| `-c / --coupled PATH ..`      | Coupled forward/reverse list                              |
| `-l / --longreads PATH ..`    | Nanopore or PacBio reads                                  |
| `-z / --longread-type TYPE`   | Read type: `ont` \| `ont_hq` \| `rs` \| `sq` \| `ccs` \| `hifi` (default: `ont`) |

### Assembly flags

| Flag                       | Description                                                        | Default    |
|----------------------------|--------------------------------------------------------------------|------------|
| `--long-read-assembler`    | `myloasm` \| `flye`                                                | `myloasm`  |
| `--use-megahit`            | Use MEGAHIT instead of metaSPAdes for short reads                  | off        |
| `--use-unicycler`          | Use Unicycler for metaSPAdes reassembly                            | off        |
| `--coassemble`             | Co-assemble multiple input read sets                               | off        |
| `-k / --kmer-sizes`        | Kmer sizes for SPAdes (odd integers <128)                          | auto       |
| `--medaka-model`           | Medaka model for long-read polishing                               | `r941_min_hac_g507` |

### Contig filtering flags (applied post-assembly)

| Flag                     | Description                                              | Default  |
|--------------------------|----------------------------------------------------------|----------|
| `--min-cov-long`         | Minimum long read coverage for contigs                   | `5`      |
| `--min-cov-short`        | Minimum short read coverage for contigs                  | `5`      |
| `--exclude-contig-cov`   | Exclude contigs with coverage ≤ this value               | `10`     |
| `--exclude-contig-size`  | Exclude contigs shorter than this bp                     | `2500`   |
| `--include-contig-size`  | Always include contigs longer than this bp               | `10000`  |

### Expected outputs (`output_dir/`)

```
assembly/
  final_contigs.fasta        # assembled contigs (use this as --assembly for recover)
www/
  assembly_stats.txt         # assembly summary statistics
  rastqc/                    # RastQC short-read QC reports
  nanoplot/                  # NanoPlot long-read QC reports
```


---

## 6. QC / READ FILTERING FLAGS (apply to assemble, recover, complete)

QC runs automatically. Use these flags to tune or skip it.

| Flag                           | Description                                           | Default |
|--------------------------------|-------------------------------------------------------|---------|
| `--skip-qc`                    | Skip all QC steps                                     | off     |
| `-r / --host-filter PATH ..`   | Host reference fasta(s) to remove host reads          | —       |
| `--min-read-size INT`          | Minimum long read size for Filtlong                   | `100`   |
| `--min-mean-q FLOAT`           | Minimum long read mean quality                        | `10`    |
| `--min-short-read-length INT`  | Minimum short read length to keep                     | `15`    |
| `--max-short-read-length INT`  | Maximum short read length (0 = no max)                | `0`     |
| `--disable-adapter-trimming`   | Disable fastp adapter trimming                        | off     |
| `--unqualified-percent-limit FLOAT` | Max % bases allowed unqualified                  | `40`    |
| `--quality-cutoff INT`         | Short read quality threshold                          | `15`    |
| `--extra-fastp-params STR`     | Extra params passed to fastp as string                | —       |


---

## 7. RECOVER — extract MAGs from an assembly using multiple binners

### Default binners (always run unless skipped)

| Binner     | Type        | GPU? | Notes                                              |
|------------|-------------|------|----------------------------------------------------|
| Rosella    | short+long  | No   | CoverM-based differential coverage binner (default lead binner) |
| MetaBAT2   | short       | No   | Coverage-based                                     |
| SemiBin2   | short+long  | Yes  | Deep learning, self-supervised                     |
| VAMB       | short       | No   | Variational autoencoder coverage binner            |

### DAS_Tool is always run to produce a final consensus bin set from all binner outputs.

```bash
# Short reads only
aviary recover \
  --assembly assembly.fasta \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o output_dir/ \
  --n-cores 16 \
  --max-threads 32

# With long reads for coverage (improves binning)
aviary recover \
  --assembly assembly.fasta \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -l sample.fastq.gz \
  -z ont \
  -o output_dir/ \
  --n-cores 16

# GPU-accelerated binners (taxvamb, comebin, semibin use CUDA)
aviary recover \
  --assembly assembly.fasta \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o output_dir/ \
  --n-cores 16 \
  --max-threads 32 \
  --request-gpu

# Add extra binners beyond defaults
aviary recover \
  --assembly assembly.fasta \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o output_dir/ \
  --n-cores 16 \
  --extra-binners taxvamb comebin

# Skip specific binners
aviary recover \
  --assembly assembly.fasta \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o output_dir/ \
  --n-cores 16 \
  s vamb semibin

# Binning only (skip annotation — useful for quick MAG recovery)
aviary recover \
  --assembly assembly.fasta \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o output_dir/ \
  --n-cores 16 \
  --binning-only

# Target a specific snakemake rule (run only up to this step)
aviary recover \
  --assembly assembly.fasta \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o output_dir/ \
  -w rosella \
  --max-threads 12
```

### Recover flags

| Flag                              | Description                                                         | Default    |
|-----------------------------------|---------------------------------------------------------------------|------------|
| `--assembly PATH`                 | Input assembled contigs (.fasta) — required                         | —          |
| `-1 / -2`                         | Short reads for coverage calculation                                | —          |
| `-l / --longreads`                | Long reads for coverage calculation                                 | —          |
| `-z / --longread-type`            | Long read type (see Section 5)                                      | `ont`      |
| `--request-gpu`                   | Enable GPU-accelerated binners (taxvamb, comebin, semibin)          | off        |
| `--binning-only`                  | Stop after binning; skip quality check and annotation               | off        |
| `--strict`                        | Fail immediately if any binner errors                               | off        |
| `--extra-binners LIST`            | Add extra binners: `maxbin` \| `maxbin2` \| `concoct` \| `comebin` \| `taxvamb` \| `quickbin` | — |
| `--skip-binners LIST`             | Skip default binners: `rosella` \| `semibin` \| `metabat1` \| `metabat2` \| `metabat` \| `vamb` | — |
| `--min-contig-size INT`           | Minimum contig size for binning (bp)                                | `1500`     |
| `--min-bin-size INT`              | Minimum MAG size (bp)                                               | `200000`   |
| `--min-completeness FLOAT`        | Minimum CheckM2 completeness % for bins passed to annotation        | `50.0`     |
| `--max-contamination FLOAT`       | Maximum CheckM2 contamination % for bins passed to annotation       | `5.0`      |
| `--semibin-model NAME`            | SemiBin2 environment model                                          | `global`   |
| `--refinery-max-iterations INT`   | Rosella refinery max iterations                                     | `5`        |
| `--refinery-max-retries INT`      | Rosella refinery max retries                                        | `3`        |
| `--coverage-job-strategy STR`     | Strategy for coverage jobs: `default` \| `never` \| `always`       | `default`  |
| `--coverage-samples-per-job INT`  | Samples per coverage job                                            | `5`        |
| `--skip-abundances`               | Skip CoverM post-binning abundance calculations                     | off        |
| `--skip-taxonomy`                 | Skip GTDB-tk taxonomy assignment                                    | off        |
| `--skip-singlem`                  | Skip SingleM recovery assessment                                    | **on**     |
| `--min-percent-read-identity-short FLOAT` | Min % identity for short-read CoverM                      | `95`       |
| `--min-percent-read-identity-long FLOAT`  | Min % identity for long-read CoverM                       | `85`       |
| `-w / --workflow RULE`            | Target a specific Snakemake rule (runs all prerequisites first)     | —          |

### Expected outputs (`output_dir/`)

```
bins/
  final_bins/                  # final consensus MAG FASTA files (use these downstream)
  bin_info.tsv                 # taxonomy, QC (completeness/contamination), size, N50 for all MAGs
  checkm_minimal.tsv           # minimal CheckM2 summary
  coverm_abundances.tsv        # CoverM relative abundance of each MAG per sample
diversity/
  singlem_out/                 # SingleM read-based community profiling output
taxonomy/                      # GTDB-tk raw output
benchmarks/                    # Snakemake rule benchmarks (time, memory per rule)
logs/                          # Rule-level logs — check here first on errors
```


---

## 8. ANNOTATE — functional (EggNOG) + taxonomic (GTDB-tk) annotation

```bash
aviary annotate \
  -d recover/bins/final_bins/ \
  -o output_dir/ \
  --n-cores 16
```

| Flag                            | Description                                           | Default |
|---------------------------------|-------------------------------------------------------|---------|
| `-d / --genome-fasta-directory` | Directory containing MAG FASTA files                  | —       |
| `-x / --fasta-extension`        | File extension of MAG files (**default: `fna`**, not `.fasta`) | `fna` |

> **Note:** The default extension is `.fna`. If your bins have a different extension (e.g. `.fa`,
> `.fasta`), always set `-x` accordingly.

### Expected outputs (`output_dir/`)

```
annotation/
  eggnog/                      # EggNOG-mapper functional annotation output
  gtdbtk/                      # GTDB-tk taxonomic classification output
```


---

## 9. COMPLETE — assemble + recover + annotate in one command

`complete` does NOT include `cluster` (dereplication). Run that separately after all samples.

```bash
# Short reads only
aviary complete \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o output_dir/ \
  --n-cores 16

# Hybrid with myloasm (default)
aviary complete \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -l sample.fastq.gz \
  -z ont \
  -o output_dir/ \
  --n-cores 16

# Hybrid with Flye
aviary complete \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -l sample.fastq.gz \
  -z ont \
  --long-read-assembler flye \
  -o output_dir/ \
  --n-cores 16

# With GPU binners
aviary complete \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -l sample.fastq.gz \
  -z ont \
  -o output_dir/ \
  --n-cores 16 \
  --request-gpu

# Preview the workflow without running
aviary complete \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o output_dir/ \
  --dry-run
```


---

## 10. CLUSTER — dereplicate MAGs across multiple Aviary runs

Run this **once, after all samples are processed**, to dereplicate the combined MAG set.

```bash
aviary cluster \
  -i run1/ run2/ run3/ \
  -o cluster_output/ \
  --n-cores 16
```

| Flag                      | Description                                            | Default  |
|---------------------------|--------------------------------------------------------|----------|
| `-i / --input-runs PATH..`| Space-separated previous Aviary output directories     | —        |
| `--ani FLOAT`             | ANI threshold for dereplication                        | `97`     |
| `--precluster-ani FLOAT`  | Minimum ANI for preclustering                          | `95`     |
| `--precluster-method NAME`| Rough ANI method: `dashing` \| `finch`                 | `dashing`|
| `--min-completeness FLOAT`| Exclude genomes below this completeness %              | —        |
| `--max-contamination FLOAT` | Exclude genomes above this contamination %           | —        |
| `--use-checkm2-scores`    | Use CheckM2 scores for representative selection        | off      |
| `--pggb-params STR`       | Parameters passed to pggb (surround with quotes)       | `-k 79 -G 7919,8069` |

### Expected outputs (`cluster_output/`)

```
dereplicated_bins/             # final dereplicated representative MAG FASTA files
```


---

## 11. ISOLATE — isolate (single-organism) assembly

Separate workflow from metagenomics. For isolated pure culture sequencing (e.g. bacterial isolates with Nanopore + Illumina).

### Pipeline steps

1. **Flye** — assembles long reads → `isolate/flye/assembly.fasta`
2. **Racon** (4 rounds, long reads) — polishes the Flye assembly → `isolate/isolate.pol.rac.fasta`
3. **Medaka** — further long-read polishing → `isolate/isolate.pol.med.fasta`
4. **Pilon** (short reads, if provided) — short-read polishing → `isolate/isolate.pol.pil.fasta`
5. **Racon** (1 round, short reads, if provided) — final short-read polish → `isolate/isolate.pol.fin.fasta`
6. **dnaapler** — reorients circular chromosomes to a standard start position → `isolate/completed_assembly.fasta`

Short reads are optional. If provided, Pilon and illumina Racon polishing steps are added after Medaka.

### Final output

```
isolate/completed_assembly.fasta    ← use this as the final assembly
```

```bash
# Hybrid (long + short reads) — recommended
aviary isolate \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -l sample.fastq.gz \
  -z ont \
  -o output_dir/ \
  --n-cores 16

# Long reads only
aviary isolate \
  -l sample.fastq.gz \
  -z ont \
  -o output_dir/ \
  --n-cores 16

# Run only up to dnaapler (reorientation) step
aviary isolate \
  -l sample.fastq.gz \
  -z ont \
  -o output_dir/ \
  --n-cores 16 \
  -w dnaapler
```

| Flag              | Description                                          | Default     |
|-------------------|------------------------------------------------------|-------------|
| `--guppy-model`   | Medaka model for polishing                           | `r941_min_hac_g507` |
| `--genome-size`   | Approximate genome size in bp (informational)        | `5000000`   |


---

## 12. COMMON FLAGS (apply across most subcommands)

| Flag                            | Description                                                       | Default |
|---------------------------------|-------------------------------------------------------------------|---------|
| `-n / --n-cores INT`            | Max CPU cores (total for Snakemake)                               | `16`    |
| `--max-threads INT`             | Max threads for any single process                                | `8`     |
| `-m / --max-memory INT`         | Max memory in GB                                                  | `250`   |
| `-p / --pplacer-threads INT`    | Threads for pplacer (dedicated — avoid high values to prevent deadlocks) | `8` |
| `-o / --output PATH`            | Output directory                                                  | `./`    |
| `--request-gpu`                 | Enable GPU-accelerated steps (taxvamb, comebin, semibin)          | off     |
| `--dry-run`                     | Show workflow steps without executing                             | off     |
| `-w / --workflow RULE`          | Target a specific Snakemake rule (all prerequisites still run)    | —       |
| `--snakemake-cmds STR`          | Extra arguments passed directly to Snakemake                      | —       |
| `--snakemake-profile NAME/PATH` | Snakemake cluster profile for HPC submission                      | —       |
| `--local-cores INT`             | Max cores for local-only steps                                    | `16`    |
| `--cluster-retries INT`         | Retry failed cluster jobs N times (memory/time auto-scaled on retry) | `0`  |
| `--tmpdir PATH`                 | Directory for temporary files                                     | `/tmp`  |
| `--clean`                       | Clean up temp files after run                                     | `True`  |
| `--rerun-triggers VALUE`        | What triggers rule reruns                                         | `mtime` |
| `--default-resources STR`       | Snakemake default resource config string                          | —       |

### Thread control explained

- `-t / --max-threads` — threads **per program** (how many each tool gets)
- `-n / --n-cores` — total cores Snakemake can use (if > `--max-threads`, multiple rules run in parallel)
- `-p / --pplacer-threads` — pplacer gets its own limit (deadlocks on too many threads)

### Unlock a directory after a crash

```bash
aviary complete \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o output_dir/ \
  --snakemake-cmds '--unlock'
```

### Resume an interrupted run

Snakemake automatically skips completed steps. Just re-run the exact same command.


---

## 13. HPC / CLUSTER JOB SUBMISSION

### PBS/qsub

```bash
# Inline qsub submission
aviary complete \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o output_dir/ \
  --n-cores 16 \
  --snakemake-cmds '--cluster "qsub -V" '
  # NOTE: space after the closing quote is required (argparse quirk)

# Using a named Snakemake profile (e.g. ~/.config/snakemake/aqua/config.yaml)
aviary complete \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o output_dir/ \
  --n-cores 16 \
  --snakemake-profile aqua \
  --cluster-retries 2
```

### Example PBS `~/.config/snakemake/aqua/config.yaml`

```yaml
cluster: qsub
cluster-status: qstat
jobs: 10000
cluster-cancel: qdel
```

> **PBS quirk:** Add `-V` to your qsub command to pass environment variables to jobs.
> Without `-V`, pysam fails to activate properly in submitted jobs.

### SLURM

```bash
aviary complete \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o output_dir/ \
  --n-cores 16 \
  --snakemake-profile /path/to/slurm-profile/
```

### Example SLURM wrapper script (`submit_aviary.sh`)

```bash
#!/bin/bash
#SBATCH --job-name=aviary
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=aviary_%j.log

conda activate aviary
aviary complete \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -o output_dir/ \
  --n-cores 16
```

```bash
sbatch submit_aviary.sh
```

### Cluster retry behaviour

With `--cluster-retries N`, failed jobs are automatically resubmitted up to N times.
Each retry increases the requested memory and wall time. Set `--max-memory` and
`--max-threads` to cap these values.

> **Note on job resources:** Job CPU and memory requirements were empirically derived from
> 1,000 Aviary runs. Binner rules are standardised at 24 CPUs. Memory is set to the nearest
> power of 2 above max RSS. These are capped by `--max-memory` / `--max-threads`.


---

## 14. FULL WORKFLOW WALKTHROUGH (end-to-end)

### Step-by-step (maximum control)

```bash
# Step 1 — Configure databases (once per system)
aviary configure \
  -o logs/ \
  --eggnog-db-path /path/to/db/eggnog/ \
  --gtdb-path /path/to/db/gtdb/ \
  --checkm2-db-path /path/to/db/checkm2/ \
  --singlem-metapackage-path /path/to/db/singlem/ \
  --download

# Step 2 — QC + Assembly
aviary assemble \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -l sample.fastq.gz \
  -z ont \
  -o 01_assemble/ \
  --n-cores 16

# Step 3 — MAG Recovery
aviary recover \
  --assembly 01_assemble/assembly/final_contigs.fasta \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -l sample.fastq.gz \
  -z ont \
  -o 02_recover/ \
  --n-cores 16 \
  --max-threads 32 \
  --request-gpu       # only if GPU node available

# Step 4 — Annotate MAGs
aviary annotate \
  -d 02_recover/bins/final_bins/ \
  -o 03_annotate/ \
  --n-cores 16

# Step 5 — Dereplicate across all samples (run once all samples processed)
aviary cluster \
  -i 02_recover/ run2/ run3/ \
  -o 04_cluster/ \
  --n-cores 16
```

### Single-command alternative (steps 2–4 only)

```bash
aviary complete \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -l sample.fastq.gz \
  -z ont \
  -o output_dir/ \
  --n-cores 16 \
  --request-gpu

# Then run cluster (step 5) separately after all samples are done
aviary cluster -i output_dir/ run2/ run3/ -o cluster/ --n-cores 16
```


---

## 15. KEY OUTPUT FILES

```
# After assemble:
output_dir/assembly/final_contigs.fasta     ← use as --assembly in recover

# After recover:
output_dir/bins/final_bins/*.fna            ← final MAG FASTA files (use these downstream)
output_dir/bins/bin_info.tsv                ← taxonomy, completeness, contamination, size, N50
output_dir/bins/checkm_minimal.tsv          ← minimal CheckM2 QC summary
output_dir/bins/coverm_abundances.tsv       ← relative abundance per MAG per sample
output_dir/diversity/singlem_out/           ← SingleM community profiling
output_dir/taxonomy/                        ← raw GTDB-tk output
output_dir/benchmarks/                      ← per-rule resource benchmarks
output_dir/logs/                            ← per-rule logs — always check here first on errors

# After annotate:
output_dir/annotation/eggnog/               ← EggNOG-mapper functional annotation
output_dir/annotation/gtdbtk/              ← GTDB-tk taxonomic annotation

# After cluster:
output_dir/dereplicated_bins/               ← final dereplicated representative MAGs
```


---

## 16. WEB INTERFACE

A live log viewer for monitoring running Aviary jobs in the browser. Parses the Snakemake log in real time and displays rule progress, resource usage, errors, assembly stats, and (if available) SingleM phylo tree and GTDB-tk results.

Requires the `web` pixi environment (Flask).

```bash
# Recommended — via pixi task
pixi run -e web server

# Or pin to a specific output directory
python -m aviary.web.server --output-dir /path/to/aviary/output
```

| Flag            | Description                                              | Default       |
|-----------------|----------------------------------------------------------|---------------|
| `--output-dir`  | Pin to a specific Aviary output directory                | cwd           |
| `--port`        | Port to listen on (auto-increments if in use)            | `8090`        |
| `--host`        | Interface to bind to                                     | `127.0.0.1`   |
| `--reload`      | Auto-reload server on source file changes (dev only)     | off           |

> **Note:** The server binds to `127.0.0.1` (loopback only). On a shared HPC you must use an
> SSH tunnel to view it locally.

### SSH tunnel (HPC)

**Terminal 1** — open the SSH tunnel (keep this running):
```bash
ssh -L 8090:localhost:8090 <username>@<hpc-address>
```

**Terminal 2** — start the server on the HPC:
```bash
python -m aviary.web.server --output-dir /path/to/aviary/output --port 8090
```

**Browser**: open `http://localhost:8090`

If port 8090 is already in use the server automatically tries 8091, 8092, … up to 8109.
Update the `-L` tunnel port to match whichever port the server reports.


---

## 17. TIPS, GOTCHAS & COMMON ERRORS

1. **`/tmp` filling up** is the most common error (`Error in prepare_binning_files`). CoverM
   mapping during binning uses `/tmp`. Fix: `--tmpdir /scratch/tmp/` or
   `aviary configure --tmpdir /scratch/tmp/`.

2. **SPAdes "Error code: -9"** almost always means out of memory. Increase with `-m` (GB).

3. **qsub + pysam ModuleNotFoundError:** Add `-V` to your qsub command to pass environment
   variables: `--snakemake-cmds '--cluster "qsub -V" '`

4. **Note the trailing space** in `--snakemake-cmds '--cluster qsub '`. The space after
   `qsub` is required due to a Python argparse quirk.

5. **`--request-gpu` on a local machine:** First install the `cuda` package into your conda
   environment. GPU-capable programs will detect it automatically.

6. **Resuming a run:** Snakemake automatically skips completed rules. Just re-run the exact
   same command with the same `-o` output directory.

7. **Default fasta extension for `annotate` is `.fna`**, not `.fasta`. If your bins have a
   different extension, always specify `-x`.

8. **`aviary complete` does NOT run `cluster`.** You must always run `aviary cluster` separately
   after all samples have been processed.

9. **`--binning-only`** in recover skips CheckM2 QC and GTDB-tk — useful when you just want
   fast MAG binning without annotation and plan to run `aviary annotate` separately.

10. **`--cluster-retries`** is highly recommended on HPC. Failed jobs are retried with
    progressively more memory and time, up to `--max-memory` / `--max-threads` caps.

11. **Setting GTDB path:** Point `GTDBTK_DATA_PATH` to the `db/` folder *inside* the GTDB
    download directory, not the download root itself.

12. **GPU binners require a GPU node.** If using `--request-gpu` with a PBS/SLURM profile,
    ensure the profile is configured to submit GPU jobs to a GPU partition/queue.

13. **`--skip-singlem` is on by default.** SingleM recovery assessment is skipped unless you
    explicitly pass `--skip-singlem False`.


---

## 18. CITATION

If you use Aviary in your research, please cite:

> Newell RJP, Aroney STN, Zaugg J, Sternes P, Tyson GW, Woodcroft BJ.
> **Aviary: Hybrid assembly and genome recovery from metagenomes with Aviary.**
> Zenodo (2024). https://doi.org/10.5281/zenodo.10806928

```
@software{aviary,
  author  = {Newell, Rhys J. P. and Aroney, Samuel T. N. and Zaugg, Julian and
             Sternes, Peter and Tyson, Gene W. and Woodcroft, Ben J.},
  title   = {Aviary: Hybrid assembly and genome recovery from metagenomes with Aviary},
  year    = {2024},
  doi     = {10.5281/zenodo.10806928},
  url     = {https://doi.org/10.5281/zenodo.10806928}
}
```

You should also cite the individual tools used in your run. A full up-to-date list is at:
https://rhysnewell.github.io/aviary/citations
