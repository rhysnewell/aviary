---
title: Installation
---

# Installation

Aviary targets Linux and uses isolated dependency environments for its
bioinformatics tools. Three install routes are documented below — most users
only need one:

!!! tip "Which install route?"
    - **Most users:** [Pixi](#install-from-source-with-pixi) — the recommended, supported route.
    - **Already managing tool dependencies yourself:** [pip](#install-from-pip).
    - **Prefer a Conda-style environment:** [Bioconda](#install-from-bioconda).

## Install from source with Pixi

From a local checkout:

```bash
cd aviary
pixi run postinstall
pixi run aviary --version
```

This installs Aviary in editable mode. Dependency definitions live in
`aviary/pixi.toml`, while `aviary/pixi.lock` pins the resolved environments.
The `postinstall` task prepares both the main and development environments.

For a shared development system, database paths can be represented by symlinks
in the repository's `db/` directory; `admin/set_env_vars.sh` documents the
expected local names used by the activation hook.

## Install from pip

The Python package name is `aviary-genome`, but pip alone does not provision
the complete collection of external bioinformatics tools. Use this route only
when you are deliberately managing those dependencies yourself:

```bash
conda env create -n aviary -f admin/environment.yml
conda activate aviary
pip install aviary-genome
```

## Install from Bioconda

Create a dedicated environment:

```bash
conda create -n aviary -c conda-forge -c bioconda aviary
conda activate aviary
```

If your Conda channels are configured globally, keep this priority:

```yaml
channels:
  - conda-forge
  - bioconda
  - defaults
```

Installing into a new environment avoids dependency conflicts with unrelated
analysis software.

## Verify the command

```bash
aviary --version
aviary complete --full-help
```

Before downloading large databases, use `aviary complete --full-help` to confirm
that the installed release exposes the options used by this documentation.

## Build analysis environments

After installation, Aviary can prepare its per-tool environments:

```bash
aviary build
```

GPU environments are optional and require compatible hardware and drivers:

```bash
aviary build --gpu
```

## Configure reference data

Aviary can use local data for GTDB-Tk, EggNOG-mapper, CheckM2, SingleM and
Metabuli. Which resources are required depends on the workflow stages and
options selected. Place shared databases outside individual run directories;
on an HPC system, coordinate their location with the system administrator.

Configure existing resources explicitly:

```bash
aviary configure \
  --gtdb-path /shared/db/gtdb/release232 \
  --eggnog-db-path /shared/db/eggnog \
  --checkm2-db-path /shared/db/checkm2/uniref100.KO.1.dmnd \
  --singlem-metapackage-path /shared/db/singlem/package.smpkg.zb \
  --metabuli-db-path /shared/db/metabuli
```

The corresponding environment variables are:

| Variable | Consumer |
| --- | --- |
| `GTDBTK_DATA_PATH` | GTDB-Tk taxonomy |
| `EGGNOG_DATA_DIR` | EggNOG functional annotation |
| `CHECKM2DB` | CheckM2 genome quality assessment |
| `SINGLEM_METAPACKAGE_PATH` | SingleM community/recovery analysis |
| `METABULI_DB_PATH` | Metabuli classification used by applicable workflows |
| `TMPDIR` | Temporary-file location |

To ask Aviary to download configured resources, add `--download` followed by
one or more of `gtdb`, `eggnog`, `singlem`, `checkm2` and `metabuli`.

```bash
aviary configure \
  --gtdb-path /shared/db/gtdb \
  --checkm2-db-path /shared/db/checkm2 \
  --download gtdb checkm2
```

!!! warning "Large downloads"
    Database downloads are large and should not be duplicated per user or run.
    Confirm available storage and release requirements before starting them.
    **Supplying `--download` with no values requests all five databases at
    once** — always list the specific databases you need.

Continue with the [quickstart](getting-started/quickstart.md). For production
clusters, see [HPC and cluster submission](guides/hpc.md).
