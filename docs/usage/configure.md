---
title: aviary configure
---

# aviary configure

Set conda environment variables for database paths used by future aviary runs. Settings are persisted across sessions.

```bash
aviary configure --gtdb-path ~/gtdbtk/release232/ --tmpdir /path/to/tmp/
```

> `configure` inherits the shared parser options, including `-t`, `-p`, `-n`,
> `-m`, `-o`, execution controls and the hidden `-w`/`--workflow`. In normal
> configuration mode they are accepted but ignored; `--tmpdir`, database paths
> and `--download` are the effective options. See
> [Shared options](centralised_commands.md).

## Examples

Configure GTDB and temp directory:
```bash
aviary configure --gtdb-path ~/gtdbtk/release232/ --tmpdir /scratch/tmp/
```

Download all databases:
```bash
aviary configure --download
```

Download specific databases:
```bash
aviary configure --download gtdb checkm2
```

View current configuration (run with no arguments):
```bash
aviary configure
```

## Database path options

**`--gtdb-path`** PATH

  Path to the local GTDB database files (used by GTDB-tk).

**`--checkm2-db-path`** PATH

  Path to the CheckM2 database.

**`--eggnog-db-path`** PATH

  Path to the local EggNOG database files.

**`--singlem-metapackage-path`** PATH

  Path to the local SingleM metapackage.

**`--metabuli-db-path`** PATH

  Path to the local Metabuli database.

**`--busco-db-path`** PATH

  Path to the local BUSCO database files.

**`--tmpdir`** PATH

  Path used for temporary files. Overrides the `TMPDIR` environment variable.

## Downloading databases

Use `--download` to fetch databases automatically:

```bash
aviary configure --download gtdb eggnog singlem checkm2 metabuli
```

Available databases: `gtdb`, `eggnog`, `singlem`, `checkm2`, `metabuli`. If no arguments are given, all databases are downloaded.

### Database versions

`--download` doesn't take a version — it fetches whatever this release of aviary is pinned to,
or (for eggNOG/SingleM) whatever the installed tool's own default is at download time:

| Database | Version this release targets | Where it's pinned |
|---|---|---|
| GTDB-tk | Release 232 | Hardcoded URL in the `download_gtdb` rule. Must stay in step with the `gtdbtk` conda package version in `pixi.toml` (its parser is release-specific) — if you bump one, bump the other. |
| CheckM2 | Whatever `checkm2 database --download` currently serves | Not pinned by aviary; follows the installed `checkm2` package (`>=1.1.0`). |
| EggNOG | Whatever `eggnog-mapper`'s installed version reports as its DB version | Read dynamically from `eggnogmapper.version.__DB_VERSION__` at download time — pin `eggnog-mapper` in `pixi.toml` if you need a specific DB version reproducibly. |
| SingleM metapackage | Whatever `singlem data` fetches for the installed `singlem` version | Not pinned by aviary; follows the installed `singlem` package. |
| Metabuli | An archived GTDB r214.1 + human T2T index | Fetched directly from an S3 archive path (upstream's `metabuli databases GTDB` command points at a relocated/404ing tarball, so aviary bypasses it). |

If you already have databases from a previous aviary version, check this table before assuming
they still match — GTDB-tk in particular will misclassify (or refuse to run) against a release
it wasn't built for.

## Gotchas

**`--checkm2-db-path` must point to the `.dmnd` file itself, not a directory.** This is the one
database path that doesn't follow the "point at a directory" pattern the other flags use (and
that aviary's own generic missing-database prompt tells you to follow) — CheckM2 reads
`CHECKM2DB` as a direct path to its diamond database file. `aviary configure --download checkm2`
handles this for you: `checkm2 database --download` normally creates
`<path>/CheckM2_database/*.dmnd`, and the download rule moves the `.dmnd` up to `<path>/`
afterwards. If you're pointing at a database you downloaded or symlinked yourself, make sure
`--checkm2-db-path` (or the `CHECKM2DB` environment variable) resolves all the way to the
`.dmnd` file, not the `CheckM2_database/` folder it usually ships in. By convention this file or
symlink is named `CheckM2_database` in aviary's own database layout.

**`--download gtdb` refuses to run into a non-empty directory.** It checks the target has at
most one entry before extracting, and exits with an error rather than merging into or
overwriting an existing GTDB installation. Point `--gtdb-path` at a fresh empty directory, or
clear out the old one first if you're intentionally replacing it.

**A mismatched GTDB-tk version against your GTDB release silently gives wrong or missing
taxonomy**, rather than a clear error, since older `gtdbtk` releases can't parse newer reference
package formats. If you supply your own `--gtdb-path` instead of using `--download`, match the
release to the `gtdbtk` version this aviary release pins in `pixi.toml`.

**Metabuli's database must land at `<path>/gtdb`.** The classify step (in `binning.smk`)
expects that exact subdirectory name under `--metabuli-db-path`; the archive aviary downloads
already extracts to that layout, but a manually-built database needs to match it too.

**Re-running `--download` for a database that already succeeded is a no-op**, since each
download rule is marked done via a `.download.done` sentinel file next to the database
directory. Delete that sentinel (or the whole database directory) to force a re-download.
