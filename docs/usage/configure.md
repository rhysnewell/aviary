---
title: aviary configure
---

# aviary configure

Set conda environment variables for database paths used by future aviary runs. Settings are persisted across sessions.

```
aviary configure --gtdb-path ~/gtdbtk/release207/ --tmpdir /path/to/tmp/
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

```
aviary configure --download gtdb eggnog singlem checkm2 metabuli
```

Available databases: `gtdb`, `eggnog`, `singlem`, `checkm2`, `metabuli`. If no arguments are given, all databases are downloaded.

## Examples

Configure GTDB and temp directory:
```
aviary configure --gtdb-path ~/gtdbtk/release207/ --tmpdir /scratch/tmp/
```

Download all databases:
```
aviary configure --download
```

Download specific databases:
```
aviary configure --download gtdb checkm2
```

View current configuration (run with no arguments):
```
aviary configure
```
