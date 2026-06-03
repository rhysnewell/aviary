# Aviary Docker image

The image includes **SingleM** and **CheckM2** databases (~5GB combined), enabling
binning-only and assembly workflows out of the box. GTDB-Tk, EggNOG, and Metabuli
are not bundled due to their size (~140GB combined) but can be mounted at runtime.

## Included databases

| Database | Path in image | Required for |
|---|---|---|
| SingleM metapackage | `/db/singlem` | Read fraction estimation |
| CheckM2 | `/db/checkm2` | Bin quality assessment |

## Commands that work without extra databases

```bash
aviary assemble
aviary recover --binning-only
```

## Commands requiring mounted databases

| Command | Requires |
|---|---|
| `aviary recover` (full) | GTDB-Tk (`--gtdb-path`) |
| `aviary annotate` | EggNOG (`--eggnog-db-path`) |
| Taxonomic classification | Metabuli (`--metabuli-db-path`) |

## Build

From this directory:

```bash
pixi run bash ./build.sh
```

## Runtime — no extra databases

```bash
docker run --rm \
  -v "$PWD":"$PWD" -w "$PWD" \
  ghcr.io/snh-star/aviary:0.13.0 \
  recover --binning-only -1 reads.1.fq.gz -2 reads.2.fq.gz -o output/
```

## Runtime — mounting additional databases

```bash
docker run --rm \
  -v /shared/db/gtdb:/db/gtdb:ro \
  -v /shared/db/eggnog:/db/eggnog:ro \
  -v /shared/db/metabuli:/db/metabuli:ro \
  -v "$PWD":"$PWD" -w "$PWD" \
  ghcr.io/snh-star/aviary:0.13.0 \
  recover -1 reads.1.fq.gz -2 reads.2.fq.gz -o output/ \
  --gtdb-path /db/gtdb \
  --eggnog-db-path /db/eggnog \
  --metabuli-db-path /db/metabuli
```
