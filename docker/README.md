# Aviary Docker image

The image is intentionally software-only. Aviary databases are large and should
normally be mounted at runtime from a shared location.

## Build

From this directory:

```bash
pixi run bash ./build.sh
```

Or manually:

```bash
AVIARY_VERSION=0.13.0
sed "s/AVIARY_VERSION/$AVIARY_VERSION/g" Dockerfile.in > Dockerfile
DOCKER_BUILDKIT=1 docker build -t "wwood/aviary:$AVIARY_VERSION" .
```

## Runtime Databases

The image sets default database paths under `/db`:

```bash
GTDBTK_DATA_PATH=/db/gtdb
EGGNOG_DATA_DIR=/db/eggnog
SINGLEM_METAPACKAGE_PATH=/db/singlem/singlem.smpkg
CHECKM2DB=/db/checkm2/uniref100.KO.1.dmnd
METABULI_DB_PATH=/db/metabuli
```

Override them with `-e` or mount databases into those paths:

```bash
docker run --rm \
  -v /shared/db:/db:ro \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  wwood/aviary:$AVIARY_VERSION --help
```

For a full workflow, pass normal Aviary arguments after the image name.
