#!/bin/bash -eo pipefail
# The -eo pipefail above only takes effect when this script is run directly
# as ./build.sh. If invoked as `bash ./build.sh` (e.g. via `pixi run bash
# ./build.sh`), the shebang's flags are ignored by bash, so we set them
# explicitly here too. We ALSO check exit codes explicitly below rather than
# relying solely on set -e propagating through chained/multi-line commands,
# since that chaining has proven unreliable here (e.g. under pixi run
# wrapping, or if the build is cancelled mid-way) - this guarantees the
# "Seems good" success message can never print unless every step actually
# succeeded.
set -eo pipefail

# To execute this script, ensure the version tag has been pushed to GitHub. Then:
#   pixi run bash ./build.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

export AVIARY_VERSION=`aviary --version | awk '{print $NF}'`
export AVIARY_DOCKER_VERSION=ghcr.io/snh-star/aviary:$AVIARY_VERSION

fail() {
    echo "" >&2
    echo "BUILD FAILED: $1" >&2
    echo "No tags were applied and nothing needs to be pushed." >&2
    exit 1
}

# Download test reads if not already present
if [ ! -f wgsim.1.fq.gz ]; then
    curl -fsSL -o wgsim.1.fq.gz https://github.com/rhysnewell/aviary/raw/main/test/data/wgsim.1.fq.gz \
        || fail "could not download wgsim.1.fq.gz"
fi
if [ ! -f wgsim.2.fq.gz ]; then
    curl -fsSL -o wgsim.2.fq.gz https://github.com/rhysnewell/aviary/raw/main/test/data/wgsim.2.fq.gz \
        || fail "could not download wgsim.2.fq.gz"
fi

sed "s/AVIARY_VERSION/$AVIARY_VERSION/g" Dockerfile.in > Dockerfile \
    || fail "sed substitution into Dockerfile failed"

DOCKER_BUILDKIT=1 docker build -t $AVIARY_DOCKER_VERSION . \
    || fail "docker build did not complete successfully"

docker run --rm \
    -e GTDBTK_DATA_PATH=/tmp \
    -e EGGNOG_DATA_DIR=/tmp \
    -e METABULI_DB_PATH=/tmp \
    -v `pwd`:`pwd` $AVIARY_DOCKER_VERSION assemble \
    -1 `pwd`/wgsim.1.fq.gz \
    -2 `pwd`/wgsim.2.fq.gz \
    --use-megahit \
    --skip-qc \
    -o `pwd`/aviary_docker_test_out \
    -n 4 -t 4 \
    || fail "test run (aviary assemble inside the built image) did not complete successfully"

# Only reachable if build AND test run both genuinely succeeded.
docker tag $AVIARY_DOCKER_VERSION ghcr.io/snh-star/aviary:v$AVIARY_VERSION \
    || fail "docker tag (v$AVIARY_VERSION) failed"
docker tag $AVIARY_DOCKER_VERSION ghcr.io/snh-star/aviary:latest \
    || fail "docker tag (latest) failed"

echo ""
echo "Build and test run both succeeded - now you just need to:"
echo "  docker push $AVIARY_DOCKER_VERSION"
echo "  docker push ghcr.io/snh-star/aviary:v$AVIARY_VERSION"
echo "  docker push ghcr.io/snh-star/aviary:latest"