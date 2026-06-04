#!/bin/bash -eo pipefail
# To execute this script, ensure the version tag has been pushed to GitHub. Then:
#   pixi run bash ./build.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
export AVIARY_VERSION=`aviary --version | awk '{print $NF}'`
export AVIARY_DOCKER_VERSION=ghcr.io/snh-star/aviary:$AVIARY_VERSION
sed "s/AVIARY_VERSION/$AVIARY_VERSION/g" Dockerfile.in > Dockerfile && \
DOCKER_BUILDKIT=1 docker build -t $AVIARY_DOCKER_VERSION . && \
docker run --rm -v `pwd`:`pwd` $AVIARY_DOCKER_VERSION assemble \
    -1 `pwd`/test/data/wgsim.1.fq.gz \
    -2 `pwd`/test/data/wgsim.2.fq.gz \
    --use-megahit \
    --skip-qc \
    -o `pwd`/aviary_docker_test_out \
    -n 4 -t 4
docker tag $AVIARY_DOCKER_VERSION ghcr.io/snh-star/aviary:v$AVIARY_VERSION
docker tag $AVIARY_DOCKER_VERSION ghcr.io/snh-star/aviary:latest
echo "Seems good - now you just need to:"
echo "  docker push $AVIARY_DOCKER_VERSION"
echo "  docker push ghcr.io/snh-star/aviary:v$AVIARY_VERSION"
echo "  docker push ghcr.io/snh-star/aviary:latest"
