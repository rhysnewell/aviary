#!/bin/bash -eo pipefail

# To execute this script, ensure the version tag has been pushed to GitHub. Then:
#   pixi run bash ./build.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

export AVIARY_VERSION=`aviary --version | awk '{print $NF}'`
export AVIARY_DOCKER_VERSION=wwood/aviary:$AVIARY_VERSION

# Download test reads if not already present
if [ ! -f wgsim.1.fq.gz ]; then
    wget --quiet https://github.com/rhysnewell/aviary/raw/main/test/data/wgsim.1.fq.gz
fi
if [ ! -f wgsim.2.fq.gz ]; then
    wget --quiet https://github.com/rhysnewell/aviary/raw/main/test/data/wgsim.2.fq.gz
fi

sed "s/AVIARY_VERSION/$AVIARY_VERSION/g" Dockerfile.in > Dockerfile && \
DOCKER_BUILDKIT=1 docker build -t $AVIARY_DOCKER_VERSION . && \
docker run --rm -v `pwd`:`pwd` $AVIARY_DOCKER_VERSION assemble \
    -1 `pwd`/wgsim.1.fq.gz \
    -2 `pwd`/wgsim.2.fq.gz \
    --use-megahit \
    --skip-qc \
    -o `pwd`/aviary_docker_test_out \
    -n 4 -t 4

echo "Seems good - now you just need to 'docker push $AVIARY_DOCKER_VERSION' to upload the image to Docker Hub"
