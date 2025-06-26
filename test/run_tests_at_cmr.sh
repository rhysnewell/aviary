#!/bin/bash

set -euo pipefail

# Run from base dir of repo
# $ bash test/run_tests_at_cmr.sh
#
# For now, need to make sure manually that all the mqsub'd jobs finish without error at the end.

MQSUB_OUT="mqsub.out"

# GPU testing - taxvamb
echo mqsub GPU testing - taxvamb ..
mqsub -t32 --A100 --segregated-log-files --bg -- pixi run --frozen --manifest-path aviary/pixi.toml -e dev "'export TEST_REQUEST_GPU=1 && python test/test_integration.py Tests.test_short_read_recovery_taxvamb'" &>> $MQSUB_OUT

# GPU testing - comebin
echo mqsub GPU testing - comebin ..
mqsub -t32 --A100 --segregated-log-files --bg -- pixi run --frozen --manifest-path aviary/pixi.toml -e dev "'export TEST_REQUEST_GPU=1 && python test/test_integration.py Tests.test_short_read_recovery_comebin'" &>> $MQSUB_OUT

# GPU testing - semibin
echo mqsub GPU testing - semibin ..
mqsub -t32 --A100 --segregated-log-files --bg -- pixi run --frozen --manifest-path aviary/pixi.toml -e dev "'export TEST_REQUEST_GPU=1 && python test/test_integration.py Tests.test_short_read_recovery_semibin'" &>> $MQSUB_OUT

# Expensive tests
echo mqsub Expensive tests ..
mqsub -t 16 -m 128 --segregated-log-files --bg -- pixi run --frozen --manifest-path aviary/pixi.toml -e dev pytest --run-expensive test/test_integration.py &>> $MQSUB_OUT


echo Running qsub tests locally, which qsub themselves ..
# TODO: The non-expensive tests are actually run twice, with --run-qsub, and
# with --run-expensive above But this isn't a lot of time extra time.
pixi run --frozen --manifest-path aviary/pixi.toml -e dev pytest test --run-qsub 2>&1 | tee test/qsub-test-log.txt
echo "Qsub tests appear to have run successfully, but check the log file test/qsub-test-log.txt for details."

echo
echo "Need to mqwait for these jobs to finish and check they succeeded. TODO: Make this automatic."
cat $MQSUB_OUT
