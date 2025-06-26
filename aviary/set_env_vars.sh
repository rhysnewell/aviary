#!/usr/bin/env bash

SINGLEM_METAPACKAGE_PATH_RELATIVE="$PIXI_PROJECT_ROOT/../db/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb"
CHECKM_DATA_PATH_RELATIVE="$PIXI_PROJECT_ROOT/../db/2015_01_16_v2"
GTDBTK_DATA_PATH_RELATIVE="$PIXI_PROJECT_ROOT/../db/release220"
CHECKM2DB_RELATIVE="$PIXI_PROJECT_ROOT/../db/CheckM2_database"
EGGNOG_DATA_DIR_RELATIVE="$PIXI_PROJECT_ROOT/../db/2.1.3"
METABULI_DB_PATH_RELATIVE="$PIXI_PROJECT_ROOT/../db/2024-3-28-GTDB214.1+humanT2T"

# For each of the above variables, check if the file exists in the relative path
# and if it does, set the corresponding environment variable to the absolute path.

rm set_env_vars.log
for var in SINGLEM_METAPACKAGE_PATH \
           CHECKM_DATA_PATH \
           GTDBTK_DATA_PATH \
           CHECKM2DB \
           EGGNOG_DATA_DIR \
           METABULI_DB_PATH; do
    relative_path="${var}_RELATIVE"
    if [[ -e "${!relative_path}" ]]; then
        export "$var"="$(realpath "${!relative_path}")"
    elif [[ -n "${!var}" ]]; then
        # If the variable is already set, use it, all good
        echo "Using preset for $var: ${!var}"
    else
        # We could croak here so that any errors are caught, but the error from
        # pixi is very vague, and unsure how to report the actual error here, so
        # would cause too much confusion.
        echo PWD: `pwd` >>set_env_vars.log
        echo "File ${!relative_path} not found" >>set_env_vars.log
    fi
done
