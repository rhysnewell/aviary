#!/usr/bin/env bash

SINGLEM_METAPACKAGE_PATH_RELATIVE="db/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb"
CHECKM_DATA_PATH_RELATIVE="db/2015_01_16_v2"
GTDBTK_DATA_PATH_RELATIVE="db/release220"
CHECKM2DB_RELATIVE="db/CheckM2_database"
EGGNOG_DATA_DIR_RELATIVE="db/2.1.3"
METABULI_DB_PATH_RELATIVE="db/2024-3-28-GTDB214.1+humanT2T"

# For each of the above variables, check if the file exists in the relative path
# and if it does, set the corresponding environment variable to the absolute path.
for var in SINGLEM_METAPACKAGE_PATH \
           CHECKM_DATA_PATH \
           GTDBTK_DATA_PATH \
           CHECKM2DB \
           EGGNOG_DATA_DIR \
           METABULI_DB_PATH; do
    relative_path="${var}_RELATIVE"
    if [[ -e "${!relative_path}" ]]; then
        export "$var"="$(realpath "${!relative_path}")"
    else
        echo "File $(basename "${!relative_path}") not found in db/" >&2
    fi
done
