__version__ = "0.12.0"


# CONSTANTS
LONG_READ_TYPES = ["ont", "ont_hq", "rs", "sq", "ccs", "hifi"]
MEDAKA_MODELS = [
    "r103_fast_g507", "r103_fast_snp_g507", "r103_fast_variant_g507", "r103_hac_g507", "r103_hac_snp_g507",
    "r103_hac_variant_g507", "r103_min_high_g345", "r103_min_high_g360", "r103_prom_high_g360", "r103_prom_snp_g3210",
    "r103_prom_variant_g3210", "r103_sup_g507", "r103_sup_snp_g507", "r103_sup_variant_g507", "r1041_e82_260bps_fast_g632",
    "r1041_e82_260bps_fast_variant_g632", "r1041_e82_260bps_hac_g632", "r1041_e82_260bps_hac_variant_g632", "r1041_e82_260bps_sup_g632",
    "r1041_e82_260bps_sup_variant_g632", "r1041_e82_400bps_fast_g615", "r1041_e82_400bps_fast_g632",
    "r1041_e82_400bps_fast_variant_g615", "r1041_e82_400bps_fast_variant_g632", "r1041_e82_400bps_hac_g615",
    "r1041_e82_400bps_hac_g632", "r1041_e82_400bps_hac_variant_g615", "r1041_e82_400bps_hac_variant_g632", "r1041_e82_400bps_sup_g615",
    "r1041_e82_400bps_sup_variant_g615", "r104_e81_fast_g5015", "r104_e81_fast_variant_g5015", "r104_e81_hac_g5015",
    "r104_e81_hac_variant_g5015", "r104_e81_sup_g5015", "r104_e81_sup_g610", "r104_e81_sup_variant_g610", "r10_min_high_g303",
    "r10_min_high_g340", "r941_e81_fast_g514", "r941_e81_fast_variant_g514", "r941_e81_hac_g514", "r941_e81_hac_variant_g514",
    "r941_e81_sup_g514", "r941_e81_sup_variant_g514", "r941_min_fast_g303", "r941_min_fast_g507", "r941_min_fast_snp_g507",
    "r941_min_fast_variant_g507", "r941_min_hac_g507", "r941_min_hac_snp_g507", "r941_min_hac_variant_g507", "r941_min_high_g303",
    "r941_min_high_g330", "r941_min_high_g340_rle", "r941_min_high_g344", "r941_min_high_g351", "r941_min_high_g360", "r941_min_sup_g507",
    "r941_min_sup_snp_g507", "r941_min_sup_variant_g507", "r941_prom_fast_g303", "r941_prom_fast_g507", "r941_prom_fast_snp_g507",
    "r941_prom_fast_variant_g507", "r941_prom_hac_g507", "r941_prom_hac_snp_g507", "r941_prom_hac_variant_g507", "r941_prom_high_g303",
    "r941_prom_high_g330", "r941_prom_high_g344", "r941_prom_high_g360", "r941_prom_high_g4011", "r941_prom_snp_g303", "r941_prom_snp_g322",
    "r941_prom_snp_g360", "r941_prom_sup_g507", "r941_prom_sup_snp_g507", "r941_prom_sup_variant_g507", "r941_prom_variant_g303",
    "r941_prom_variant_g322", "r941_prom_variant_g360", "r941_sup_plant_g610", "r941_sup_plant_variant_g610"
]
