# Written by setuptools_scm from the git tag (see [tool.setuptools_scm] in
# pyproject.toml). It IS tracked in git, unlike the usual advice in its own
# header comment, because bird-release-python's flow depends on that: it
# force-writes this file, tolerates exactly " M aviary/version.py" in its
# pre-tag clean-tree guard, and then runs `git commit -a` -- which stages only
# tracked modifications, so an ignored version.py would leave nothing to commit
# and abort the release. singlem and bird_tool_utils track theirs for the same
# reason.
#
# The fallback still matters: the file can be absent in a fresh checkout before
# any build, and aviary.py imports __version__ from here for `--version` and its
# startup log line, so a bare ImportError would break the CLI outright.
try:
    from aviary.version import __version__
except ImportError:  # pragma: no cover - only hit in an unbuilt checkout
    __version__ = "0.0.0.dev0"


# CONSTANTS
LONG_READ_TYPES = ["ont", "ont_hq", "rs", "sq", "ccs", "hifi"]
LONG_READ_ASSEMBLERS = ["myloasm", "flye"]
# Aligner family used for long reads. rammap is a minimap2-compatible Rust
# implementation and accepts the same -x presets, so switching families does
# not change which preset is chosen for a given --long-read-type. minibwa has
# no CoverM preset suffix for long reads (unlike rammap/minimap2); its long-read
# behaviour is instead controlled via the raw --minibwa-params passthrough
# (e.g. "-x lr"), which the caller is responsible for setting appropriately.
LONG_READ_MAPPERS = ["rammap", "minimap2", "minibwa"]
# Short-read aligner family. CoverM defaults to strobealign from v0.7.0;
# minimap2 selects the '-x sr' preset instead. bwa-mem/bwa-mem2 are single
# fixed CoverM presets, same as strobealign/minibwa -- see MAPPERS_WITHOUT_MODELS.
#
# strobealign-aemb is not a real sibling of the others: it is not a `-p`
# value at all (CoverM's clap parser rejects it there -- verified against
# `coverm contig -p strobealign-aemb`, which lists 16 valid -p values, none
# of them this). It is a `-m`/--methods value that internally forces
# MappingProgram::STROBEALIGN regardless of -p (see
# strobealign_aemb_coverage() in CoverM's src/bin/coverm.rs) and shells out to
# `strobealign --aemb` directly, bypassing CoverM's normal BAM/pileup pipeline
# entirely. It is also `coverm contig`-only: `coverm genome` never special-
# cases it and panics on it (mosdepth_genome_coverage_estimators.rs). It is
# listed here anyway, as a --short-read-mapper CLI value, so get_coverage.py
# (the one coverm-contig call site that does per-contig binning coverage) can
# special-case it; every coverm-genome call site must instead use
# short_read_mapper_for_alignment() below, never this value directly.
SHORT_READ_MAPPERS = ["strobealign", "minimap2", "rammap", "minibwa", "bwa-mem", "bwa-mem2", "strobealign-aemb"]
# User-facing family -> the value CoverM's -p expects for short reads. Long
# reads resolve per --long-read-type instead, so they are not listed here.
SHORT_READ_MAPPER_TO_COVERM = {
    "strobealign": "strobealign",
    "minimap2": "minimap2-sr",
    "rammap": "rammap-sr",
    "minibwa": "minibwa",
    "bwa-mem": "bwa-mem",
    "bwa-mem2": "bwa-mem2",
    # Not a real -p value -- see the SHORT_READ_MAPPERS comment above. Kept
    # here as an identity entry so config["short_read_mapper"] stays the
    # literal string get_coverage.py switches on; short_read_mapper_for_alignment()
    # is what actually goes on any coverm-genome command line.
    "strobealign-aemb": "strobealign-aemb",
}

# Anywhere a real per-read alignment is needed -- coverm genome (per-genome
# abundance), a raw `coverm contig -m metabat -p ...` coverage call, or
# polish.py building a racon PAF -- cannot use strobealign-aemb (see the
# SHORT_READ_MAPPERS comment above: it is coverm-contig -m-only and not a -p
# value at all). Fall back to the real aligner it runs under the hood instead
# of the pseudo-mapper name.
GENOME_COVERM_MAPPER_OVERRIDE = {"strobealign-aemb": "strobealign"}


def short_read_mapper_for_alignment(short_read_mapper_coverm: str) -> str:
    """Translate a resolved config["short_read_mapper"] value (already a
    CoverM -p-shaped string, e.g. from SHORT_READ_MAPPER_TO_COVERM) into one
    that is safe to use anywhere a real alignment (BAM or PAF) is produced --
    coverm genome, a raw coverm-contig call outside get_coverage.py, or
    polish.py. Identity for every real -p value; only strobealign-aemb is
    remapped, since none of those can run it.
    """
    return GENOME_COVERM_MAPPER_OVERRIDE.get(short_read_mapper_coverm, short_read_mapper_coverm)


# Mapper families that expose a real model/preset axis through CoverM's -p
# (e.g. rammap-ont, rammap-hifi). Verified against `coverm contig --full-help`:
# only minimap2 and rammap have more than one preset -- strobealign, minibwa,
# bwa-mem/bwa-mem2 and strobealign-aemb are each a single fixed preset (or, for
# strobealign-aemb, not a preset at all) with no model to choose.
MAPPER_MODELS = {
    "minimap2": ["sr", "lr-hq", "ont", "pb", "hifi", "no-preset"],
    "rammap": ["sr", "lr-hq", "ont", "pb", "hifi", "no-preset"],
}
MAPPERS_WITHOUT_MODELS = ("strobealign", "bwa-mem", "bwa-mem2", "minibwa", "strobealign-aemb")

# Which models make sense for which read length. CoverM will happily accept
# `-p minimap2-ont` for short reads: it returns a well-formed coverage table of
# near-zero depths rather than failing, so the mistake surfaces as wrong numbers
# rather than an error. Splitting the sets lets processor.py reject the
# combination up front instead. "no-preset" is in both because it means "pass no
# -x at all", which is a deliberate choice at either read length.
SHORT_READ_MODELS = ("sr", "no-preset")
LONG_READ_MODELS = ("lr-hq", "ont", "pb", "hifi", "no-preset")

# Defaults for --short-read-mapper / --long-read-mapper. These live here rather
# than in argparse's `default=` so that argparse can default to None, which is
# what lets processor.py tell "user asked for this mapper" apart from "nobody
# said anything" -- required to error on a mapper selected without the reads it
# would map. Resolved to these values immediately afterwards.
DEFAULT_SHORT_READ_MAPPER = "strobealign"
DEFAULT_LONG_READ_MAPPER = "rammap"

# CoverM's -p suffix is NOT the aligner's own -x preset name -- CoverM
# translates internally. Anything that builds a raw minimap2/rammap command
# line (polish.py, for racon's PAF) has to do that translation itself. Verified
# against `coverm contig --full-help` and both aligners' own error messages:
# only "sr" happens to be spelled the same on both sides, so a naive
# removeprefix() of the CoverM name produces `-x lr-hq` / `-x ont` / ... which
# minimap2 and rammap both reject with "unknown preset".
# None means "emit no -x flag at all", which is what CoverM's "no-preset" does.
COVERM_MODEL_TO_ALIGNER_PRESET = {
    "sr": "sr",
    "lr-hq": "lr:hq",
    "ont": "map-ont",
    "pb": "map-pb",
    "hifi": "map-hifi",
    "no-preset": None,
}

# Default model per --long-read-type when --long-read-mapper-model is not
# given, used by resolve_mapper_model(). Previously each of get_coverage.py,
# get_abundances.py and fraction_recovered.py independently lumped ont_hq in
# with ont and hifi in with pb (rs/sq/ccs/hifi -> "-pb"); CoverM's lr-hq and
# hifi presets existed but were unreachable. This table is the single source
# of truth now, and fixes that gap as a side effect.
#
# "ont" -> "lr-hq" (not the legacy "ont"/map-ont model) and "ccs" -> "hifi"
# (not the legacy "pb"/map-pb model): minimap2's own docs describe map-ont as
# assuming a ~10% error rate, and map-pb as "effectively deprecated by HiFi"
# for anything but very old data -- both legacy presets built for accuracy
# levels well below current chemistry. lr:hq is what ONT's own developers
# recommend for chemistry v14 reads (~99% accuracy) and minimap2's docs say it
# outperforms map-hifi there; CCS reads are, by PacBio's current terminology,
# HiFi reads, so routing them through the deprecated CLR preset was the same
# class of bug. rs/sq (RSII/Sequel) stay on "pb" -- that is genuine older CLR
# chemistry, not deprecated-by-mistake. The legacy "ont"/"pb" models remain
# fully selectable via --long-read-mapper-model for anyone who wants them.
LONG_READ_TYPE_TO_MODEL = {
    "ont": "lr-hq",
    "ont_hq": "lr-hq",
    "rs": "pb",
    "sq": "pb",
    "ccs": "hifi",
    "hifi": "hifi",
}

# spades_assembly.py's --long-read-type only accepts spades.py's own
# --nanopore/--pacbio vocabulary ("ont", "ont_hq", "pacbio", "pacbio_hifi"),
# narrower than aviary's --longread-type (LONG_READ_TYPES above). Passing one
# of the other four values straight through crashes it with "invalid choice"
# -- verified by actually running `aviary assemble --longread-type ccs`.
# rs/sq are older PacBio CLR chemistry -> "pacbio"; ccs/hifi are HiFi ->
# "pacbio_hifi", matching the rs/sq -> "pb", ccs/hifi -> "hifi" split above.
LONG_READ_TYPE_TO_SPADES = {
    "ont": "ont",
    "ont_hq": "ont_hq",
    "rs": "pacbio",
    "sq": "pacbio",
    "ccs": "pacbio_hifi",
    "hifi": "pacbio_hifi",
}


def resolve_mapper_model(mapper: str, model: str | None, long_read_type: str | None = None) -> str:
    """Combine a mapper family with an explicit model or a --long-read-type default
    into the CoverM -p value. `model` (from --short-read-mapper-model /
    --long-read-mapper-model) takes precedence when given; otherwise
    `long_read_type` is used to look up the default via LONG_READ_TYPE_TO_MODEL.
    Callers are expected to have already validated mapper/model compatibility
    (see processor.py) -- this raises only as a last-resort guard.
    """
    if mapper not in MAPPER_MODELS:
        raise ValueError(f"Mapper {mapper!r} has no selectable model.")
    if model is None:
        if long_read_type is None:
            raise ValueError(f"No model given and no long_read_type to derive a default for mapper {mapper!r}.")
        model = LONG_READ_TYPE_TO_MODEL[long_read_type]
    if model not in MAPPER_MODELS[mapper]:
        raise ValueError(f"Invalid model {model!r} for mapper {mapper!r}. Valid models: {MAPPER_MODELS[mapper]}")
    return f"{mapper}-{model}"


def short_read_mapper_extra_params(
    short_read_mapper_coverm: str,
    bwa_params: str | None = None,
    strobealign_params: str | None = None,
    minimap2_params: str | None = None,
    rammap_params: str | None = None,
) -> list[str]:
    """Pick CoverM's raw per-aligner passthrough flag (--bwa-params/
    --strobealign-params/--minimap2-params/--rammap-params) for an
    already-resolved config["short_read_mapper"]-shaped value (e.g.
    "bwa-mem2", "minimap2-sr", "strobealign"), as extra argv tokens. Shared by
    get_coverage.py, get_abundances.py and fraction_recovered.py so the
    family-matching logic isn't duplicated three times.

    strobealign-aemb is deliberately excluded: it shells out to
    `strobealign --aemb` directly inside CoverM (see SHORT_READ_MAPPERS'
    comment above) and does not accept --strobealign-params.
    """
    if short_read_mapper_coverm in ("bwa-mem", "bwa-mem2") and bwa_params:
        return ["--bwa-params", bwa_params]
    if short_read_mapper_coverm == "strobealign" and strobealign_params:
        return ["--strobealign-params", strobealign_params]
    if short_read_mapper_coverm.startswith("minimap2") and minimap2_params:
        return ["--minimap2-params", minimap2_params]
    if short_read_mapper_coverm.startswith("rammap") and rammap_params:
        return ["--rammap-params", rammap_params]
    return []


def aligner_preset_args(model: str) -> list[str]:
    """The `-x` arguments for a raw minimap2/rammap command line, given the
    model half of a CoverM -p name (e.g. "lr-hq" -> ["-x", "lr:hq"]). Returns
    an empty list for "no-preset", which means the aligner runs with no -x.
    See COVERM_MODEL_TO_ALIGNER_PRESET for why this translation is needed.
    """
    if model not in COVERM_MODEL_TO_ALIGNER_PRESET:
        raise ValueError(
            f"Unknown mapper model {model!r}. Valid models: "
            f"{sorted(COVERM_MODEL_TO_ALIGNER_PRESET)}"
        )
    preset = COVERM_MODEL_TO_ALIGNER_PRESET[model]
    return [] if preset is None else ["-x", preset]


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
COVERAGE_JOB_STRATEGIES = ["default", "always", "never"]
COVERAGE_JOB_CUTOFF = 10
