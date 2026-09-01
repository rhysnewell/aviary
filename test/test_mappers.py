#!/usr/bin/env python
"""Fast unit tests for the read-mapper selection added by --short-read-mapper
and --long-read-mapper.

These deliberately avoid running the pipeline: they cover the pure logic that
decides *which* aligner is invoked and *how* its command line is spelled. That
logic is easy to get wrong and hard to notice, because a wrong mapper name
surfaces only when a real coverage or polishing job runs -- and a short-read
aligner given long reads can return a well-formed table of near-zero depths
rather than failing.

The corresponding end-to-end coverage lives in test_integration.py.
"""

import importlib.util
import os
import tempfile
import unittest
from unittest import mock

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load(module_name, relative_path):
    """Import a script by path.

    The coverage/polishing scripts are executed by Snakemake as standalone
    files inside their own pixi environments, so they are not importable as
    part of the aviary package.
    """
    spec = importlib.util.spec_from_file_location(
        module_name, os.path.join(REPO_ROOT, relative_path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


aviary_init = _load("aviary_init", "aviary/__init__.py")
polish = _load("polish", "aviary/modules/assembly/scripts/polish.py")
get_abundances = _load(
    "get_abundances", "aviary/modules/binning/scripts/get_abundances.py")

# Short-read mappers with no PAF mode: polish.py routes bwa-mem/bwa-mem2
# through bwa_mem_paf_cmd() (index + mem + sam2paf pipeline) instead of
# short_read_paf_cmd(), which only knows single-command aligners.
# strobealign-aemb has no PAF mode at all -- it isn't a real alignment (see
# aviary/__init__.py's SHORT_READ_MAPPERS comment), so it can never be used
# for polishing; processor.py routes it to plain strobealign there via
# short_read_mapper_for_alignment() before it ever reaches polish.py.
MAPPERS_WITHOUT_PAF_COMMAND = ("bwa-mem", "bwa-mem2", "strobealign-aemb")


class TestMapperConstants(unittest.TestCase):
    def test_every_short_read_mapper_resolves_to_a_coverm_name(self):
        # A choice without a mapping raises KeyError in processor.py at
        # config-write time, i.e. after the user has already started a run.
        for mapper in aviary_init.SHORT_READ_MAPPERS:
            self.assertIn(mapper, aviary_init.SHORT_READ_MAPPER_TO_COVERM)

    def test_no_orphan_entries_in_the_mapping(self):
        for mapper in aviary_init.SHORT_READ_MAPPER_TO_COVERM:
            self.assertIn(mapper, aviary_init.SHORT_READ_MAPPERS)

    def test_defaults_are_the_new_mappers(self):
        self.assertEqual(aviary_init.SHORT_READ_MAPPERS[0], "strobealign")
        self.assertEqual(aviary_init.LONG_READ_MAPPERS[0], "rammap")

    def test_short_read_coverm_names(self):
        self.assertEqual(aviary_init.SHORT_READ_MAPPER_TO_COVERM, {
            "strobealign": "strobealign",
            "minimap2": "minimap2-sr",
            "rammap": "rammap-sr",
            "minibwa": "minibwa",
            "bwa-mem": "bwa-mem",
            "bwa-mem2": "bwa-mem2",
            # Not a real -p value -- kept as an identity entry so
            # config["short_read_mapper"] stays the literal string
            # get_coverage.py switches on. See short_read_mapper_for_alignment().
            "strobealign-aemb": "strobealign-aemb",
        })

    def test_strobealign_aemb_falls_back_to_strobealign_for_alignment(self):
        # coverm genome / polish.py / any raw coverm-contig call outside
        # get_coverage.py cannot run strobealign-aemb -- see
        # short_read_mapper_for_alignment()'s docstring.
        self.assertEqual(
            aviary_init.short_read_mapper_for_alignment("strobealign-aemb"), "strobealign")

    def test_short_read_mapper_for_alignment_is_identity_for_real_mappers(self):
        for coverm_name in aviary_init.SHORT_READ_MAPPER_TO_COVERM.values():
            if coverm_name == "strobealign-aemb":
                continue
            self.assertEqual(aviary_init.short_read_mapper_for_alignment(coverm_name), coverm_name)


class TestShortReadMapperExtraParams(unittest.TestCase):
    """short_read_mapper_extra_params() picks CoverM's raw per-aligner
    passthrough (--bwa-params/--strobealign-params/--minimap2-params/
    --rammap-params), by family, from an already-resolved
    config["short_read_mapper"]-shaped value."""

    ALL_PARAMS = dict(
        bwa_params="-x ont2d", strobealign_params="-U", minimap2_params="-k 21", rammap_params="-H")

    def test_picks_bwa_params_for_bwa_mem_and_bwa_mem2(self):
        for mapper in ("bwa-mem", "bwa-mem2"):
            self.assertEqual(
                aviary_init.short_read_mapper_extra_params(mapper, **self.ALL_PARAMS),
                ["--bwa-params", "-x ont2d"])

    def test_picks_strobealign_params_for_strobealign_only(self):
        self.assertEqual(
            aviary_init.short_read_mapper_extra_params("strobealign", **self.ALL_PARAMS),
            ["--strobealign-params", "-U"])

    def test_strobealign_aemb_never_gets_strobealign_params(self):
        # strobealign-aemb shells out to `strobealign --aemb` directly inside
        # CoverM and does not accept --strobealign-params.
        self.assertEqual(
            aviary_init.short_read_mapper_extra_params("strobealign-aemb", **self.ALL_PARAMS), [])

    def test_picks_minimap2_params_for_any_minimap2_model(self):
        for mapper in ("minimap2-sr", "minimap2-no-preset"):
            self.assertEqual(
                aviary_init.short_read_mapper_extra_params(mapper, **self.ALL_PARAMS),
                ["--minimap2-params", "-k 21"])

    def test_picks_rammap_params_for_any_rammap_model(self):
        for mapper in ("rammap-sr", "rammap-no-preset"):
            self.assertEqual(
                aviary_init.short_read_mapper_extra_params(mapper, **self.ALL_PARAMS),
                ["--rammap-params", "-H"])

    def test_minibwa_gets_nothing_here(self):
        # minibwa's passthrough is --minibwa-params, handled separately in
        # every caller (it also applies to --long-read-mapper minibwa, which
        # has no config["short_read_mapper"]-shaped value at all).
        self.assertEqual(
            aviary_init.short_read_mapper_extra_params("minibwa", **self.ALL_PARAMS), [])

    def test_no_params_given_returns_empty(self):
        self.assertEqual(aviary_init.short_read_mapper_extra_params("strobealign"), [])


class TestMapperModels(unittest.TestCase):
    """--short-read-mapper-model / --long-read-mapper-model only apply to
    mappers with more than one CoverM preset (minimap2, rammap)."""

    def test_mappers_without_models_have_no_entry(self):
        for mapper in aviary_init.MAPPERS_WITHOUT_MODELS:
            self.assertNotIn(mapper, aviary_init.MAPPER_MODELS)

    def test_every_mapper_is_classified_one_way_or_the_other(self):
        for mapper in aviary_init.SHORT_READ_MAPPERS:
            in_models = mapper in aviary_init.MAPPER_MODELS
            in_without = mapper in aviary_init.MAPPERS_WITHOUT_MODELS
            self.assertTrue(in_models or in_without, mapper)
            self.assertFalse(in_models and in_without, mapper)

    def test_resolve_uses_explicit_model_over_long_read_type(self):
        self.assertEqual(
            aviary_init.resolve_mapper_model("rammap", "no-preset", "ont"),
            "rammap-no-preset")

    def test_resolve_derives_default_from_long_read_type(self):
        # ont defaults to the modern lr-hq preset, not the legacy noisy "ont"
        # (map-ont, ~10% error) preset -- see LONG_READ_TYPE_TO_MODEL's
        # comment in aviary/__init__.py. rs/sq (genuine older PacBio CLR
        # chemistry, not deprecated-by-mistake) are unaffected and stay on pb.
        self.assertEqual(aviary_init.resolve_mapper_model("rammap", None, "ont"), "rammap-lr-hq")
        self.assertEqual(aviary_init.resolve_mapper_model("minimap2", None, "rs"), "minimap2-pb")

    def test_resolve_still_reaches_the_legacy_ont_and_pb_presets_explicitly(self):
        # Demoted from the default, not removed -- still selectable via
        # --long-read-mapper-model.
        self.assertEqual(aviary_init.resolve_mapper_model("rammap", "ont", "ont"), "rammap-ont")
        self.assertEqual(aviary_init.resolve_mapper_model("rammap", "pb", "ont"), "rammap-pb")

    def test_resolve_fixes_the_previously_dead_lr_hq_and_hifi_presets(self):
        # ont_hq used to be lumped in with ont, and hifi with pb (rs/sq/ccs);
        # CoverM's lr-hq/hifi presets existed but were unreachable via any
        # --long-read-type. This is the regression test for that fix.
        self.assertEqual(aviary_init.resolve_mapper_model("rammap", None, "ont_hq"), "rammap-lr-hq")
        self.assertEqual(aviary_init.resolve_mapper_model("rammap", None, "hifi"), "rammap-hifi")

    def test_resolve_ccs_defaults_to_hifi_not_legacy_pb(self):
        # CCS reads are HiFi reads by current PacBio terminology; map-pb is
        # the deprecated-by-HiFi CLR preset (see minimap2's own docs).
        self.assertEqual(aviary_init.resolve_mapper_model("minimap2", None, "ccs"), "minimap2-hifi")

    def test_resolve_rejects_a_model_for_a_modelless_mapper(self):
        with self.assertRaises(ValueError):
            aviary_init.resolve_mapper_model("strobealign", "sr", "ont")

    def test_resolve_rejects_an_invalid_model(self):
        with self.assertRaises(ValueError):
            aviary_init.resolve_mapper_model("rammap", "not-a-real-model", None)


class TestShortReadPafCommand(unittest.TestCase):
    """polish.py generates racon's PAF itself, not via CoverM.

    Each aligner spells PAF output differently, so these are the exact command
    shapes rather than a pattern: minimap2/rammap take a preset as the value of
    -x, strobealign's -x is a boolean switch, and minibwa needs the 'map'
    subcommand plus -f.
    """

    def test_strobealign_uses_boolean_x(self):
        self.assertEqual(
            polish.short_read_paf_cmd("strobealign", "ref.fa", "r.fq.gz", 8),
            "strobealign -x -t 8 ref.fa r.fq.gz")

    def test_minimap2_uses_sr_preset(self):
        self.assertEqual(
            polish.short_read_paf_cmd("minimap2-sr", "ref.fa", "r.fq.gz", 8),
            "minimap2 -x sr -t 8 ref.fa r.fq.gz")

    def test_rammap_uses_sr_preset(self):
        self.assertEqual(
            polish.short_read_paf_cmd("rammap-sr", "ref.fa", "r.fq.gz", 8),
            "rammap -x sr -t 8 ref.fa r.fq.gz")

    def test_minibwa_uses_map_subcommand_and_f(self):
        self.assertEqual(
            polish.short_read_paf_cmd("minibwa", "ref.fa", "r.fq.gz", 8),
            "minibwa map -f -x sr -t 8 ref.fa r.fq.gz")

    def test_minimap2_model_override_selects_the_given_preset(self):
        # short_read_paf_cmd keys on the fully-resolved CoverM -p value, which
        # includes a --short-read-mapper-model override when one was given.
        # "no-preset" means the aligner runs with no -x at all.
        self.assertEqual(
            polish.short_read_paf_cmd("minimap2-no-preset", "ref.fa", "r.fq.gz", 8),
            "minimap2 -t 8 ref.fa r.fq.gz")

    def test_coverm_model_names_are_translated_to_aligner_presets(self):
        # CoverM's -p suffix is not the aligner's own -x value: only "sr" is
        # spelled the same on both sides. Passing the CoverM spelling straight
        # through produced `-x lr-hq` / `-x ont` / ..., which real minimap2 and
        # rammap both reject with "unknown preset", failing racon mid-polish.
        for model, expected_preset in (
            ("sr", "sr"), ("lr-hq", "lr:hq"), ("ont", "map-ont"),
            ("pb", "map-pb"), ("hifi", "map-hifi"),
        ):
            for family in ("minimap2", "rammap"):
                self.assertEqual(
                    polish.short_read_paf_cmd(
                        f"{family}-{model}", "ref.fa", "r.fq.gz", 8),
                    f"{family} -x {expected_preset} -t 8 ref.fa r.fq.gz")

    def test_unknown_model_suffix_raises(self):
        with self.assertRaises(ValueError):
            polish.short_read_paf_cmd("minimap2-bogus", "ref.fa", "r.fq.gz", 8)

    def test_every_paf_capable_short_read_mapper_has_a_paf_command(self):
        # Any --short-read-mapper choice must be usable for racon polishing,
        # otherwise a valid flag combination fails partway through a run.
        # bwa-mem/bwa-mem2 are excluded: they have no PAF mode and are routed
        # through bwa_mem_paf_cmd() instead (see TestBwaMemPafCommand).
        for mapper, coverm_name in aviary_init.SHORT_READ_MAPPER_TO_COVERM.items():
            if mapper in MAPPERS_WITHOUT_PAF_COMMAND:
                continue
            command = polish.short_read_paf_cmd(
                coverm_name, "ref.fa", "r.fq.gz", 4)
            self.assertTrue(command.startswith(coverm_name.split('-')[0]))
            self.assertIn("-t 4", command)

    def test_unknown_mapper_raises(self):
        # Fail loudly rather than falling back to a default, which would
        # silently polish against the wrong aligner. bwa-mem is a real mapper
        # but has no PAF mode, so short_read_paf_cmd specifically must not
        # know it -- callers are expected to route it to bwa_mem_paf_cmd
        # instead (see run_polish()).
        with self.assertRaises(ValueError):
            polish.short_read_paf_cmd("bwa-mem", "ref.fa", "r.fq.gz", 8)
        with self.assertRaises(ValueError):
            polish.short_read_paf_cmd("bwa-mem2", "ref.fa", "r.fq.gz", 8)


class TestMinibwaLongReadPafCommand(unittest.TestCase):
    """minibwa's long-read PAF command, used by run_polish() when
    --long-read-mapper minibwa. Mirrors short_read_paf_cmd()'s minibwa branch
    (same 'map' subcommand, same -f) but with -x lr instead of -x sr -- -x
    belongs to minibwa map's own option loop, not a top-level flag, so
    `minibwa -x lr ...` (the minimap2 calling convention) is not equivalent
    and fails with "unknown command '-x'".
    """

    def test_minibwa_uses_map_subcommand_f_and_lr_preset(self):
        self.assertEqual(
            polish.minibwa_long_read_paf_cmd("ref.fa", "r.fq.gz", 8),
            "minibwa map -f -x lr -t 8 ref.fa r.fq.gz")


class TestBwaMemPafCommand(unittest.TestCase):
    """bwa mem/bwa-mem2 mem only emit SAM, so racon PAF generation needs an
    index step and a sam2paf conversion that short_read_paf_cmd's
    single-command model can't express.
    """

    def _run(self, mapper, paired_interleaved=False):
        with tempfile.TemporaryDirectory() as tmp:
            output_paf = os.path.join(tmp, "out.paf")
            with mock.patch.object(polish, "run") as mocked_run, \
                 mock.patch.object(polish, "Popen") as mocked_popen:
                mocked_run.return_value = mock.Mock(returncode=0)
                mem_proc = mock.Mock(stdout=mock.Mock(), returncode=0)
                sam2paf_proc = mock.Mock(returncode=0)
                mocked_popen.side_effect = [mem_proc, sam2paf_proc]
                polish.bwa_mem_paf_cmd(
                    mapper, "ref.fa", "r1.fq.gz r2.fq.gz", 8, output_paf, "/dev/null",
                    paired_interleaved=paired_interleaved)
                return mocked_run, mocked_popen

    def test_bwa_mem_indexes_before_mapping(self):
        mocked_run, mocked_popen = self._run("bwa-mem")
        index_cmd = mocked_run.call_args[0][0]
        self.assertEqual(index_cmd, ["bwa", "index", "ref.fa"])
        mem_cmd = mocked_popen.call_args_list[0].args[0]
        self.assertEqual(mem_cmd[:2], ["bwa", "mem"])

    def test_bwa_mem2_uses_its_own_binary_name(self):
        mocked_run, mocked_popen = self._run("bwa-mem2")
        index_cmd = mocked_run.call_args[0][0]
        self.assertEqual(index_cmd, ["bwa-mem2", "index", "ref.fa"])
        mem_cmd = mocked_popen.call_args_list[0].args[0]
        self.assertEqual(mem_cmd[:2], ["bwa-mem2", "mem"])

    def test_conversion_step_uses_sam2paf(self):
        _, mocked_popen = self._run("bwa-mem")
        sam2paf_cmd = mocked_popen.call_args_list[1].args[0]
        self.assertEqual(sam2paf_cmd[:2], ["paftools.js", "sam2paf"])

    def test_index_failure_raises(self):
        with mock.patch.object(polish, "run") as mocked_run:
            mocked_run.return_value = mock.Mock(returncode=1)
            with self.assertRaises(RuntimeError):
                polish.bwa_mem_paf_cmd("bwa-mem", "ref.fa", "r.fq.gz", 8, "out.paf", "/dev/null")

    def test_paired_interleaved_passes_p_to_bwa_mem(self):
        # Regression test: verified empirically against a real bwa binary
        # that without -p, bwa mem treats a genuinely interleaved input as
        # unpaired single-end reads and never sets FLAG 0x40/0x80 at all
        # (FLAG comes back 0 for every record), which makes sam2paf's own
        # mate-suffix restoration silently do nothing and brings back the
        # exact empty-read-set bug this whole mechanism exists to fix. -p is
        # what makes bwa mem actually set those FLAG bits.
        _, mocked_popen = self._run("bwa-mem", paired_interleaved=True)
        mem_cmd = mocked_popen.call_args_list[0].args[0]
        self.assertIn("-p", mem_cmd)

    def test_non_interleaved_reads_do_not_get_p(self):
        # run_polish()'s other reads source (paired R1/R2 concatenated
        # block-then-block via clean_short_reads, not interleaved) must never
        # get -p: it would wrongly pair up two R1 reads as mates.
        _, mocked_popen = self._run("bwa-mem", paired_interleaved=False)
        mem_cmd = mocked_popen.call_args_list[0].args[0]
        self.assertNotIn("-p", mem_cmd)

    def test_mem_failure_raises(self):
        with tempfile.TemporaryDirectory() as tmp, \
             mock.patch.object(polish, "run") as mocked_run, \
             mock.patch.object(polish, "Popen") as mocked_popen:
            mocked_run.return_value = mock.Mock(returncode=0)
            mem_proc = mock.Mock(stdout=mock.Mock(), returncode=1)
            sam2paf_proc = mock.Mock(returncode=0)
            mocked_popen.side_effect = [mem_proc, sam2paf_proc]
            with self.assertRaises(RuntimeError):
                polish.bwa_mem_paf_cmd(
                    "bwa-mem", "ref.fa", "r.fq.gz", 8, os.path.join(tmp, "out.paf"), "/dev/null")


class TestLongReadMapperGuard(unittest.TestCase):
    """get_abundances.run_coverm refuses short-read mappers for long reads.

    CoverM defaults to strobealign, which is short-read only. Long reads mapped
    with it do not necessarily error -- they can yield a well-formed abundance
    table of near-zero values, so the mistake would surface as wrong numbers
    rather than a crash.
    """

    BASE = dict(reads="r.fq", output_file="o.tsv", read_type="--single",
                threads=1, strain_analysis=False, output_dir="d/",
                log="/dev/null", bins_dir="b/")

    def test_rejects_short_read_mappers_for_long_reads(self):
        # minibwa is deliberately excluded here: unlike minimap2-sr/rammap-sr,
        # bare "minibwa" IS the valid long-read value too (it has no preset
        # suffix at all), steered instead via --minibwa-params.
        for bad in (None, "strobealign", "minimap2-sr", "rammap-sr"):
            with self.assertRaises(ValueError):
                get_abundances.run_coverm(
                    minimap2_type=bad, long_reads=True, **self.BASE)

    def test_accepts_long_read_mappers(self):
        for good in ("minibwa", "minimap2-ont", "minimap2-pb", "minimap2-lr-hq", "minimap2-hifi",
                     "rammap-ont", "rammap-pb", "rammap-lr-hq", "rammap-hifi"):
            self.assertIn(good, get_abundances.LONG_READ_MAPPERS)

    def test_minibwa_accepted_for_long_reads_without_a_preset_suffix(self):
        # minibwa has no CoverM preset suffix for long reads; passing bare
        # "minibwa" must not be rejected the way other bare short-read mapper
        # names are.
        try:
            get_abundances.run_coverm(
                minimap2_type="minibwa", long_reads=True, **self.BASE)
        except ValueError as exception:
            self.fail(f"minibwa wrongly rejected for long reads: {exception}")
        except (FileNotFoundError, OSError):
            pass  # reached the CoverM exec, which is as far as we need to get

    def test_short_read_calls_are_unaffected(self):
        # long_reads defaults to False, so short-read calls must not be gated.
        try:
            get_abundances.run_coverm(minimap2_type="strobealign", **self.BASE)
        except ValueError as exception:
            self.fail(f"short-read path wrongly blocked: {exception}")
        except (FileNotFoundError, OSError):
            pass  # reached the CoverM exec, which is as far as we need to get


class TestAlignerPresetTranslation(unittest.TestCase):
    def test_no_preset_emits_no_x_flag(self):
        self.assertEqual(aviary_init.aligner_preset_args("no-preset"), [])

    def test_every_model_has_a_translation(self):
        # A model in MAPPER_MODELS without a translation would build an
        # aligner command line the aligner rejects.
        for models in aviary_init.MAPPER_MODELS.values():
            for model in models:
                self.assertIn(model, aviary_init.COVERM_MODEL_TO_ALIGNER_PRESET)

    def test_models_are_partitioned_by_read_length(self):
        # Every selectable model must be valid for exactly one read length,
        # except no-preset which is deliberately valid for both.
        for models in aviary_init.MAPPER_MODELS.values():
            for model in models:
                in_short = model in aviary_init.SHORT_READ_MODELS
                in_long = model in aviary_init.LONG_READ_MODELS
                if model == "no-preset":
                    self.assertTrue(in_short and in_long)
                else:
                    self.assertNotEqual(in_short, in_long, f"{model!r} misclassified")

    def test_unknown_model_raises(self):
        with self.assertRaises(ValueError):
            aviary_init.aligner_preset_args("bogus")


class TestMapperCliValidation(unittest.TestCase):
    """End-to-end CLI guards: a mapper named without the reads it would map,
    or with a preset for the wrong read length, must fail before the run
    starts rather than partway through a mapping job.
    """

    DATA = os.path.join(REPO_ROOT, "test", "data")

    def _run(self, extra_args):
        import subprocess
        env = dict(os.environ, GTDBTK_DATA_PATH=".", CHECKM2DB=".",
                   EGGNOG_DATA_DIR=".", METABULI_DB_PATH=".",
                   SINGLEM_METAPACKAGE_PATH=".")
        with tempfile.TemporaryDirectory() as tmpdir:
            return subprocess.run(
                ["aviary", "recover",
                 "--assembly", os.path.join(self.DATA, "assembly.fasta"),
                 "--output", os.path.join(tmpdir, "out"),
                 "--dryrun", "--tmpdir", tmpdir] + extra_args,
                capture_output=True, text=True, env=env)

    def test_short_read_mapper_without_short_reads_is_rejected(self):
        result = self._run(["--short-read-mapper", "bwa-mem"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no reads were supplied", result.stdout + result.stderr)

    def test_long_read_mapper_without_long_reads_is_rejected(self):
        result = self._run(["--long-read-mapper", "minimap2"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no reads were supplied", result.stdout + result.stderr)

    def test_long_read_preset_for_short_reads_is_rejected(self):
        result = self._run([
            "-1", os.path.join(self.DATA, "wgsim.1.fq.gz"),
            "-2", os.path.join(self.DATA, "wgsim.2.fq.gz"),
            "--short-read-mapper", "minimap2",
            "--short-read-mapper-model", "ont"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not a short-read preset", result.stdout + result.stderr)

    def test_matching_mapper_and_reads_is_accepted(self):
        result = self._run([
            "-1", os.path.join(self.DATA, "wgsim.1.fq.gz"),
            "-2", os.path.join(self.DATA, "wgsim.2.fq.gz"),
            "--short-read-mapper", "bwa-mem"])
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_defaults_need_no_flags_and_stay_strobealign_and_rammap(self):
        # The argparse default is None purely so "given" is distinguishable
        # from "unset"; the resolved defaults must not change.
        self.assertEqual(aviary_init.DEFAULT_SHORT_READ_MAPPER, "strobealign")
        self.assertEqual(aviary_init.DEFAULT_LONG_READ_MAPPER, "rammap")
        result = self._run([
            "-1", os.path.join(self.DATA, "wgsim.1.fq.gz"),
            "-2", os.path.join(self.DATA, "wgsim.2.fq.gz")])
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == '__main__':
    unittest.main()
