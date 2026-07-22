#!/usr/bin/env python3

import unittest
import os
import tempfile
from aviary.modules.binning.scripts.convert_metabuli import (
    create_conversion_dict,
    process_classifications,
    read_classifications,
)
import pandas as pd
import pandas.testing as pdt
import numpy as np

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

# These are fast, self-contained unit tests of convert_metabuli's parsing --
# no aviary run, no cluster, no database. They deliberately do NOT write into
# example/: that tree holds real pipeline output, and fabricating an
# aviary_out/logs/... layout here would imply a run happened when none did.
#
# End-to-end coverage of the same bug lives in
# test_integration.py::test_short_read_recovery_metabuli_unclassified, which
# runs the actual pipeline over unclassifiable contigs and produces a genuine
# example/ tree. These remain alongside it because they are ~0.5s rather than
# a 150GB/32-core job, and because they do not depend on Metabuli's behaviour:
# if that integration test's assumption (random sequence stays unclassified)
# ever breaks, its precondition assertion fires and it becomes inconclusive --
# these still guard the parser.


class Tests(unittest.TestCase):
    def test_read_classifications_real_fixture_with_unclassified_rows(self):
        # Same defect as the constructed fixture below, but exercised against a
        # captured sample of REAL Metabuli output (test/data/, 30 rows taken
        # verbatim from a taxvamb run) with 5 rows rewritten into Metabuli's
        # unclassified shape. Kept alongside the constructed one because they
        # fail differently if they ever diverge: this catches "our synthetic
        # bytes drifted from what Metabuli actually emits", the constructed one
        # documents the exact shape being defended against.
        #
        # The example/ test outputs cannot serve this purpose -- they are
        # gitignored generated output, regenerated on every run, and every
        # metabuli-producing test (queue_submission_gpus, taxvamb,
        # taxvamb_gpu) shares one fixture in which all 7070 reads classify,
        # i.e. 0 unclassified rows. That is why this bug shipped green.
        path = os.path.join(DATA_DIR, "metabuli_tax_classifications_with_unclassified.tsv")

        df = read_classifications(path)

        # Header row + 25 classified + 5 unclassified, uniform 8 columns.
        self.assertEqual(df.shape, (31, 8))
        self.assertEqual((df[0] == "0").sum(), 5)
        self.assertEqual((df[0] == "1").sum(), 25)

        # Unclassified rows keep their real contig name (col 1) and taxID 0
        # (col 2) -- i.e. the trailing tab was absorbed rather than shifting
        # the columns rightward.
        unclassified = df[df[0] == "0"]
        self.assertTrue(all(n.startswith("NODE_") for n in unclassified[1]))
        self.assertEqual(set(unclassified[2]), {"0"})

    def test_read_classifications_handles_unclassified_trailing_tab(self):
        # Reproduces the real Metabuli tax_classifications.tsv shape that
        # crashed convert_metabuli: unclassified (is_classified=0) rows carry a
        # stray trailing tab, giving them 9 tab-fields, while classified (1)
        # rows have 8. A plain pd.read_csv(header=None) raises
        # "Expected 8 fields, saw 9" on the first 0-row. read_classifications
        # must parse both row types into a uniform 8-column frame.
        #
        # The prior test fixture (test_short_read_recovery...gpus) never
        # exercised this because it happened to contain only is_classified=1
        # rows -- so the 0-row path was completely untested until this.
        lines = [
            "#is_classified\tname\ttaxID\tquery_length\tscore\te_value\trank\ttaxID:match_count",
            "1\tNODE_1_len\t370388\t138942\t0.00113\t1.2e+16\tspecies\t10051421:75 ",
            "1\tNODE_2_len\t1\t134220\t0.00086\t1.1e+16\tno rank\t",       # empty match_count, 8 fields
            "0\tk141_54460\t0\t1701\t0\t-\t-\t-\t",                        # trailing tab -> 9 raw fields
            "0\tk141_239609\t0\t2124\t0\t-\t-\t-\t",                       # trailing tab -> 9 raw fields
        ]
        fd, path = tempfile.mkstemp(suffix=".tsv")
        try:
            with os.fdopen(fd, "w") as fh:
                fh.write("\n".join(lines) + "\n")

            df = read_classifications(path)
        finally:
            os.remove(path)

        # 8 columns for every row, header + 2 classified + 2 unclassified
        self.assertEqual(df.shape, (5, 8))
        self.assertEqual(list(df.columns), list(range(8)))

        # The two unclassified rows parsed with their real 8 fields intact
        # (name in col 1, taxID 0 in col 2), the phantom 9th field absorbed.
        unclassified = df[df[0] == "0"]
        self.assertEqual(len(unclassified), 2)
        self.assertEqual(unclassified[1].tolist(), ["k141_54460", "k141_239609"])
        self.assertEqual(unclassified[2].tolist(), ["0", "0"])

        # Classified rows still read correctly alongside them.
        classified = df[df[0] == "1"]
        self.assertEqual(classified[1].tolist(), ["NODE_1_len", "NODE_2_len"])

    def test_create_conversion_dict(self):
        ids = [
            0,
            1,
            14524,
            101055,
            101056,
            171647,
            219704,
            260237,
            309837,
            302751,
            384239,
            378555,
            491973,
            171648,
            171649,
            171650,
            187121,
            304979,
            239098,
            346849,
            496373,
            243063,
            338255,
            338256,
            338257,
            10000000,
            10000001,
            10000002,
            10000003,
            10000004,
            10000005,
            10000006,
        ]
        labels = [
            "unclassified",
            "root",
            "  d__Bacteria",
            "    p__Actinomycetota",
            "      c__Actinomycetia",
            "        o__Actinomycetales",
            "          f__Microbacteriaceae",
            "            g__Microbacterium",
            "              s__Microbacterium aurum",
            "            g__Rhodoglobus",
            "              s__Rhodoglobus sp004118935",
            "            g__Gryllotalpicola",
            "              s__Gryllotalpicola protaetiae",
            "          f__Bifidobacteriaceae",
            "            g__Bifidobacterium",
            "              s__Bifidobacterium longum",
            "          f__Micrococcaceae",
            "            g__Arthrobacter",
            "          f__Actinomycetaceae",
            "            g__Flaviflexus",
            "              s__Flaviflexus salsibiostraticola",
            "          f__Dermatophilaceae",
            "            g__Allobranchiibius",
            "              s__Allobranchiibius huperziae",
            "                RS_GCF_021391335.1",
            "  d__Eukaryota",
            "    p__Chordata",
            "      c__Mammalia",
            "        o__Primates",
            "          f__Hominidae",
            "            g__Homo",
            "              s__Homo sapiens",
        ]


        expected = {
            14524: "d__Bacteria",
            101055: "d__Bacteria,p__Actinomycetota",
            101056: "d__Bacteria,p__Actinomycetota,c__Actinomycetia",
            171647: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales",
            219704: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Microbacteriaceae",
            260237: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Microbacteriaceae,g__Microbacterium",
            309837: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Microbacteriaceae,g__Microbacterium,s__Microbacterium aurum",
            302751: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Microbacteriaceae,g__Rhodoglobus",
            384239: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Microbacteriaceae,g__Rhodoglobus,s__Rhodoglobus sp004118935",
            378555: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Microbacteriaceae,g__Gryllotalpicola",
            491973: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Microbacteriaceae,g__Gryllotalpicola,s__Gryllotalpicola protaetiae",
            171648: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Bifidobacteriaceae",
            171649: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Bifidobacteriaceae,g__Bifidobacterium",
            171650: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Bifidobacteriaceae,g__Bifidobacterium,s__Bifidobacterium longum",
            187121: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Micrococcaceae",
            304979: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Micrococcaceae,g__Arthrobacter",
            239098: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Actinomycetaceae",
            346849: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Actinomycetaceae,g__Flaviflexus",
            496373: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Actinomycetaceae,g__Flaviflexus,s__Flaviflexus salsibiostraticola",
            243063: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Dermatophilaceae",
            338255: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Dermatophilaceae,g__Allobranchiibius",
            338256: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Dermatophilaceae,g__Allobranchiibius,s__Allobranchiibius huperziae",
            338257: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Dermatophilaceae,g__Allobranchiibius,s__Allobranchiibius huperziae",
            10000000: "d__Eukaryota",
            10000001: "d__Eukaryota,p__Chordata",
            10000002: "d__Eukaryota,p__Chordata,c__Mammalia",
            10000003: "d__Eukaryota,p__Chordata,c__Mammalia,o__Primates",
            10000004: "d__Eukaryota,p__Chordata,c__Mammalia,o__Primates,f__Hominidae",
            10000005: "d__Eukaryota,p__Chordata,c__Mammalia,o__Primates,f__Hominidae,g__Homo",
            10000006: "d__Eukaryota,p__Chordata,c__Mammalia,o__Primates,f__Hominidae,g__Homo,s__Homo sapiens",
        }

        observed = create_conversion_dict(ids, labels)
        self.assertEqual(expected, observed)

    def test_process_classifications(self):
        df = pd.DataFrame({
            0: [1,1,1,1,1,1,1,1,1],
            1: [
                "NODE_1_length_138944_cov_12.417585",
                "NODE_2_length_134223_cov_12.412460",
                "NODE_3_length_120075_cov_12.461090",
                "NODE_4_length_113818_cov_12.404648",
                "NODE_5_length_74510_cov_12.470123",
                "NODE_6_length_71626_cov_5.654567",
                "NODE_7_length_69828_cov_12.473579",
                "NODE_7_length_1498_cov_12.473579",
                "NODE_7_length_275_cov_12.473579",
                ],
            2: [101055,1,10000006,10000006,10000000,338257,384239,1,1],
            3: [138942,134220,120072,113814,74508,71622,69825,69825,69825],
            4: [0.00115876,0.000865494,0.00132837,0.00104996,0.00171123,0.00321828,0.00134622,0.00134622,0.00134622],
            5: ["species","no rank","species","species","species","species","species","species","species"],
            6: ["10051421:75 ","","10000006:1109 10000007:915 ","10000006:1160 10000007:904 ","260438:20 10031517:1 10031520:12 ","48381:203 10006244:1 10006245:1 ","10024879:37 ","10024879:37 ","10024879:37 "],
        })
        conversion_dict = {
            101055: "d__Bacteria,p__Actinomycetota",
            384239: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Microbacteriaceae,g__Rhodoglobus,s__Rhodoglobus sp004118935",
            338257: "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Dermatophilaceae,g__Allobranchiibius,s__Allobranchiibius huperziae",
            10000000: "d__Eukaryota",
            10000006: "d__Eukaryota,p__Chordata,c__Mammalia,o__Primates,f__Hominidae,g__Homo,s__Homo sapiens",
        }
        coverm_out = pd.DataFrame({
            "contigs": [
                "NODE_1_length_138944_cov_12.417585",
                "NODE_2_length_134223_cov_12.412460",
                "NODE_3_length_120075_cov_12.461090",
                "NODE_4_length_113818_cov_12.404648",
                "NODE_5_length_74510_cov_12.470123",
                "NODE_6_length_71626_cov_5.654567",
                "NODE_7_length_69828_cov_12.473579",
                ]
        })

        expected = pd.DataFrame({
            "contigs": [
                "NODE_1_length_138944_cov_12.417585",
                "NODE_2_length_134223_cov_12.412460",
                "NODE_3_length_120075_cov_12.461090",
                "NODE_4_length_113818_cov_12.404648",
                "NODE_5_length_74510_cov_12.470123",
                "NODE_6_length_71626_cov_5.654567",
                "NODE_7_length_69828_cov_12.473579"
                ],
            "predictions": [
                "d__Bacteria,p__Actinomycetota",
                np.nan,
                "d__Eukaryota,p__Chordata,c__Mammalia,o__Primates,f__Hominidae,g__Homo,s__Homo sapiens",
                "d__Eukaryota,p__Chordata,c__Mammalia,o__Primates,f__Hominidae,g__Homo,s__Homo sapiens",
                "d__Eukaryota",
                "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Dermatophilaceae,g__Allobranchiibius,s__Allobranchiibius huperziae",
                "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Microbacteriaceae,g__Rhodoglobus,s__Rhodoglobus sp004118935"
                ],
        })

        observed = process_classifications(df, conversion_dict, coverm_out)
        # Normalize column index dtype across pandas versions (string vs object)
        expected.columns = expected.columns.astype(object)
        observed.columns = observed.columns.astype(object)
        pdt.assert_frame_equal(expected, observed)


if __name__ == '__main__':
    unittest.main()
