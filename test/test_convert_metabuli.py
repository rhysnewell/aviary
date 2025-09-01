#!/usr/bin/env python3

import unittest
from aviary.modules.binning.scripts.convert_metabuli import create_conversion_dict, process_classifications
import pandas as pd
import pandas.testing as pdt

class Tests(unittest.TestCase):
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
                None,
                "d__Eukaryota,p__Chordata,c__Mammalia,o__Primates,f__Hominidae,g__Homo,s__Homo sapiens",
                "d__Eukaryota,p__Chordata,c__Mammalia,o__Primates,f__Hominidae,g__Homo,s__Homo sapiens",
                "d__Eukaryota",
                "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Dermatophilaceae,g__Allobranchiibius,s__Allobranchiibius huperziae",
                "d__Bacteria,p__Actinomycetota,c__Actinomycetia,o__Actinomycetales,f__Microbacteriaceae,g__Rhodoglobus,s__Rhodoglobus sp004118935"
                ],
        })

        observed = process_classifications(df, conversion_dict, coverm_out)
        pdt.assert_frame_equal(expected, observed)


if __name__ == '__main__':
    unittest.main()
