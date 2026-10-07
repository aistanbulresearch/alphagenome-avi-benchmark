"""Score parsers and the lookup over a synthetic tabix file (no real scores)."""
import sys
import unittest
import unittest.mock
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import gnomad_remote as gr  # noqa: E402
import tool_scores as ts  # noqa: E402
from test_gnomad_remote import build  # noqa: E402


class Parsers(unittest.TestCase):
    def test_spliceai_max_over_scores_and_entries(self):
        info = "SpliceAI=G|GENE1|0.10|0.00|0.35|0.02|-3|10|2|-40,G|GENE2|0.50|0.01|0.00|0.00|1|2|3|4,T|GENE1|0.9|0|0|0|1|1|1|1"
        self.assertEqual(ts.spliceai_max(info, "G"), "0.5")
        self.assertEqual(ts.spliceai_max(info, "C"), "")
        parts = ["1", "100", ".", "A", "G", ".", ".", info]
        self.assertEqual(ts.parse_spliceai(parts, "A", "G"), {"spliceai_max": "0.5"})
        self.assertIsNone(ts.parse_spliceai(parts, "C", "G"))

    def test_cadd_columns(self):
        parts = ["1", "100", "A", "G", "0.523", "7.812\n"]
        self.assertEqual(ts.parse_cadd(parts, "A", "G"), {"cadd_raw": "0.523", "cadd_phred": "7.812"})
        self.assertIsNone(ts.parse_cadd(parts, "A", "T"))

    def test_avi_columns_found_by_header_name(self):
        parse = ts.parse_avi(["#chrom", "pos", "ref", "alt", "avi_score", "avi_phred"])
        self.assertEqual(parse(["chr1", "5", "C", "T", "0.91", "23.4"], "C", "T"), {"avi_raw": "0.91", "avi_phred": "23.4"})
        self.assertIsNone(parse(["chr1", "5", "C", "T", "0.91", "23.4"], "C", "A"))
        with self.assertRaises(ValueError):
            ts.parse_avi(["chrom", "pos", "ref", "alt", "x"])


class Lookup(unittest.TestCase):
    def test_lookup_over_synthetic_tabix(self):
        data, tbi, _ = build(500)
        index = gr.parse_tbi(tbi)
        fetch = lambda url, a, b: data[a:b]  # noqa: E731
        parse = lambda parts, ref, alt: {"x": parts[4]} if parts[3] == ref and parts[4] == alt else None  # noqa: E731
        variants = [{"chrom": "21", "pos": "70", "ref": "A", "alt": "T", "variation_id": "1"},
                    {"chrom": "21", "pos": "70", "ref": "A", "alt": "C", "variation_id": "2"},
                    {"chrom": "21", "pos": "75", "ref": "A", "alt": "G", "variation_id": "3"},
                    {"chrom": "22", "pos": "10", "ref": "A", "alt": "G", "variation_id": "4"}]
        got = {v["variation_id"]: h for v, h in ts.lookup_all(variants, "u", index, lambda c: f"chr{c}", parse, fetch, 2)}
        self.assertEqual(got, {"1": {"x": "T"}, "2": None, "3": None, "4": None})
        variants.append({"chrom": "21", "pos": "4000", "ref": "A", "alt": "G", "variation_id": "5"})
        for limit in (0, ts.WHOLE_CHUNK):                     # searched one by one, or one whole-chunk read
            with unittest.mock.patch.object(ts, "WHOLE_CHUNK", limit):
                got = {v["variation_id"]: h for v, h in ts.lookup_all(variants, "u", index, lambda c: f"chr{c}",
                                                                      parse, fetch, 2)}
            self.assertEqual(got, {"1": {"x": "T"}, "2": None, "3": None, "4": None, "5": {"x": "G"}})


if __name__ == "__main__":
    unittest.main()
