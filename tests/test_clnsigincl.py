"""CLNSIGINCL parsing: '|' separates pairs and also the terms of one classification.
Every case is checked in the three paths that read the field, so that they cannot diverge:
the INFO check (step1_index.check_info), the exposure rule (clinvar_records.vcf_seen_strict) and the
record built from a VCF line (strict token and unknown-term collection)."""
import os
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
os.environ["AVI_EVAL_SYNTHETIC"] = "1"

import clinvar_records as cr  # noqa: E402
import clinvar_index as idx  # noqa: E402


def record(clnsigincl: str, clnsig: str = "."):
    (rec,) = cr.from_vcf_line(f"3\t100\t1\tA\tG\t.\t.\tALLELEID=11;CLNSIG={clnsig};CLNSIGINCL={clnsigincl}", "E-dev")
    return rec


CASES = {
    # value: (pairs, seen, strict_token, unknown_terms)
    "8152:Pathogenic|protective|risk_factor":                                   # a real line of the 2025-06-15 VCF
        ([("8152", "Pathogenic|protective|risk_factor")], True, True, ()),
    "8152:Uncertain_significance|Benigne":                                      # unknown term in a continuation
        ([("8152", "Uncertain_significance|Benigne")], True, False, ("benigne",)),
    "1:Benign|2:Pathogenic|protective":                                         # continuation after the last pair
        ([("1", "Benign"), ("2", "Pathogenic|protective")], True, True, ()),
    "1:Uncertain_significance|risk_factor|2:not_provided|protective":           # two pairs, each with continuations
        ([("1", "Uncertain_significance|risk_factor"), ("2", "not_provided|protective")], False, False, ()),
    "424752:Pathogenic|424753:Benign":                                          # ordinary pairs
        ([("424752", "Pathogenic"), ("424753", "Benign")], True, True, ()),
}
MALFORMED = ("protective|1:Benign", "1:", "1:Benign|", "1:Benign||2:Pathogenic", "x:Benign", ":Benign",
             # a colon-bearing token that is not a valid pair is malformed in any position
             "1:Uncertain_significance|x:Benign", "1:Uncertain_significance|:Benign",
             "1:Uncertain_significance|risk_factor:x")


class Parser(unittest.TestCase):
    def test_pairs(self):
        for value, (pairs, *_rest) in CASES.items():
            with self.subTest(value=value):
                self.assertEqual(cr.incl_pairs(value), pairs)

    def test_malformed_values(self):
        for value in MALFORMED:
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    cr.incl_pairs(value)


class ThreePaths(unittest.TestCase):
    def test_info_check(self):
        for value in CASES:
            with self.subTest(value=value):
                idx.check_info(f"ALLELEID=11;CLNSIGINCL={value}", "E-dev", 1)
        for value in MALFORMED:
            with self.subTest(value=value):
                with self.assertRaises(idx.InputFormatError):
                    idx.check_info(f"ALLELEID=11;CLNSIGINCL={value}", "E-dev", 1)

    def test_exposure_rule(self):
        for value, (_pairs, seen, *_rest) in CASES.items():
            with self.subTest(value=value):
                self.assertEqual(cr.vcf_seen_strict(".", value), seen)

    def test_malformed_values_are_refused_by_the_exposure_rule_and_the_record(self):
        for value in MALFORMED:
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    cr.vcf_seen_strict(".", value)
                with self.assertRaises(ValueError):
                    record(value)

    def test_record_from_vcf_line(self):
        for value, (_pairs, seen, strict, unknown) in CASES.items():
            with self.subTest(value=value):
                rec = record(value)
                self.assertEqual((rec.seen, rec.strict_token, rec.unknown_terms), (seen, strict, unknown))


if __name__ == "__main__":
    unittest.main()
