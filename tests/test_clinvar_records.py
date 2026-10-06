import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import clinvar_records as cr  # noqa: E402


def vs_row(**over):
    row = {"#AlleleID": "15000", "Type": "single nucleotide variant",
           "Name": "NM_000001.1(GENE1):c.100+60A>G", "GeneID": "101", "GeneSymbol": "GENE1",
           "ClinicalSignificance": "Pathogenic", "ReviewStatus": "criteria provided, single submitter",
           "Assembly": "GRCh38", "Chromosome": "3", "VariationID": "9000", "PositionVCF": "12345",
           "ReferenceAlleleVCF": "A", "AlternateAlleleVCF": "G"}
    row.update(over)
    return row


class Classification(unittest.TestCase):
    def test_strict_eligibility_is_exact(self):
        for value, group in [("Pathogenic", "P/LP"), ("Likely pathogenic", "P/LP"),
                             ("Pathogenic/Likely pathogenic", "P/LP"), ("Benign", "B/LB"),
                             ("Likely benign", "B/LB"), ("Benign/Likely benign", "B/LB"),
                             ("Likely pathogenic, low penetrance", None), ("Pathogenic; risk factor", None),
                             ("pathogenic", None), ("Conflicting classifications of pathogenicity", None),
                             ("Uncertain significance", None), ("", None), (None, None)]:
            with self.subTest(value=value):
                self.assertEqual(cr.strict_group(value), group)

    def test_exposure_tokens(self):
        for value, seen in [("Likely pathogenic, low penetrance", True), ("Pathogenic; risk factor", True),
                            ("Pathogenic/Likely pathogenic/Pathogenic, low penetrance", True),
                            ("Benign/Likely benign|other", True), ("risk factor", False),
                            ("Conflicting classifications of pathogenicity", False),
                            ("Uncertain significance", False), ("Likely risk allele", False),
                            ("Pathogenic, low penetrance", True), ("", False), (None, False)]:
            with self.subTest(value=value):
                self.assertEqual(cr.seen_strict(value), seen)

    def test_unknown_tokens_count_as_seen(self):
        self.assertTrue(cr.seen_strict("Pathogneic"))                 # a misspelling is unknown
        self.assertEqual(cr.unknown_tokens("Pathogneic; risk factor"), ["pathogneic"])
        self.assertFalse(cr.seen_strict("Uncertain significance; risk factor"))
        self.assertFalse(cr.seen_strict("-"))
        self.assertFalse(cr.seen_strict("."))                         # no classification
        self.assertTrue(cr.strict_present("Likely_pathogenic,_low_penetrance"))
        self.assertFalse(cr.strict_present("Pathogneic"))
        self.assertTrue(cr.vcf_seen_strict("Uncertain_significance", "424752:Benigne"))

    def test_vcf_tokens(self):
        self.assertTrue(cr.vcf_seen_strict("Likely_pathogenic,_low_penetrance", None))
        self.assertTrue(cr.vcf_seen_strict("Pathogenic/Likely_pathogenic|_risk_factor", ""))
        self.assertFalse(cr.vcf_seen_strict("Uncertain_significance", None))
        self.assertTrue(cr.vcf_seen_strict("Uncertain_significance", "424752:Pathogenic|424753:Benign"))
        self.assertFalse(cr.vcf_seen_strict("", "424752:Uncertain_significance"))
        self.assertFalse(cr.vcf_seen_strict("Conflicting_classifications_of_pathogenicity", None))


class StarMapping(unittest.TestCase):
    def test_every_listed_string(self):
        for text, value in cr.STAR_MAP.items():
            with self.subTest(text=text):
                self.assertEqual(cr.stars(text), value)
                self.assertEqual(cr.stars(text.replace(" ", "_")), value)

    def test_unknown_strings(self):
        for text in ["", None, "criteria provided", "Reviewed by expert panel", "four stars"]:
            with self.subTest(text=text):
                self.assertEqual(cr.stars(text), cr.UNKNOWN)


class Chromosomes(unittest.TestCase):
    def test_groups_cover_autosomes_and_x(self):
        all_chroms = [c for chroms in cr.CHROM_GROUPS.values() for c in chroms]
        self.assertEqual(sorted(all_chroms), sorted([str(i) for i in range(1, 23)] + ["X"]))

    def test_group_and_scopes(self):
        self.assertEqual(cr.chromosome_group("3"), "untouched")
        self.assertEqual(cr.chromosome_group("chr2"), "validation")
        self.assertEqual(cr.chromosome_group("x"), "validation")
        self.assertEqual(cr.chromosome_group("1"), "training")
        for c in ["Y", "MT", "Un", "", None, "NW_123"]:
            self.assertEqual(cr.chromosome_group(c), "excluded")
        self.assertTrue(cr.in_scope("untouched", "S1"))
        self.assertFalse(cr.in_scope("validation", "S1"))
        self.assertTrue(cr.in_scope("validation", "S2"))
        self.assertFalse(cr.in_scope("training", "S2"))
        self.assertTrue(cr.in_scope("training", "S3"))
        self.assertFalse(cr.in_scope("excluded", "S3"))


class Alleles(unittest.TestCase):
    def test_normalization(self):
        self.assertEqual(cr.normalize_allele("chr3", "100", "a", "g"), "3:100:A:G")
        for chrom in ("X", "x", "chrX", "chrx", "CHRX"):
            self.assertEqual(cr.normalize_allele(chrom, 100, "A", "G"), "X:100:A:G")
        self.assertEqual(cr.normalize_allele("3", 100, "CAT", "CGT"), "3:101:A:G")
        self.assertEqual(cr.normalize_allele("3", 100, "CA", "C"), "3:100:CA:C")
        self.assertEqual(cr.normalize_allele("3", 100, "CAG", "CG"), "3:100:CA:C")
        for args in [("3", -1, "A", "G"), ("3", "na", "A", "G"), ("3", 100, "na", "G"),
                     ("3", 100, "A", "A"), ("3", 100, "-", "G"), ("", 100, "A", "G")]:
            self.assertIsNone(cr.normalize_allele(*args))

    def test_variant_type(self):
        self.assertEqual(cr.variant_type("3:100:A:G"), "snv")
        self.assertEqual(cr.variant_type("3:100:CA:C"), "small_indel")
        self.assertEqual(cr.variant_type("3:100:C:C" + "A" * 10), "small_indel")
        self.assertEqual(cr.variant_type("3:100:C:C" + "A" * 11), "other")
        self.assertEqual(cr.variant_type("3:100:AC:GT"), "other")
        self.assertEqual(cr.variant_type("3:100:ACG:T"), "other")
        self.assertEqual(cr.variant_type(None), "other")


class Records(unittest.TestCase):
    def test_variant_summary_grch38(self):
        r = cr.from_variant_summary(vs_row(), "T1")
        self.assertEqual((r.key, r.allele_id, r.variation_id, r.stars, r.seen, r.gene),
                         ("3:12345:A:G", "15000", "9000", 1, True, "101"))

    def test_variant_summary_grch37_has_key37_only(self):
        r = cr.from_variant_summary(vs_row(Assembly="GRCh37"), "T1")
        self.assertIsNone(r.key)
        self.assertEqual((r.key37, r.allele_id), ("3:12345:A:G", "15000"))
        self.assertIsNone(cr.from_variant_summary(vs_row(), "T1").key37)

    def test_missing_ids(self):
        r = cr.from_variant_summary(vs_row(**{"#AlleleID": "-1", "GeneID": "-1", "VariationID": ""}), "T1")
        self.assertIsNone(r.allele_id)
        self.assertIsNone(r.variation_id)
        self.assertEqual(r.gene, "missing")

    def test_vcf_line(self):
        line = ("3\t12345\t9000\tA\tG\t.\t.\tALLELEID=15000;CLNSIG=Likely_pathogenic;"
                "CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;GENEINFO=GENE1:101\n")
        (r,) = cr.from_vcf_line(line, "E-dev")
        self.assertEqual((r.key, r.allele_id, r.variation_id, r.stars, r.seen, r.gene),
                         ("3:12345:A:G", "15000", "9000", 2, True, "101"))
        (r,) = cr.from_vcf_line("3\t1\t5\tA\t.\t.\t.\tALLELEID=1;CLNSIG=Pathogenic", "E-dev")
        self.assertEqual((r.key, r.allele_id, r.variation_id, r.seen), (None, "1", "5", True))

    def test_identity_tuple(self):
        a = cr.from_variant_summary(vs_row(), "T1")
        b = cr.from_variant_summary(vs_row(ReviewStatus="reviewed by expert panel"), "T1")
        self.assertNotEqual(a.identity(), b.identity())


if __name__ == "__main__":
    unittest.main()
