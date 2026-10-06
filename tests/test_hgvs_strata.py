import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import hgvs_strata as hs  # noqa: E402

TX = "NM_000001.1(GENE1):"


def a(expr, prot=None):
    name = TX + expr + (f" ({prot})" if prot else "")
    return hs.assign(name)


class IntronicBoundaries(unittest.TestCase):
    CASES = [
        ("c.100+1G>A", "canonical", 1), ("c.100+2T>C", "canonical", 2),
        ("c.101-1G>A", "canonical", -1), ("c.101-2A>G", "canonical", -2),
        ("c.100+3A>G", "donor_region", 3), ("c.100+6T>A", "donor_region", 6),
        ("c.100+7G>C", "near_intronic", 7), ("c.100+50C>T", "near_intronic", 50),
        ("c.100+51C>T", "deep_intronic", 51),
        ("c.101-3C>G", "acceptor_region", -3), ("c.101-20T>C", "acceptor_region", -20),
        ("c.101-21T>C", "near_intronic", -21), ("c.101-50A>G", "near_intronic", -50),
        ("c.101-51A>G", "deep_intronic", -51), ("c.101-5000A>G", "deep_intronic", -5000),
    ]

    def test_every_boundary(self):
        for expr, stratum, offset in self.CASES:
            with self.subTest(expr=expr):
                r = a(expr)
                self.assertEqual((r.family, r.stratum, r.signed_offset), ("intronic", stratum, offset))
                self.assertFalse(r.utr_anchored)

    def test_utr_anchored_offsets(self):
        r = a("c.-45+4C>T")
        self.assertEqual((r.stratum, r.signed_offset, r.utr_anchored), ("donor_region", 4, True))
        r = a("c.*20-60G>A")
        self.assertEqual((r.stratum, r.signed_offset, r.utr_anchored), ("deep_intronic", -60, True))

    def test_zero_offset_rejected(self):
        with self.assertRaises(ValueError):
            hs.intronic_stratum("+", 0)
        self.assertEqual(a("c.100+0G>A").family, "unparsed")


class UtrAndCoding(unittest.TestCase):
    def test_utr(self):
        self.assertEqual((a("c.-15C>T").family, a("c.-15C>T").stratum), ("utr", "utr5"))
        self.assertEqual(a("c.*30A>G").stratum, "utr3")

    def test_protein_classes(self):
        cases = [
            ("p.Arg12Ter", "nonsense"), ("p.Arg12*", "nonsense"), ("p.(Arg12Ter)", "nonsense"),
            ("p.Arg12His", "missense"), ("p.Arg12=", "synonymous"), ("p.(=)", "synonymous"),
            ("p.Ter400=", "synonymous"), ("p.Gln1756fs", "frameshift"), ("p.Gln5ProfsTer3", "frameshift"),
            ("p.Ter110GlnextTer17", "extension"), ("p.Met1ext-5", "extension"),
            ("p.Lys5del", "coding_indel"), ("p.Lys5dup", "coding_indel"),
            ("p.Lys5_Leu6insGly", "coding_indel"), ("p.Lys5delinsAsnGlu", "coding_indel"),
            ("p.Met1?", "coding_unknown"), ("p.?", "coding_unknown"), ("p.0", "coding_other"),
            ("p.Arg12Arg", "coding_other"),
        ]
        for prot, stratum in cases:
            with self.subTest(prot=prot):
                self.assertEqual(hs.protein_class(prot), stratum)

    def test_coding_from_name(self):
        r = a("c.36G>A", "p.Lys12=")
        self.assertEqual((r.family, r.stratum, r.flags), ("coding", "synonymous", ("annotation_pending",)))
        r = a("c.36G>A")
        self.assertEqual(r.stratum, "coding_unknown")
        r = a("c.5266dup", "p.Gln1756fs")
        self.assertEqual(r.stratum, "frameshift")


class Ranges(unittest.TestCase):
    def test_same_intron_side(self):
        r = a("c.100+4_100+9del")
        self.assertEqual((r.stratum, r.signed_offset, r.is_range), ("donor_region", 4, True))
        r = a("c.101-25_101-22dup")
        self.assertEqual((r.stratum, r.signed_offset), ("near_intronic", -22))

    def test_crossing_into_exon_is_canonical(self):
        self.assertEqual(a("c.95_100+3del").stratum, "canonical")
        self.assertEqual(a("c.95_100+3del").signed_offset, 1)
        self.assertEqual(a("c.101-4_103del").stratum, "canonical")
        self.assertEqual(a("c.101-4_103del").signed_offset, -1)
        self.assertEqual(a("c.100+60_102del").stratum, "canonical")

    def test_insertion_flanks(self):
        self.assertEqual(a("c.100_100+1insA").stratum, "canonical")
        self.assertEqual(a("c.100+2_100+3insT").stratum, "canonical")
        self.assertEqual(a("c.100+7_100+8insT").stratum, "near_intronic")
        self.assertEqual(a("c.100_101insGGC", "p.Lys34_Leu35insGly").stratum, "coding_indel")

    def test_within_one_intron(self):
        r = a("c.100+60_101-70del")
        self.assertEqual((r.stratum, r.signed_offset), ("deep_intronic", 60))
        r = a("c.100+8_101-80del")
        self.assertEqual(r.stratum, "near_intronic")
        r = a("c.-10+60_-9-60del")
        self.assertEqual((r.stratum, r.utr_anchored), ("deep_intronic", True))
        r = a("c.-1+60_1-60del")
        self.assertEqual(r.stratum, "deep_intronic")

    def test_ties(self):
        self.assertEqual(a("c.100+4_101-5del").family, "annotation_ambiguous")
        self.assertEqual(a("c.-3_3del").family, "annotation_ambiguous")
        self.assertEqual(a("c.1500+5_*1-5del").family, "annotation_ambiguous")
        # Same stratum on both sides is not a tie.
        self.assertEqual(a("c.100+30_101-30del").stratum, "near_intronic")

    def test_far_offset_keeps_range_extent(self):
        self.assertEqual((a("c.100+1_100+2del").signed_offset, a("c.100+1_100+2del").far_offset), (1, 2))
        self.assertEqual((a("c.100+1_100+3del").signed_offset, a("c.100+1_100+3del").far_offset), (1, 3))
        self.assertEqual(a("c.101-25_101-22dup").far_offset, -25)
        self.assertEqual(a("c.95_100+3del").far_offset, 3)
        self.assertIsNone(a("c.100+60_101-70del").far_offset)
        self.assertEqual(a("c.100+7A>G").far_offset, 7)

    def test_range_ends_kept(self):
        self.assertEqual(a("c.100+7_101-9del").ends, (7, -9))
        self.assertEqual(a("c.100+9_101-9del").ends, (9, -9))
        self.assertEqual(a("c.95_100+3del").ends, (None, 3))
        self.assertIsNone(a("c.100+7A>G").ends)

    def test_exonic_range(self):
        r = a("c.10_12del", "p.Lys4del")
        self.assertEqual((r.family, r.stratum), ("coding", "coding_indel"))
        self.assertEqual(a("c.-20_-18del").stratum, "utr5")


class Unparsed(unittest.TestCase):
    def test_reasons(self):
        for name in [None, "", "   ", "NC_012920.1:m.3243A>G", "NR_000002.1(GENE2):n.200A>G",
                     "GRCh38/hg38 1p36.33(chr1:1-100)x1", TX + "c.[100A>G;200C>T]",
                     TX + "c.(?_-30)_(*20_?)del", TX + "c.(100_200)del", TX + "c.100A=",
                     TX + "c.100_200foo", TX + "c.0A>G", TX + "c.100_?del",
                     "NM_000001.1(GENE1):c.100A>G (p.Lys34Glu); NM_2:c.1A>G"]:
            with self.subTest(name=name):
                r = hs.assign(name)
                self.assertEqual(r.family, "unparsed")
                self.assertTrue(r.reason)

    def test_name_without_gene_symbol(self):
        self.assertEqual(hs.assign("NM_000001.1:c.100+1G>A").stratum, "canonical")


if __name__ == "__main__":
    unittest.main()
