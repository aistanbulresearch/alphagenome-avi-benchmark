"""Distance of Part 2 variants to the assayed exon, per dataset (synthetic rows)."""
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import part2_exons as px  # noqa: E402


class Exons(unittest.TestCase):
    def test_hgvs_offset(self):
        self.assertEqual(px.hgvs_distance("Intronic; c.292+131; NM_001122757.2"), 131)
        self.assertEqual(px.hgvs_distance("intronic; c.1447+40T>G"), 40)
        self.assertEqual(px.hgvs_distance("intronic; c.1306-12A>G"), 12)
        self.assertEqual(px.hgvs_distance("Missense; c.250A>G; NM_x"), 0)
        self.assertIsNone(px.hgvs_distance("no name here"))

    def test_interval_distance(self):
        self.assertEqual(px.interval_distance(150, [(100, 200)]), 0)
        self.assertEqual(px.interval_distance(95, [(100, 200), (300, 400)]), 5)
        self.assertEqual(px.interval_distance(260, [(100, 200), (300, 400)]), 40)

    def test_mfass_exon_from_construct(self):
        # construct [48981364, 48981534) 0-based, left flank 51, exon 57: exon 48981416..48981472 (1-based)
        row = {"start_hg38_0based": 48981364, "intron1_len": 51, "exon_len": 57}
        self.assertEqual(px.mfass_exon(row), (48981416, 48981472))
        self.assertEqual(px.interval_distance(48981506, [px.mfass_exon(row)]), 34)

    def test_exon_intervals_from_labelled_positions(self):
        d = pd.DataFrame({"key": ["e1"] * 4 + ["e2"] * 2, "pos": [10, 11, 12, 40, 100, 105],
                          "is_exon": [True, True, True, False, True, True]})
        self.assertEqual(px.intervals_from_labels(d, "key", "pos", "is_exon"), {"e1": [(10, 12)], "e2": [(100, 105)]})

    def test_opensplice_span_includes_the_exonic_bases_of_the_splice_site_regions(self):
        """The 5'SS region holds the last 3 exonic bases (-3..+6), the 3'SS region the first exonic base
        (-4..+1); which genomic side each lies on depends on the strand."""
        def rows(exon_id, regions, length):
            out = []
            for region, (a, b) in regions.items():
                out += [{"exon_id": exon_id, "region": region, "POS": p, "exon_length": length} for p in range(a, b + 1)]
            return out
        plus = rows("p", {"Intron up": (60, 95), "3'SS": (96, 100), "Exon": (101, 196), "5'SS": (197, 205),
                          "Intron down": (206, 224)}, 100)            # true exon 100..199
        minus = rows("m", {"Intron down": (61, 79), "5'SS": (80, 88), "Exon": (89, 184), "3'SS": (185, 189),
                           "Intron up": (190, 255)}, 100)             # true exon 86..185
        wrong = rows("w", {"3'SS": (96, 100), "Exon": (101, 196), "5'SS": (197, 205)}, 120)
        no_length = rows("n", {"3'SS": (96, 100), "Exon": (101, 250), "5'SS": (251, 259)}, None)
        no_exon = rows("x", {"3'SS": (96, 100), "5'SS": (101, 109)}, 10)
        spans = px.opensplice_spans(pd.DataFrame(plus + minus + wrong + no_length + no_exon))
        self.assertEqual(spans, {"p": [(100, 199)], "m": [(86, 185)]})

    def test_multi_exon_construct_other_exons_set_apart_intronic_distance_to_nearest_exon(self):
        """Inside the assayed exon: exonic; inside another exon of the construct: "other exon"; in an intron:
        bases to the nearest exon of the construct (Part 1's definition of the intronic regions)."""
        exons = [(100, 150), (300, 350), (500, 560)]                 # e.g. MST1R exons 10, 11, 12; 11 assayed
        self.assertEqual(px.construct_position(320, exons, 1), (0, None))
        self.assertEqual(px.construct_position(298, exons, 1), (2, None))
        self.assertEqual(px.construct_position(152, exons, 1), (2, None))     # 2 bases from exon 10
        self.assertEqual(px.construct_position(225, exons, 1), (75, None))
        self.assertEqual(px.construct_position(120, exons, 1), (None, "other exon"))
        self.assertEqual(px.construct_position(560, exons, 1), (None, "other exon"))
        self.assertEqual(px.construct_position(120, exons, 0), (0, None))

    def test_outer_sides_of_first_and_last_exons(self):
        gtf = [
            "chr1\tX\texon\t100\t200\t.\t+\t.\ttranscript_id \"A\"; tag \"MANE_Select\";",
            "chr1\tX\texon\t300\t400\t.\t+\t.\ttranscript_id \"A\"; tag \"MANE_Select\";",
            "chr1\tX\texon\t500\t600\t.\t+\t.\ttranscript_id \"A\"; tag \"MANE_Select\";",
            "chr2\tX\texon\t100\t200\t.\t-\t.\ttranscript_id \"B\"; tag \"MANE_Select\";",
            "chr2\tX\texon\t300\t400\t.\t-\t.\ttranscript_id \"B\"; tag \"MANE_Select\";",
            "chr3\tX\texon\t100\t200\t.\t+\t.\ttranscript_id \"C\";",
        ]
        low, high = px.outer_sides(gtf)
        # plus strand: first exon (100-200) has its low side outside the gene, last (500-600) its high side
        self.assertEqual(low["chr1"], {200})
        self.assertEqual(high["chr1"], {500})
        # minus strand: first exon is 300-400 (high side outside), last is 100-200 (low side outside)
        self.assertEqual(low["chr2"], {200})
        self.assertEqual(high["chr2"], {300})
        self.assertNotIn("chr3", low)
        self.assertTrue(px.is_flank(90, (100, 200), low["chr1"], high["chr1"]))
        self.assertFalse(px.is_flank(210, (100, 200), low["chr1"], high["chr1"]))
        self.assertFalse(px.is_flank(150, (100, 200), low["chr1"], high["chr1"]))
        self.assertTrue(px.is_flank(610, (500, 600), low["chr1"], high["chr1"]))
        self.assertFalse(px.is_flank(290, (300, 400), low["chr1"], high["chr1"]))

    def test_class_overrides_distance(self):
        self.assertEqual(px.position_class(5, None), "intronic 3-50")
        self.assertEqual(px.position_class(5, "flank"), "flank")
        self.assertEqual(px.position_class(None, "other exon"), "other exon")
        self.assertEqual(px.position_class(float("nan"), float("nan")), "unmapped")    # missing values from a frame
        self.assertEqual(px.position_class(4.0, float("nan")), "intronic 3-50")

    def test_compass_span_must_be_one_annotated_exon(self):
        exons = {(100, 150), (180, 200)}
        self.assertTrue(px.annotated_exon((100, 150), exons))
        self.assertFalse(px.annotated_exon((100, 200), exons))       # ends of two different exons
        self.assertFalse(px.annotated_exon(None, exons))

    def test_classes(self):
        self.assertEqual([px.distance_class(x) for x in (0, 1, 2, 3, 50, 51, None)],
                         ["exonic", "intronic 1-2", "intronic 1-2", "intronic 3-50", "intronic 3-50", "intronic >50",
                          "unmapped"])


if __name__ == "__main__":
    unittest.main()
