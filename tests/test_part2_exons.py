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

    def test_classes(self):
        self.assertEqual([px.distance_class(x) for x in (0, 1, 2, 3, 50, 51, None)],
                         ["exonic", "intronic 1-2", "intronic 1-2", "intronic 3-50", "intronic 3-50", "intronic >50",
                          "unmapped"])


if __name__ == "__main__":
    unittest.main()
