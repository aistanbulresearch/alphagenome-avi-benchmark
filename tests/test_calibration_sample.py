"""The calibration sample: all P/LP, a proportional stratified B/LB sample, reproducible with its seed."""
import sys
import unittest
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import calibration_sample as cs  # noqa: E402


def row(i, label, region, group="untouched"):
    return {"variation_id": str(i), "label": label, "region": region, "chromosome_group": group,
            "chrom": "3", "stars": "1", "family": "intronic", "new_label": "no"}


ROWS = ([row(i, "P/LP", "deep_intronic") for i in range(30)]
        + [row(100 + i, "B/LB", "near_intronic") for i in range(600)]
        + [row(1000 + i, "B/LB", "deep_intronic", "validation") for i in range(200)])


class Sample(unittest.TestCase):
    def test_all_pathogenic_and_proportional_benign(self):
        s = cs.draw(ROWS, 80, 7)
        c = Counter((r["label"], r["region"]) for r in s)
        self.assertEqual(c[("P/LP", "deep_intronic")], 30)
        self.assertEqual(c[("B/LB", "near_intronic")], 60)
        self.assertEqual(c[("B/LB", "deep_intronic")], 20)

    def test_same_seed_same_sample(self):
        ids = lambda s: [r["variation_id"] for r in s]  # noqa: E731
        self.assertEqual(ids(cs.draw(ROWS, 80, 7)), ids(cs.draw(list(reversed(ROWS)), 80, 7)))
        self.assertNotEqual(ids(cs.draw(ROWS, 80, 7)), ids(cs.draw(ROWS, 80, 8)))

    def test_eligibility(self):
        self.assertTrue(cs.eligible(row(1, "B/LB", "near_intronic")))
        self.assertFalse(cs.eligible(row(1, "P/LP", "canonical")))
        self.assertFalse(cs.eligible(dict(row(1, "B/LB", "near_intronic"), new_label="yes")))
        self.assertFalse(cs.eligible(dict(row(1, "B/LB", "near_intronic"), stars="0")))


if __name__ == "__main__":
    unittest.main()
