"""Counts per dataset and position class (classes come from src/part2_exons.py)."""
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import part2_positions as pp  # noqa: E402


class Positions(unittest.TestCase):
    def test_table(self):
        d = pd.DataFrame({"dataset": ["x", "x", "x"], "y": [1, 0, 1]})
        t = pp.table(d, pd.Series(["exonic", "intronic >50", "exonic"]))
        row = t[(t["dataset"] == "x") & (t["class"] == "exonic")].iloc[0]
        self.assertEqual((row["variants"], row["splice_altering"]), (2, 2))


if __name__ == "__main__":
    unittest.main()
