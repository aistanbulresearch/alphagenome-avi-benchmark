"""Paired AUPRC and sensitivity-at-specificity tables on the synthetic Part 1 table (no real scores)."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import part1_supplement as ps  # noqa: E402
from test_part1_analysis import table  # noqa: E402


class Supplement(unittest.TestCase):
    def test_paired_rows(self):
        v = table(3000, 4)
        mask = pd.Series(True, index=v.index)
        rows = ps.paired(v, mask, "AVI", "SpliceAI", 50)
        by = {r["metric"]: r for r in rows}
        self.assertEqual(set(by), {"AUPRC", "sensitivity at 95% specificity", "sensitivity at 99% specificity"})
        for r in rows:
            self.assertAlmostEqual(r["difference"], r["value_a"] - r["value_b"])
            lo, hi = r["ci"]
            self.assertLessEqual(lo, hi)
        self.assertGreater(by["AUPRC"]["difference"], 0)          # AVI is built with the stronger signal
        self.assertGreater(by["AUPRC"]["ci"][0], 0)

    def test_weights_enter_paired_metrics(self):
        v = table(3000, 4)
        v["Pangolin"] = v["SpliceAI"]
        v["w_Pangolin"] = np.where(v["y"] == 1, 1.0, 3.0)
        mask = pd.Series(True, index=v.index)
        weighted = {r["metric"]: r for r in ps.paired(v, mask, "AVI", "Pangolin", 5)}
        plain = {r["metric"]: r for r in ps.paired(v, mask, "AVI", "SpliceAI", 5)}
        self.assertNotAlmostEqual(weighted["AUPRC"]["value_a"], plain["AUPRC"]["value_a"])
        self.assertAlmostEqual(weighted["sensitivity at 99% specificity"]["value_b"],
                               plain["sensitivity at 99% specificity"]["value_b"], delta=0.02)

    def test_tables(self):
        v = table(3000, 4)
        out = ps.tables(v, ["AVI", "SpliceAI", "CADD"], 10)
        self.assertIn("main", set(out["set"]))
        self.assertIn("AVI - CADD", set(out["pair"]))


if __name__ == "__main__":
    unittest.main()
