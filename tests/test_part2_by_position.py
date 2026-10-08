"""Part 2 tables by position class on synthetic data (no real scores)."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import part2_by_position as bp  # noqa: E402


def synthetic(seed=0):
    rng = np.random.default_rng(seed)
    parts = []
    for ds, n in (("compass", 600), ("opensplice", 600), ("mfass", 300)):
        y = (rng.random(n) < 0.3).astype(int)
        cls = rng.choice(["exonic", "intronic 3-50", "intronic >50", "intronic 1-2"], n, p=[0.4, 0.3, 0.28, 0.02])
        parts.append(pd.DataFrame({"dataset": ds, "gene": [f"G{i % 30}" for i in range(n)], "y": y,
                                   "effect": np.abs(y + rng.normal(size=n)), "class": cls,
                                   "AVI": y + rng.normal(size=n), "SpliceAI": 1.5 * y + rng.normal(size=n)}))
    return pd.concat(parts, ignore_index=True)


class ByPosition(unittest.TestCase):
    def test_minimum_counts_and_pooling(self):
        d = synthetic()
        sets = {ds: d["dataset"] == ds for ds in ("compass", "opensplice", "mfass")}
        rows, diffs, pooled = bp.by_position(d, sets, ["AVI", "SpliceAI"], 20, main_sets=["compass", "opensplice"])
        classes = {r["class"] for r in rows}
        self.assertNotIn("intronic 1-2", classes)                      # under 10 positives: not analysed
        self.assertTrue({"exonic", "intronic >50"} <= classes)
        p = pooled["intronic >50"]["AVI - SpliceAI"]
        self.assertEqual(p["datasets"], ["compass", "opensplice"])        # mfass is not a main dataset
        self.assertLess(p["estimate"], 0)


if __name__ == "__main__":
    unittest.main()
