"""Part 1 tables on a synthetic variant table (no real scores)."""
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import part1_analysis as pa  # noqa: E402
import part1_data as pdata  # noqa: E402


def table(n=1200, seed=0):
    rng = np.random.default_rng(seed)
    y = (rng.random(n) < 0.2).astype(int)
    df = pd.DataFrame({
        "variation_id": [str(i) for i in range(n)], "gene_id": [str(i % 90) for i in range(n)], "y": y,
        "label": np.where(y == 1, "P/LP", "B/LB"), "stars_n": 1, "mane_select": "yes",
        "new_label": np.where(np.arange(n) < n // 2, "yes", "no"),
        "chromosome_group": rng.choice(["untouched", "validation", "training"], n),
        "region": rng.choice(["acceptor_region", "near_intronic", "deep_intronic", "utr3"], n),
        "rare": True, "in_gnomad": rng.random(n) < 0.7, "inside_long_unknown_record": "no",
        "mentions_any_splice_predictor": rng.choice(["yes", "no"], n),
        "calibration": np.arange(n) >= n // 2,
    })
    df["AVI"] = y * 1.5 + rng.normal(size=n)
    df["SpliceAI"] = y * 1.0 + rng.normal(size=n)
    df["CADD"] = np.where(rng.random(n) < 0.05, np.nan, y * 0.5 + rng.normal(size=n))
    return df


class Part1(unittest.TestCase):
    def test_sets_and_regions(self):
        v = table()
        s = pdata.sets(v)
        self.assertTrue((s["main"] <= (v["new_label"] == "yes")).all())
        self.assertFalse((s["main"] & (v["chromosome_group"] == "training")).any())
        r = pdata.regions(v)
        self.assertTrue((r["intronic >=3 bp (main)"] == v["region"].isin(pdata.MAIN_REGION)).all())

    def test_discrimination_and_pairs(self):
        v = table()
        rows, diffs = pa.discrimination(v, pdata.sets(v)["main"], ["AVI", "SpliceAI", "CADD"], 50)
        by = {r["tool"]: r for r in rows}
        self.assertGreater(by["AVI"]["auroc"], by["SpliceAI"]["auroc"])
        lo, hi = by["AVI"]["auroc_ci"]
        self.assertLess(lo, by["AVI"]["auroc"])
        self.assertGreater(hi, by["AVI"]["auroc"])
        pair = {d["pair"]: d for d in diffs}["AVI - SpliceAI"]
        self.assertGreater(pair["difference"], 0)

    def test_coverage(self):
        v = table()
        c = pa.coverage(v, {"all": pd.Series(True, index=v.index)}, ["AVI", "CADD"])
        self.assertEqual(int(c.loc[c.tool == "AVI", "scored"].iloc[0]), len(v))
        self.assertLess(int(c.loc[c.tool == "CADD", "scored"].iloc[0]), len(v))

    def test_calibration_tables_run(self):
        v = table(3000, 2)
        out = pa.calibration_tables(v, pdata.sets(v), ["AVI", "SpliceAI"], 20, 50)
        self.assertIn("SpliceAI", out)
        self.assertIn("published PP3 (>=0.2) pathogenic", out["SpliceAI"]["evaluated_on_main"])
        self.assertGreater(out["AVI"]["derivation"]["n_pos"], 0)


if __name__ == "__main__":
    unittest.main()
