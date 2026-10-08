"""Part 2 analysis helpers on synthetic data (no real scores)."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import part2_analysis as p2  # noqa: E402


class Part2(unittest.TestCase):
    def test_dersimonian_laird(self):
        same = p2.dersimonian_laird([0.1, 0.1, 0.1], [0.05, 0.05, 0.05])
        self.assertAlmostEqual(same["estimate"], 0.1)
        self.assertEqual(same["tau2"], 0.0)
        spread = p2.dersimonian_laird([-0.2, 0.0, 0.3], [0.02, 0.02, 0.02])
        self.assertGreater(spread["tau2"], 0)
        self.assertLess(spread["ci"][0], spread["estimate"])

    def test_heterogeneity_and_hartung_knapp(self):
        est, se = [-0.2, 0.0, 0.3, -0.05], [0.02, 0.03, 0.02, 0.04]
        r = p2.dersimonian_laird(est, se)
        self.assertAlmostEqual(r["i2"], max(0.0, (r["q"] - 3) / r["q"]))
        # Hartung-Knapp: same centre, t-based interval with the weighted residual variance
        w = 1 / (np.asarray(se) ** 2 + r["tau2"])
        q_hk = (w * (np.asarray(est) - r["estimate"]) ** 2).sum() / 3
        half = 3.182446305284263 * (q_hk / w.sum()) ** 0.5            # t(0.975, 3)
        self.assertAlmostEqual(r["ci_hartung_knapp"][0], r["estimate"] - half)
        self.assertAlmostEqual(r["ci_hartung_knapp"][1], r["estimate"] + half)
        same = p2.dersimonian_laird([0.1, 0.1, 0.1], [0.05, 0.05, 0.05])
        self.assertEqual(same["i2"], 0.0)

    def test_pooled_skips_comparison_and_developer_datasets(self):
        diffs = []
        for ds, diff in (("compass", -0.1), ("opensplice", 0.5), ("mfass", 0.9), ("kitzman_wt1_ex9", -0.1)):
            for pair in ("AVI - SpliceAI", "AVI - SpliceAI2"):
                diffs.append({"dataset": ds, "pair": pair, "difference": diff, "se": 0.02})
        main_sets = ["compass", "opensplice", "kitzman_wt1_ex9"]
        out = p2.pooled_differences(diffs, main_sets, ["AVI", "SpliceAI", "SpliceAI2"])
        self.assertEqual(out["AVI - SpliceAI"]["datasets"], ["compass", "opensplice", "kitzman_wt1_ex9"])
        self.assertEqual(out["AVI - SpliceAI2"]["datasets"], ["compass", "kitzman_wt1_ex9"])
        self.assertAlmostEqual(out["AVI - SpliceAI2"]["estimate"], -0.1)

    def test_development_hits(self):
        alleles = pd.DataFrame({"allele": ["1:100:A:G", "1:205:C:T", "2:50:G:A"], "chrom": ["1", "1", "2"],
                                "pos": ["100", "205", "50"]})
        with tempfile.TemporaryDirectory() as t:
            path = Path(t) / "pos.tsv"
            path.write_text("chrom\tstart\tend\tsource\n1\t90\t110\tX\n1\t300\t400\tY\n", encoding="utf-8")
            self.assertEqual(p2.development_hits(alleles, path), {"1:100:A:G"})

    def test_dataset_tables_single_gene_and_pairs(self):
        rng = np.random.default_rng(0)
        n = 400
        y = (rng.random(n) < 0.3).astype(int)
        d = pd.DataFrame({"dataset": "kitzman_wt1_ex9", "gene": "WT1", "y": y,
                          "effect": np.abs(y * 2 + rng.normal(size=n)),
                          "AVI": y + rng.normal(size=n), "SpliceAI": 0.5 * y + rng.normal(size=n)})
        rows, diffs = p2.dataset_tables(d, pd.Series(True, index=d.index), "kitzman_wt1_ex9", ["AVI", "SpliceAI"], 50)
        by = {r["tool"]: r for r in rows}
        self.assertGreater(by["AVI"]["auroc"], by["SpliceAI"]["auroc"])
        self.assertGreater(by["AVI"]["spearman_abs_effect"], 0)
        self.assertEqual(diffs[0]["pair"], "AVI - SpliceAI")
        self.assertGreater(diffs[0]["se"], 0)


if __name__ == "__main__":
    unittest.main()
