"""Weighted AUROC/AUPRC against direct definitions; the cluster bootstrap weights."""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import metrics as m  # noqa: E402


def auroc_pairs(s, y, w):
    num = den = 0.0
    for i in range(len(s)):
        for j in range(len(s)):
            if y[i] == 1 and y[j] == 0:
                ww = w[i] * w[j]
                den += ww
                num += ww * (1.0 if s[i] > s[j] else 0.5 if s[i] == s[j] else 0.0)
    return num / den


def ap_direct(s, y, w):
    total = sum(wi for wi, yi in zip(w, y) if yi == 1)
    out = 0.0
    for t in sorted(set(s), reverse=True):
        sel = [i for i in range(len(s)) if s[i] >= t]
        tp = sum(w[i] for i in sel if y[i] == 1)
        prec = tp / sum(w[i] for i in sel)
        new_pos = sum(w[i] for i in range(len(s)) if s[i] == t and y[i] == 1)
        out += prec * new_pos
    return out / total


class Metrics(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(1)
        self.s = np.round(rng.normal(size=60), 1)                  # rounding creates ties
        self.y = (rng.random(60) < 0.3).astype(int)
        self.w = rng.integers(0, 4, 60)

    def test_auroc(self):
        self.assertAlmostEqual(m.auroc(self.s, self.y), auroc_pairs(self.s, self.y, np.ones(60)))
        self.assertAlmostEqual(m.auroc(self.s, self.y, self.w), auroc_pairs(self.s, self.y, self.w))
        self.assertEqual(m.auroc(np.array([1.0, 2.0]), np.array([1, 0])), 0.0)
        self.assertTrue(np.isnan(m.auroc(np.array([1.0, 2.0]), np.array([1, 1]))))

    def test_auprc(self):
        self.assertAlmostEqual(m.auprc(self.s, self.y), ap_direct(self.s, self.y, np.ones(60)))
        self.assertAlmostEqual(m.auprc(self.s, self.y, self.w), ap_direct(self.s, self.y, self.w))

    def test_sensitivity_at_specificity(self):
        def direct(s, y, w, spec):
            neg, pos = sum(w[y == 0]), sum(w[y == 1])
            best = 0.0
            for t in sorted(set(s)):
                if sum(w[(s >= t) & (y == 0)]) <= (1 - spec) * neg:
                    best = sum(w[(s >= t) & (y == 1)]) / pos
                    break
            return best
        for spec in (0.5, 0.8, 0.95, 0.99):
            with self.subTest(spec=spec):
                self.assertAlmostEqual(m.sensitivity_at_specificity(self.s, self.y, spec),
                                       direct(self.s, self.y, np.ones(60), spec))
                self.assertAlmostEqual(m.sensitivity_at_specificity(self.s, self.y, spec, self.w),
                                       direct(self.s, self.y, self.w, spec))
        s, y = np.array([0.1, 0.2, 0.9, 0.8]), np.array([0, 0, 1, 1])
        self.assertEqual(m.sensitivity_at_specificity(s, y, 0.99), 1.0)
        self.assertEqual(m.sensitivity_at_specificity(np.array([1.0, 0.5]), np.array([0, 1]), 0.99), 0.0)

    def test_cluster_weights(self):
        clusters = np.array(["a", "a", "b", "c", "c", "c"])
        for w in m.cluster_weights(clusters, replicates=20, seed=3):
            self.assertEqual(w[0], w[1])
            self.assertEqual(len(set(w[3:])), 1)
            self.assertEqual(sum({c: wi for c, wi in zip(clusters, w)}.values()), 3)

    def test_interval_with_infinite_values(self):
        lo, hi = m.interval([1.0, 2.0, 3.0] * 40 + [float("inf")] * 10)
        self.assertFalse(np.isnan(lo) or np.isnan(hi))
        self.assertEqual(hi, float("inf"))
        lo, hi = m.interval(list(np.linspace(0, 1, 101)))
        self.assertAlmostEqual(lo, 0.025)                    # finite values: unchanged interpolation

    def test_clusters_of(self):
        c = m.clusters_of(["101", "", "-1", "101", "5;6"])
        self.assertEqual(c[0], c[3])
        self.assertEqual(len(set(c)), 4)


if __name__ == "__main__":
    unittest.main()
