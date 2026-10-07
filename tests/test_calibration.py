"""Calibration on simulated scores with a known likelihood ratio."""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import calibration as c  # noqa: E402


class Calibration(unittest.TestCase):
    def test_local_lr_known(self):
        # benign uniform on [0, 1); pathogenic: half uniform on [0, 1), half uniform on [0.8, 1)
        rng = np.random.default_rng(0)
        neg = rng.random(20000)
        pos = np.concatenate([rng.random(10000), 0.8 + 0.2 * rng.random(10000)])
        s = np.concatenate([pos, neg])
        y = np.concatenate([np.ones(len(pos), int), np.zeros(len(neg), int)])
        grid = np.array([0.4, 0.9])
        lr = c.local_lr(s, y, grid, np.array([0.05, 0.05]))
        self.assertAlmostEqual(lr[0], 0.5, delta=0.1)           # density ratio 0.5 below 0.8
        self.assertAlmostEqual(lr[1], 3.0, delta=0.3)           # 0.5 + 0.5 / 0.2 = 3 above 0.8

    def test_thresholds_rules(self):
        grid = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
        lower = np.array([0.1, 0.3, 2.5, 2.0, 6.0])
        upper = np.array([0.3, 0.4, 3.0, 4.0, 8.0])
        t = c.thresholds(grid, lower, upper)
        self.assertEqual(t["supporting"][0], 0.5)               # 2.0 at 0.4 breaks "every grid score >= t"
        self.assertEqual(t["moderate"][0], 0.5)
        self.assertIsNone(t["strong"][0])
        self.assertEqual(t["supporting"][1], 0.2)               # upper <= 1/2.406 = 0.416 at 0.1 and 0.2
        self.assertEqual(c.thresholds(grid, lower, np.array([0.3, 0.5, 0.1, 0.1, 0.1]))["supporting"][1], 0.1)

    def test_bounds_without_nan_when_windows_lack_benign(self):
        # pathogenic-only top of the score range: lr+ is infinite in many replicates there
        rng = np.random.default_rng(3)
        s = np.concatenate([rng.normal(0, 1, 2000), 6 + rng.random(150)])
        y = np.concatenate([np.zeros(2000, int), np.ones(150, int)])
        grid = c.grid_points(s)
        lower, upper = c.bootstrap_bounds(s, y, grid, c.windows(s, grid), replicates=200)
        self.assertFalse(np.isnan(lower).any())
        self.assertFalse(np.isnan(upper).any())
        self.assertTrue(np.isinf(lower[-1]))

    def test_level_intervals_and_interval_lr(self):
        thr = {"supporting": (0.2, 0.1), "moderate": (0.4, 0.02), "strong": (0.8, None), "very_strong": (None, None)}
        iv = c.level_intervals(thr)
        self.assertIn(("supporting", "ge", 0.2, 0.4), iv)
        self.assertIn(("strong", "ge", 0.8, float("inf")), iv)
        self.assertIn(("supporting", "le", 0.02, 0.1), iv)
        self.assertIn(("moderate", "le", float("-inf"), 0.02), iv)
        s = np.array([0.1, 0.3, 0.3, 0.9, 0.9, 0.05, 0.3])
        y = np.array([1, 1, 1, 1, 1, 0, 0])
        # supporting interval [0.2, 0.4): P 2/5, B 1/2 -> 0.8; cumulative >= 0.2: P 4/5, B 1/2 -> 1.6
        self.assertAlmostEqual(c.interval_lr(s, y, 0.2, 0.4, "ge"), 0.8)
        self.assertAlmostEqual(c.likelihood_ratio(s, y, 0.2, "ge"), 1.6)

    def test_calibrate_runs_and_lr_on_test(self):
        rng = np.random.default_rng(1)
        s = np.concatenate([rng.normal(2, 1, 600), rng.normal(0, 1, 3000)])
        y = np.concatenate([np.ones(600, int), np.zeros(3000, int)])
        out = c.calibrate(s, y, replicates=200)
        p, b = out["thresholds"]["supporting"]
        self.assertIsNotNone(p)
        self.assertIsNotNone(b)
        self.assertLess(b, p)
        self.assertGreaterEqual(c.likelihood_ratio(s, y, p, "ge"), 2.406 * 0.8)
        self.assertLessEqual(c.likelihood_ratio(s, y, b, "le"), 1 / 2.406 * 1.25)


if __name__ == "__main__":
    unittest.main()
