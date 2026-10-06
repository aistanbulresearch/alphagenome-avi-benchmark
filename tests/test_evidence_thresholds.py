import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import evidence_thresholds as et  # noqa: E402


class PosteriorTests(unittest.TestCase):
    def test_reproduces_printed_pejaver_posteriors(self):
        for level, printed in et.PEJAVER_2022_POSTERIOR_PRINTED.items():
            value = et.posterior_from_lr(et.PRIMARY.lr_pathogenic[level], et.PRIMARY.prior)
            self.assertAlmostEqual(round(value, 4), printed, places=4, msg=level)

    def test_moderate_and_strong_posteriors_fixed(self):
        table = et.posterior_table(et.PRIMARY)[et.PRIMARY.prior]
        self.assertAlmostEqual(round(table["moderate"], 4), 0.2108, places=4)
        self.assertAlmostEqual(round(table["strong"], 4), 0.6074, places=4)

    def test_pejaver_values_follow_eighth_root_scaling(self):
        # 1124^(1/8) = 2.41: the printed supporting value is consistent with the scaling (four digits).
        scaled = et.scaled_thresholds(et.PRIMARY.lr_pathogenic["very_strong"])
        self.assertAlmostEqual(scaled["supporting"], et.PRIMARY.lr_pathogenic["supporting"], delta=0.005)
        self.assertAlmostEqual(scaled["moderate"], et.PRIMARY.lr_pathogenic["moderate"], delta=0.01)
        self.assertAlmostEqual(scaled["strong"], et.PRIMARY.lr_pathogenic["strong"], delta=0.05)

    def test_tavtigian_set(self):
        self.assertAlmostEqual(et.TAVTIGIAN_2018.lr_pathogenic["supporting"], 2.08, places=2)
        self.assertAlmostEqual(et.TAVTIGIAN_2018.lr_pathogenic["moderate"], 4.33, places=2)
        self.assertAlmostEqual(et.TAVTIGIAN_2018.lr_pathogenic["strong"], 18.7, places=1)

    def test_benign_reciprocals(self):
        benign = et.PRIMARY.lr_benign()
        self.assertAlmostEqual(benign["supporting"] * et.PRIMARY.lr_pathogenic["supporting"], 1.0)

    def test_strength_of(self):
        self.assertIsNone(et.strength_of(2.40))
        self.assertEqual(et.strength_of(2.406), "supporting")
        self.assertEqual(et.strength_of(6.0), "moderate")
        self.assertEqual(et.strength_of(2000.0), "very_strong")

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            et.posterior_from_lr(0, 0.1)
        with self.assertRaises(ValueError):
            et.posterior_from_lr(2, 1.0)
        with self.assertRaises(ValueError):
            et.scaled_thresholds(1.0)


if __name__ == "__main__":
    unittest.main()
