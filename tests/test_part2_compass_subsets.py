"""COMPASS source-subset analysis helpers on synthetic data (no real scores)."""
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import part2_analysis as p2  # noqa: E402
import part2_compass_subsets as cs  # noqa: E402


class CompassSubsets(unittest.TestCase):
    def test_subset_masks_split_the_analysed_compass_set_by_the_authors_labels(self):
        d = pd.DataFrame({"dataset": ["compass"] * 5 + ["opensplice"],
                          "subset": ["random_snv", "exac", "random_snv", "exac", "random_snv", "random_snv"]})
        analysed = pd.Series([True, True, True, False, True, True])      # row 3: excluded upstream
        masks = cs.subset_masks(d, analysed & (d["dataset"] == "compass"))
        self.assertEqual(list(masks), [cs.RANDOM, cs.EXAC])
        self.assertEqual(masks[cs.RANDOM].tolist(), [True, False, True, False, True, False])
        self.assertEqual(masks[cs.EXAC].tolist(), [False, True, False, False, False, False])

    def test_pooled_replaces_compass_with_its_random_subset_only(self):
        planned = [{"dataset": ds, "pair": "AVI - SpliceAI", "difference": diff, "se": 0.02}
                   for ds, diff in (("compass", 0.3), ("kitzman_wt1_ex9", -0.1), ("braun2018_mst1r_ex11", 0.0))]
        subsets = [{"dataset": cs.RANDOM, "pair": "AVI - SpliceAI", "difference": -0.05, "se": 0.03},
                   {"dataset": cs.EXAC, "pair": "AVI - SpliceAI", "difference": 0.9, "se": 0.01}]
        main_sets = ["compass", "kitzman_wt1_ex9", "braun2018_mst1r_ex11"]
        out = cs.pooled_replacing(planned, subsets, main_sets, ["AVI - SpliceAI"])["AVI - SpliceAI"]
        expected = p2.dersimonian_laird([-0.1, 0.0, -0.05], [0.02, 0.02, 0.03])
        self.assertAlmostEqual(out["estimate"], expected["estimate"])
        self.assertEqual(out["ci_hartung_knapp"], expected["ci_hartung_knapp"])
        self.assertEqual(out["datasets"], ["kitzman_wt1_ex9", "braun2018_mst1r_ex11", cs.RANDOM])

    def test_planned_rows_are_read_with_numeric_differences(self):
        planned = cs.planned_rows(pd.DataFrame({"dataset": ["compass"], "pair": ["AVI - SpliceAI"],
                                                "difference": ["-0.1"], "se": ["0.02"]}))
        self.assertEqual(planned, [{"dataset": "compass", "pair": "AVI - SpliceAI", "difference": -0.1, "se": 0.02}])


if __name__ == "__main__":
    unittest.main()
