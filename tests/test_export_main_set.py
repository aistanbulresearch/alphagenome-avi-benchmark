"""The exported main-set list holds identities and labels only, never a tool score (synthetic rows)."""
import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import export_main_set as ems  # noqa: E402


class ExportMainSet(unittest.TestCase):
    def test_only_identity_columns_of_main_rows(self):
        v = pd.DataFrame({"variation_id": ["1", "2", "3"], "chrom": ["1", "2", "3"], "pos": ["10", "20", "30"],
                          "ref": ["A", "C", "G"], "alt": ["G", "T", "A"], "label": ["P/LP", "B/LB", "B/LB"],
                          "region": ["deep_intronic", "utr3", "synonymous"],
                          "chromosome_group": ["untouched", "validation", "untouched"],
                          "AVI": [0.9, 0.1, 0.2], "SpliceAI": [0.5, 0.0, 0.0]})
        out = ems.main_rows(v, pd.Series([True, False, True]))
        self.assertEqual(list(out.columns), list(ems.COLUMNS))
        self.assertEqual(list(out["variation_id"]), ["1", "3"])
        self.assertFalse(set(out.columns) & {"AVI", "SpliceAI"})


if __name__ == "__main__":
    unittest.main()
