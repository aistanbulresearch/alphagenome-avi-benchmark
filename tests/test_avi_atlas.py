"""AVI Atlas helpers with a fake client (no network, no real scores)."""
import sys
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import avi_atlas as aa  # noqa: E402


class FakeAnnData:
    def __init__(self, ids, raw, quant):
        import pandas as pd
        self.obs = pd.DataFrame({"variant": ids})
        self.X = np.array(raw, dtype=float).reshape(-1, 1)
        self.layers = {"quantiles": np.array(quant, dtype=float).reshape(-1, 1)}


def fake_query(rejected=()):
    def query(ids):
        if any(i in rejected for i in ids):
            raise ValueError("reference mismatch")
        return {aa.SCORER: FakeAnnData(ids, [0.1 * (k + 1) for k in range(len(ids))], [0.5] * len(ids))}
    return query


V = [{"variation_id": str(k), "chrom": "3", "pos": str(100 + k), "ref": "A", "alt": "G"} for k in range(4)]


class Atlas(unittest.TestCase):
    def test_variant_id(self):
        self.assertEqual(aa.variant_id(V[0]), "chr3:100:A>G")

    def test_chunk_scores_in_order(self):
        out = aa.run_chunk(fake_query(), V)
        self.assertEqual([round(o[1], 2) for o in out], [0.1, 0.2, 0.3, 0.4])
        self.assertEqual({o[3] for o in out}, {"ok"})

    def test_rejected_variant_falls_back_one_by_one(self):
        with mock.patch.object(aa.time, "sleep", lambda s: None):
            out = aa.run_chunk(fake_query(rejected={"chr3:101:A>G"}), V)
        status = {o[0]["variation_id"]: o[3] for o in out}
        self.assertEqual(status, {"0": "ok", "1": "error: ValueError", "2": "ok", "3": "ok"})


if __name__ == "__main__":
    unittest.main()
