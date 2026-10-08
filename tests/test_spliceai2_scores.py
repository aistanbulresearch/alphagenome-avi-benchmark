"""Combining SpliceAI2 rows into one score per allele (no real scores)."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import spliceai2_scores as s2  # noqa: E402

HEADER = ["chrom", "pos", "ref", "alt", "strand", "donor_gain_delta_score_0", "donor_gain_delta_score_1",
          "donor_loss_delta_score_0", "acceptor_gain_delta_score_0", "acceptor_loss_delta_score_0",
          "spliceai2_summary_score"]


def row(ref, alt, strand, dg, dl, ag, al, summary):
    return ["chr7", "100", ref, alt, strand, str(dg), "0.99", str(dl), str(ag), str(al), str(summary)]


class Combine(unittest.TestCase):
    def setUp(self):
        self.cols = s2.columns(HEADER)

    def test_columns(self):
        self.assertEqual(self.cols["summary"], 10)
        self.assertEqual(sorted(self.cols["delta0"]), [5, 7, 8, 9])

    def test_max_over_strands_and_agreement(self):
        rows = [row("A", "G", "+", 0.1, 0.2, 0.0, 0.0, 0.2), row("A", "G", "-", 0.5, 0.0, 0.0, 0.3, 0.5),
                row("A", "T", "+", 0.9, 0.0, 0.0, 0.0, 0.9)]
        got = s2.combine(rows, "A", "G", self.cols)
        self.assertEqual(got, {"spliceai2_summary": "0.5", "spliceai2_maxdelta": "0.5", "rows": 2, "differs": 0})

    def test_disagreement_counted_and_k1_ignored(self):
        rows = [row("A", "G", "+", 0.1, 0.0, 0.0, 0.0, 0.3)]          # summary differs from max delta_0
        got = s2.combine(rows, "A", "G", self.cols)
        self.assertEqual((got["spliceai2_summary"], got["spliceai2_maxdelta"], got["differs"]), ("0.3", "0.1", 1))

    def test_missing_allele(self):
        self.assertIsNone(s2.combine([row("A", "T", "+", 0, 0, 0, 0, 0)], "A", "G", self.cols))


class FetchRetry(unittest.TestCase):
    def test_incomplete_read_is_retried(self):
        import http.client
        from unittest import mock

        calls = {"n": 0}

        class Resp:
            def __enter__(self):
                return self

            def __exit__(self, *a):
                return False

            def read(self):
                calls["n"] += 1
                if calls["n"] == 1:
                    raise http.client.IncompleteRead(b"", 10)
                return b"x" * 10

        remote = s2.Remote("t")
        remote.signed["u"] = ("https://signed.example/f", 1e18)
        with mock.patch("urllib.request.urlopen", return_value=Resp()), mock.patch("time.sleep"):
            self.assertEqual(remote.fetch("u", 0, 10), b"x" * 10)
        self.assertEqual(calls["n"], 2)


if __name__ == "__main__":
    unittest.main()
