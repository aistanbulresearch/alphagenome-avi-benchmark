import os
import sys
import unittest
from pathlib import Path

os.environ["AVI_EVAL_SYNTHETIC"] = "1"
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import liftover  # noqa: E402

# Two chains: chr1 forward (two blocks with a gap) and chr2 mapped to the reverse strand of chr5.
CHAIN = """chain 1000 chr1 1000 + 100 300 chr1 1200 + 200 400 1
50 30 30
120

chain 500 chr2 500 + 0 100 chr5 1000 - 0 100 2
100
"""


class Chain(unittest.TestCase):
    def setUp(self):
        self.lift = liftover.parse_chain(CHAIN.splitlines())

    def test_forward_blocks(self):
        # 1-based base 101 is 0-based 100 -> q 200 -> 1-based 201.
        self.assertEqual(self.lift(("1", 101, 101)), [("1", 201, 201)])
        self.assertEqual(self.lift(("chr1", 150, 150)), [("1", 250, 250)])
        # Base in the gap (0-based 150..179) does not convert.
        self.assertEqual(self.lift(("1", 160, 160)), [])
        # Second block starts at t 180 -> q 280.
        self.assertEqual(self.lift(("1", 181, 181)), [("1", 281, 281)])

    def test_span_split_by_gap(self):
        self.assertEqual(self.lift(("1", 149, 182)), [("1", 249, 250), ("1", 281, 282)])

    def test_reverse_strand(self):
        # t 0-based 0 -> q 0 on the minus strand -> plus strand 999 (1-based 1000).
        self.assertEqual(self.lift(("2", 1, 1)), [("5", 1000, 1000)])
        self.assertEqual(self.lift(("2", 1, 10)), [("5", 991, 1000)])

    def test_unknown_chromosome(self):
        self.assertEqual(self.lift(("7", 1, 1)), [])


if __name__ == "__main__":
    unittest.main()
