"""bigWig data-block parsing on synthetic blocks (the live reader was checked against the UCSC API)."""
import struct
import sys
import unittest
import zlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import bigwig_remote as bw  # noqa: E402


def reader(blob: bytes) -> bw.BigWig:
    r = bw.BigWig.__new__(bw.BigWig)
    r.url, r.compressed, r.fetch = "u", True, (lambda url, a, b: blob[a:b])
    r.chroms, r._node_cache = {}, {}
    return r


class Blocks(unittest.TestCase):
    def test_fixed_step(self):
        head = struct.pack("<IIIIIBBH", 0, 100, 103, 1, 1, 3, 0, 3)
        blob = zlib.compress(head + struct.pack("<fff", 0.5, -1.0, 2.0))
        self.assertEqual(reader(blob).block_values(0, len(blob)), [(100, 101, 0.5), (101, 102, -1.0), (102, 103, 2.0)])

    def test_bedgraph_and_variable_step(self):
        a = struct.pack("<IIIIIBBH", 0, 10, 30, 0, 0, 1, 0, 2) + struct.pack("<IIf", 10, 20, 1.5) + struct.pack("<IIf", 20, 30, 2.5)
        b = struct.pack("<IIIIIBBH", 0, 40, 60, 0, 2, 2, 0, 2) + struct.pack("<If", 40, 3.0) + struct.pack("<If", 50, 4.0)
        blob = zlib.compress(a + b)
        self.assertEqual(reader(blob).block_values(0, len(blob)),
                         [(10, 20, 1.5), (20, 30, 2.5), (40, 42, 3.0), (50, 52, 4.0)])


if __name__ == "__main__":
    unittest.main()
