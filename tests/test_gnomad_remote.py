"""The remote tabix reader on a synthetic BGZF file and index held in memory."""
import gzip
import struct
import sys
import unittest
import zlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import gnomad_remote as gr  # noqa: E402


def bgzf_block(payload: bytes) -> bytes:
    comp = zlib.compressobj(6, zlib.DEFLATED, -15)
    body = comp.compress(payload) + comp.flush()
    bsize = 18 + len(body) + 8 - 1
    header = b"\x1f\x8b\x08\x04" + b"\x00" * 4 + b"\x00\xff" + struct.pack("<H", 6) + b"BC" + struct.pack("<HH", 2, bsize)
    return header + body + struct.pack("<II", zlib.crc32(payload), len(payload))


def build(n_records: int, info_size: int = 900, block_payload: int = 4000):
    """A VCF of n records at positions 10, 20, ... on chr21 (two alleles at every 7th position),
    split into BGZF blocks at arbitrary byte boundaries, and a tabix index with one chunk in a 128 kb bin."""
    lines = []
    for i in range(1, n_records + 1):
        pos = 10 * i
        info = f"AC_joint={i};AN_joint=1000;AF_joint={i / 1000};fafmax_faf95_max_joint=0.{i:04d};" + "X" * info_size
        lines.append(f"chr21\t{pos}\t.\tA\tG\t.\tPASS\t{info}")
        if i % 7 == 0:
            lines.append(f"chr21\t{pos}\t.\tA\tT\t.\tAC0\t{info}")
    text = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "\n".join(lines) + "\n").encode()
    first_record = text.index(b"chr21")
    data, voff_first, coff = b"", None, 0
    for start in range(0, len(text), block_payload):
        payload = text[start:start + block_payload]
        if voff_first is None and start <= first_record < start + block_payload:
            voff_first = (coff << 16) | (first_record - start)
        block = bgzf_block(payload)
        data += block
        coff += len(block)
    data += bgzf_block(b"")                                          # EOF block
    voff_end = coff << 16
    names = b"chr21\x00"
    body = b"TBI\x01" + struct.pack("<iiiiiii", 1, 2, 1, 2, 0, ord("#"), 0) + struct.pack("<i", len(names)) + names
    body += struct.pack("<i", 1) + struct.pack("<Ii", 585, 1) + struct.pack("<QQ", voff_first, voff_end)
    body += struct.pack("<i", 2) + struct.pack("<QQ", voff_first, voff_first)
    return data, gzip.compress(body), text


class Remote(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data, tbi, cls.text = build(3000)
        cls.index = gr.parse_tbi(tbi)
        cls.calls = 0

        def fetch(url, a, b):
            cls.calls += 1
            return cls.data[a:b]
        cls.fetch = staticmethod(fetch)

    def truth(self, pos):
        return sorted((p[3], p[4]) for p in (l.split("\t") for l in self.text.decode().split("\n"))
                      if len(p) > 4 and p[0] == "chr21" and int(p[1]) == pos)

    def test_index(self):
        self.assertEqual(self.index.names, ["chr21"])
        self.assertEqual(len(gr.chunks_for(self.index, "chr21", 100)), 1)

    def test_search_equals_full_read(self):
        old = (gr.PROBE, gr.SCAN)
        gr.PROBE, gr.SCAN = 3000, 6000                                # force probing on a small file
        try:
            for pos in (10, 20, 70, 15, 700, 7000, 14000, 29990, 30000, 30010, 21000, 5):
                with self.subTest(pos=pos):
                    got = sorted((r[3], r[4]) for r in gr.search("u", self.index, "21", pos, self.fetch))
                    self.assertEqual(got, self.truth(pos))
        finally:
            gr.PROBE, gr.SCAN = old

    def test_lookup_fields(self):
        got = gr.lookup_fast("u", self.index, "21", 70, [("A", "G"), ("A", "T"), ("A", "C")], self.fetch)
        self.assertEqual(got[("A", "G")]["AC_joint"], "7")
        self.assertEqual(got[("A", "T")]["FILTER"], "AC0")
        self.assertIsNone(got[("A", "C")])
        self.assertEqual(gr.as_float(got[("A", "G")]["fafmax_faf95_max_joint"]), 0.0007)

    def test_reg2bins_level5(self):
        self.assertIn(4681 + (20000 >> 14), gr.reg2bins(20000, 20001))


if __name__ == "__main__":
    unittest.main()
