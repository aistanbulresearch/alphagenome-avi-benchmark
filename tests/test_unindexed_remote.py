"""Reading a sorted BGZF file without an index: sparse map, windows, located scans; synthetic data in memory."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import unindexed_remote as ur  # noqa: E402
from test_gnomad_remote import bgzf_block  # noqa: E402

HEADER = "chrom\tpos\tref\talt\tstrand\tdg\tdl\tag\tal\tsummary"


def build(n_pos=4000, step=3, block_payload=3000):
    """Positions 4, 7, ... on chr7, three alleles each; every 5th position also on the minus strand.
    Gene gap: no rows between 6000 and 9000. Split into BGZF blocks at arbitrary byte boundaries."""
    lines = [HEADER]
    pos = 1
    for i in range(n_pos):
        pos += step if pos != 6001 else 3000
        for alt in "CGT":
            lines.append(f"chr7\t{pos}\tA\t{alt}\t+\t0.{i % 10}\t0.1\t0.2\t0.3\t0.{i % 10}")
            if i % 5 == 0:
                lines.append(f"chr7\t{pos}\tA\t{alt}\t-\t0.9\t0.1\t0.2\t0.3\t0.9")
    text = ("\n".join(lines) + "\n").encode()
    data = b"".join(bgzf_block(text[s:s + block_payload]) for s in range(0, len(text), block_payload))
    return data + bgzf_block(b""), text


class Unindexed(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data, cls.text = build()
        cls.fetch = staticmethod(lambda url, a, b: cls.data[a:b])
        cls.rows = [l.split("\t") for l in cls.text.decode().split("\n")[1:] if l]

    def truth(self, positions):
        want = set(positions)
        return sorted(tuple(r) for r in self.rows if int(r[1]) in want)

    def test_windows(self):
        self.assertEqual(ur.windows([10, 12, 500, 30, 12], gap=50), [(10, 30), (500, 500)])

    def test_scan_matches_full_read(self):
        size = len(self.data)
        smap = ur.sparse_map("u", size, "chr7", self.fetch, step=2000, probe=8000, threads=4)
        self.assertGreater(len(smap), 10)
        cases = [[103], [106, 109, 130], [6001], [9001, 9004], [5000], [12097], [11998, 12000, 12097], [7000], [14998], [4]]
        for positions in cases:
            with self.subTest(positions=positions):
                got = ur.fetch_positions("u", size, smap, "chr7", positions, self.fetch, probe=8000, scan=6000)
                self.assertEqual(sorted(tuple(r) for rs in got.values() for r in rs), self.truth(positions))

    def test_many_windows_in_one_call(self):
        size = len(self.data)
        smap = ur.sparse_map("u", size, "chr7", self.fetch, step=20000, probe=8000, threads=4)
        positions = list(range(4, 15000, 297)) + [6001, 9001, 14998]
        got = ur.fetch_positions("u", size, smap, "chr7", positions, self.fetch, probe=8000, scan=6000, gap=50)
        self.assertEqual(sorted(tuple(r) for rs in got.values() for r in rs), self.truth(positions))

    def test_records_that_begin_a_block_are_found(self):
        """Blocks cut exactly at record starts: the first line of a block is complete, and a target there must
        still be returned (the reader skips only lines before a record below the target)."""
        lines = self.text.decode().split("\n")
        header, rows = lines[0], [line for line in lines[1:] if line]
        blocks, starts, chunk, prev = [], [], [header], None
        for row in rows:
            pos = int(row.split("\t")[1])
            if len(chunk) >= 30 and pos != prev:        # cut at the first row of a position
                blocks.append("\n".join(chunk) + "\n")
                starts.append(pos)
                chunk = []
            chunk.append(row)
            prev = pos
        blocks.append("\n".join(chunk) + "\n")
        data = b"".join(bgzf_block(b.encode()) for b in blocks) + bgzf_block(b"")
        fetch = lambda url, a, b: data[a:b]  # noqa: E731
        smap = ur.sparse_map("u", len(data), "chr7", fetch, step=20000, probe=1500, threads=4)
        targets = sorted(set(starts))
        self.assertGreater(len(targets), 50)
        for gap in (0, 50):
            with self.subTest(gap=gap):
                got = ur.fetch_positions("u", len(data), smap, "chr7", targets, fetch, probe=1500, scan=700,
                                         gap=gap)
                self.assertEqual(sorted(tuple(r) for rs in got.values() for r in rs), self.truth(targets))

    def test_bracket_brackets(self):
        smap = [(0, None), (100, 50), (200, 80), (300, 80), (400, 120)]
        self.assertEqual(ur.bracket(smap, 500, 80), (100, 50, 200, 80))
        self.assertEqual(ur.bracket(smap, 500, 10), (0, None, 100, 50))
        self.assertEqual(ur.bracket(smap, 500, 999), (400, 120, 500, None))


if __name__ == "__main__":
    unittest.main()
