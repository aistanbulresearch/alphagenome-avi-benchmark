"""Distance-to-exon classes of the Part 2 variants on a synthetic annotation."""
import gzip
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import part2_positions as pp  # noqa: E402

GTF = "\n".join([
    'chr1\tHAVANA\texon\t100\t200\t.\t+\t.\tgene_id "G1"; gene_type "protein_coding"; gene_name "AAA"; transcript_type "protein_coding";',
    'chr1\tHAVANA\texon\t500\t600\t.\t+\t.\tgene_id "G1"; gene_type "protein_coding"; gene_name "AAA"; transcript_type "protein_coding";',
    'chr1\tHAVANA\texon\t300\t320\t.\t-\t.\tgene_id "G2"; gene_type "protein_coding"; gene_name "BBB"; transcript_type "protein_coding";',
    'chr1\tHAVANA\texon\t900\t950\t.\t+\t.\tgene_id "G3"; gene_type "lncRNA"; gene_name "CCC"; transcript_type "lncRNA";',
]) + "\n"


class Positions(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.gtf = Path(self.tmp.name) / "a.gtf.gz"
        with gzip.open(self.gtf, "wt") as h:
            h.write(GTF)
        self.exons = pp.read_exons(self.gtf)

    def tearDown(self):
        self.tmp.cleanup()

    def test_classes_within_gene(self):
        cases = {150: "exonic", 100: "exonic", 99: "intronic 1-2", 202: "intronic 1-2", 203: "intronic 3-50",
                 260: "intronic >50", 250: "intronic 3-50", 449: "intronic >50", 450: "intronic 3-50"}
        for pos, want in cases.items():
            with self.subTest(pos=pos):
                self.assertEqual(pp.classify(self.exons, "1", pos, "AAA"), want)

    def test_other_gene_exons_ignored_when_gene_known(self):
        self.assertEqual(pp.classify(self.exons, "1", 310, "AAA"), "intronic >50")
        self.assertEqual(pp.classify(self.exons, "1", 310, ""), "exonic")         # no gene: any coding exon
        self.assertEqual(pp.classify(self.exons, "1", 920, ""), "intronic >50")   # lncRNA exons not used

    def test_table(self):
        d = pd.DataFrame({"dataset": ["x", "x", "x"], "chrom": ["1", "1", "1"], "pos": ["150", "250", "99"],
                          "gene": ["AAA", "AAA", "AAA"], "y": [1, 0, 1]})
        t = pp.table(d, self.exons)
        row = t[(t["dataset"] == "x") & (t["class"] == "exonic")].iloc[0]
        self.assertEqual((row["variants"], row["splice_altering"]), (1, 1))


if __name__ == "__main__":
    unittest.main()
