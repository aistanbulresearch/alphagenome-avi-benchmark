"""Every counting reader against the real header of its input file, plus one synthetic row.

The headers were recorded from the downloaded files without reading any data row
(src/record_headers.py; data/derived/header_*.txt)."""
import gzip
import os
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
os.environ["AVI_EVAL_SYNTHETIC"] = "1"

import clinvar_index as idx  # noqa: E402

DERIVED = Path(__file__).resolve().parent.parent / "data" / "derived"


def header(name: str) -> str:
    return (DERIVED / name).read_text(encoding="utf-8")


class RealHeaders(unittest.TestCase):
    def test_submission_summary(self):
        text = header("header_submission_summary.txt")
        self.assertIn("\n#VariationID:", text)                        # a documentation line before the header
        columns = text.rstrip("\n").splitlines()[-1].split("\t")
        row = ["x"] * len(columns)
        row[columns.index("#VariationID")], row[columns.index("Submitter")] = "17", "Lab A"
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "submission_summary.txt.gz"
            with gzip.open(p, "wt", encoding="utf-8", newline="") as h:
                h.write(text + "\t".join(row) + "\n")
            self.assertEqual(list(idx.read_submitters(p, Counter())), [("17", "Lab A")])

    def test_variant_summary(self):
        for name in sorted(p.name for p in DERIVED.glob("header_variant_summary_*.txt")):
            with self.subTest(name=name):
                columns = header(name).rstrip("\n").split("\t")
                row = {c: "-" for c in columns}
                row.update({"#AlleleID": "11", "VariationID": "1", "Assembly": "GRCh38", "Chromosome": "3",
                            "ClinicalSignificance": "Pathogenic", "ReviewStatus": "criteria provided, single submitter",
                            "PositionVCF": "100", "ReferenceAlleleVCF": "A", "AlternateAlleleVCF": "G",
                            "Start": "100", "Stop": "100", "GeneID": "672", "Name": "NM_007294.4(BRCA1):c.5277+50A>G"})
                with tempfile.TemporaryDirectory() as d:
                    p = Path(d) / name.replace("header_", "").replace(".txt", ".txt.gz")
                    with gzip.open(p, "wt", encoding="utf-8", newline="") as h:
                        h.write("\t".join(columns) + "\n" + "\t".join(row[c] for c in columns) + "\n")
                    (rec,) = list(idx.read_variant_summary(p, "T0", Counter()))
                    self.assertEqual((rec.variation_id, rec.key, rec.seen), ("1", "3:100:A:G", True))

    def test_clinvar_vcf(self):
        text = header("header_clinvar_20250615.txt")
        self.assertTrue(text.rstrip("\n").splitlines()[-1].startswith("#CHROM"))
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "clinvar_20250615.vcf.gz"
            with gzip.open(p, "wt", encoding="utf-8", newline="") as h:
                h.write(text + "3\t100\t1\tA\tG\t.\t.\tALLELEID=11;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter\n")
            (rec,) = list(idx.read_vcf(p, "E-dev", Counter()))
            self.assertEqual((rec.variation_id, rec.key, rec.seen), ("1", "3:100:A:G", True))


if __name__ == "__main__":
    unittest.main()
