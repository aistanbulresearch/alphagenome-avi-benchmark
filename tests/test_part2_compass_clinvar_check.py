"""COMPASS ClinVar provenance check on synthetic data (labels only, no scores)."""
import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import part2_compass_clinvar_check as cc  # noqa: E402

HEADER = "#AlleleID\tAssembly\tChromosome\tPositionVCF\tReferenceAlleleVCF\tAlternateAlleleVCF\n"


class CompassClinvarCheck(unittest.TestCase):
    def test_counts_by_subset_only_grch38_exact_alleles(self):
        keys = {("1", 100, "A", "G"): "random_snv", ("1", 200, "C", "T"): "exac", ("X", 5, "G", "A"): "random_snv"}
        lines = ["1\tGRCh38\t1\t100\tA\tG\n",        # random, present
                 "2\tGRCh37\t1\t200\tC\tT\n",        # wrong assembly
                 "3\tGRCh38\t1\t200\tC\tA\n",        # other allele
                 "4\tGRCh38\tX\t5\tG\tA\n",          # random, present
                 "5\tGRCh38\t1\t200\tC\tT\n",        # exac, present
                 "6\tGRCh38\t1\tna\tC\tT\n"]         # unparsable position
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "variant_summary.txt.gz"
            with gzip.open(path, "wt", encoding="utf-8") as h:
                h.write(HEADER + "".join(lines))
            self.assertEqual(cc.present_by_subset(keys, path), {"random_snv": 2, "exac": 1})

    def test_significance_of_present_alleles_counted_once_per_allele(self):
        keys = {("1", 100, "A", "G"): "random_snv", ("1", 200, "C", "T"): "random_snv", ("2", 7, "G", "A"): "exac"}
        header = HEADER.rstrip("\n") + "\tClinicalSignificance\n"
        lines = ["1\tGRCh38\t1\t100\tA\tG\tUncertain significance\n",
                 "2\tGRCh38\t1\t100\tA\tG\tUncertain significance\n",     # same allele, second record
                 "3\tGRCh38\t1\t200\tC\tT\t-\n",
                 "4\tGRCh38\t2\t7\tG\tA\tBenign\n"]                      # other subset
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "variant_summary.txt.gz"
            with gzip.open(path, "wt", encoding="utf-8") as h:
                h.write(header + "".join(lines))
            self.assertEqual(cc.significance_of(keys, path, "random_snv"), {"Uncertain significance": 1, "-": 1})

    def test_significance_reads_every_row_and_shows_conflicting_rows(self):
        keys = {("1", 100, "A", "G"): "random_snv", ("1", 200, "C", "T"): "random_snv"}
        header = HEADER.rstrip("\n") + "\tClinicalSignificance\n"
        lines = ["1\tGRCh38\t1\t100\tA\tG\tUncertain significance\n",
                 "2\tGRCh38\t1\t100\tA\tG\tLikely benign\n",        # same allele, a later row disagrees
                 "3\tGRCh38\t1\t200\tC\tT\tUncertain significance\n"]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "variant_summary.txt.gz"
            with gzip.open(path, "wt", encoding="utf-8") as h:
                h.write(header + "".join(lines))
            self.assertEqual(cc.significance_of(keys, path, "random_snv"),
                             {"Likely benign | Uncertain significance": 1, "Uncertain significance": 1})

    def test_analysed_keys_drop_part1_and_benchmark_alleles(self):
        keys = {("1", 100, "A", "G"): "random_snv", ("1", 200, "C", "T"): "random_snv",
                ("2", 7, "G", "A"): "exac", ("3", 9, "T", "C"): "random_snv"}
        out = cc.analysed_keys(keys, part1={"1:200:C:T"}, development={"2:7:G:A"})
        self.assertEqual(out, {("1", 100, "A", "G"): "random_snv", ("3", 9, "T", "C"): "random_snv"})

    def test_chromosome_names_are_normalised(self):
        self.assertEqual(cc.norm("chr7"), "7")
        self.assertEqual(cc.norm("chrM"), "MT")
        self.assertEqual(cc.norm("X"), "X")


if __name__ == "__main__":
    unittest.main()
