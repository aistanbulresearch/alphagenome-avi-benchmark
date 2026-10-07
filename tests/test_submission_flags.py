"""Presence-only flags for splicing predictors in submission texts (synthetic texts)."""
import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import submission_flags as sf  # noqa: E402


class Flags(unittest.TestCase):
    def test_terms(self):
        self.assertEqual(sf.flags(["SpliceAI delta score supports PP3"]), (True, True))
        self.assertEqual(sf.flags(["Splice AI predicts no effect"]), (True, True))
        self.assertEqual(sf.flags(["MaxEntScan and Pangolin predict loss of the donor"]), (False, True))
        self.assertEqual(sf.flags(["in silico tools predict an impact on splicing"]), (False, True))
        self.assertEqual(sf.flags(["predicted to disrupt normal splicing"]), (False, True))
        self.assertEqual(sf.flags(["observed in homozygous state in healthy controls"]), (False, False))
        self.assertEqual(sf.flags([]), (False, False))

    def test_reader_uses_the_header_and_wanted_ids(self):
        header = ("#VariationID\tClinicalSignificance\tDateLastEvaluated\tDescription\tSubmittedPhenotypeInfo\t"
                  "ReportedPhenotypeInfo\tReviewStatus\tCollectionMethod\tOriginCounts\tSubmitter\tSCV\t"
                  "SubmittedGeneSymbol\tExplanationOfInterpretation\tSomaticClinicalImpact\tOncogenicity\t"
                  "ContributesToAggregateClassification")
        rows = ["1\tPathogenic\t-\tSpliceAI 0.9\t-\t-\t-\t-\t-\tLab\tSCV1\tG\t-\t-\t-\tyes",
                "2\tBenign\t-\t-\t-\t-\t-\t-\t-\tLab\tSCV2\tG\tcommon in population\t-\t-\tyes"]
        with tempfile.TemporaryDirectory() as t:
            path = Path(t) / "s.gz"
            with gzip.open(path, "wt", encoding="utf-8") as h:
                h.write("##Overview\n#VariationID: the identifier\n" + header + "\n" + "\n".join(rows) + "\n")
            got = list(sf.read_submissions(path, {"1"}))
        self.assertEqual(got, [("1", "SpliceAI 0.9", "-")])


if __name__ == "__main__":
    unittest.main()
