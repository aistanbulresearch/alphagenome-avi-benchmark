"""Tests for the small derived exposure files, on synthetic inputs."""
import contextlib
import gzip
import hashlib
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
os.environ["AVI_EVAL_SYNTHETIC"] = "1"

import build_exposure_files as be  # noqa: E402
import development_sources as s1  # noqa: E402
import clinvar_index as idx  # noqa: E402

ORG_HEADER = "#organization\torganization ID\tinstitution type\tcity\n"
BROAD = "Broad Center for Mendelian Genomics; Broad Institute of MIT and Harvard"


class Submitters(unittest.TestCase):
    def test_file_lists_the_organization_and_both_spellings(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            (d / "orgs.txt").write_text(ORG_HEADER + f"Other Lab\t1\tlab\tX\n{BROAD}\t506627\tlab\tCambridge\n",
                                        encoding="utf-8")
            be.write_submitters(d / "orgs.txt", d / "subs.tsv")
            text = (d / "subs.tsv").read_text(encoding="utf-8")
            self.assertEqual(text, "organization_id\tsubmitter_name\n"
                                   "506627\tBroad Center for Mendelian Genomics, Broad Institute of MIT and Harvard\n"
                                   f"506627\t{BROAD}\n")
            with contextlib.nullcontext():
                self.assertIn("506627", text)

    def test_organization_must_exist_with_the_listed_name(self):
        cases = {"absent": ORG_HEADER + "Other Lab\t1\tlab\tX\n",
                 "renamed": ORG_HEADER + "Broad CMG\t506627\tlab\tX\n",
                 "no header": f"{BROAD}\t506627\tlab\tX\n",
                 "listed twice": ORG_HEADER + f"{BROAD}\t506627\tlab\tX\n{BROAD}\t506627\tlab\tY\n"}
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            for label, text in cases.items():
                with self.subTest(label=label):
                    (d / "orgs.txt").write_text(text, encoding="utf-8")
                    with self.assertRaises(idx.InputFormatError):
                        be.write_submitters(d / "orgs.txt", d / f"subs_{label}.tsv")


class Chain(unittest.TestCase):
    def test_copied_unchanged(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            data = gzip.compress(b"chain 1000 chr1 100 + 0 100 chr1 100 + 0 100 1\n100\n\n")
            (d / "hg19ToHg38.over.chain.gz").write_bytes(data)
            digest = be.copy_chain(d / "hg19ToHg38.over.chain.gz", d / "out" / "hg19ToHg38.over.chain.gz")
            self.assertEqual((d / "out" / "hg19ToHg38.over.chain.gz").read_bytes(), data)
            self.assertEqual(digest, hashlib.sha256(data).hexdigest())
            with self.assertRaises(FileExistsError):
                be.copy_chain(d / "hg19ToHg38.over.chain.gz", d / "out" / "hg19ToHg38.over.chain.gz")


class Sources(unittest.TestCase):
    def test_coverage_file_passes_the_counting_check(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            evidence = {}
            for source in s1.REQUIRED_DEV_SOURCES:
                name = f"evidence/{source.replace(' ', '_')}.json"
                (d / name).parent.mkdir(exist_ok=True)
                (d / name).write_text("{}", encoding="utf-8")
                evidence[source] = name
            be.write_sources(evidence, d / "development_sources.tsv")
            with contextlib.nullcontext():
                sources = s1.read_sources(d / "development_sources.tsv")
            self.assertEqual(s1.coverage_problems(sources, d), [])
            self.assertEqual(sources["GREGoR"][0], "submitter rule")

    def test_every_source_needs_an_existing_evidence_file(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            evidence = {s: "missing.json" for s in s1.REQUIRED_DEV_SOURCES}
            with self.assertRaises(idx.InputFormatError):
                be.write_sources(evidence, d / "development_sources.tsv")
            del evidence["GREGoR"]
            with self.assertRaises(idx.InputFormatError):
                be.write_sources(evidence, d / "development_sources.tsv")
            self.assertFalse((d / "development_sources.tsv").exists())


if __name__ == "__main__":
    unittest.main()
