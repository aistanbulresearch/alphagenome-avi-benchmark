"""Tests for the ClinVar replacement history, on synthetic XML."""
import contextlib
import gzip
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
os.environ["AVI_EVAL_SYNTHETIC"] = "1"

import build_replacements as br  # noqa: E402
import clinvar_index as idx  # noqa: E402


def archive(vid, status="current", replaced_by=None, replaced=(), nested=()):
    """One VariationArchive; `nested` adds an RCV-level ReplacedList that must be ignored."""
    parts = [f'<VariationArchive VariationID="{vid}" Accession="VCV{int(vid):09d}" Version="1">',
             f"<RecordStatus>{status}</RecordStatus>"]
    if replaced_by:
        parts.append(f'<ReplacedBy Accession="VCV{int(replaced_by):09d}" Version="1" DateChanged="2026-01-01" '
                     f'VariationID="{replaced_by}"/>')
    if replaced:
        parts.append("<ReplacedList>" + "".join(r for r in replaced) + "</ReplacedList>")
    parts.append("<Species>Homo sapiens</Species><ClassifiedRecord><RCVList><RCVAccession>")
    if nested:
        parts.append("<ReplacedList>" + "".join(nested) + "</ReplacedList>")
    parts.append("</RCVAccession></RCVList></ClassifiedRecord></VariationArchive>")
    return "\n".join(parts)


def rep(vid=None, accession=None):
    acc = accession or f"VCV{int(vid):09d}"
    attr = f' VariationID="{vid}"' if vid is not None else ""
    return f'<Replaced Accession="{acc}" Version="2" DateChanged="2025-05-05"{attr}/>'


def write(path, archives):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as h:
        h.write('<?xml version="1.0" encoding="UTF-8" ?>\n<ClinVarVariationRelease>\n' + "\n".join(archives)
                + "\n</ClinVarVariationRelease>\n")


class Replacements(unittest.TestCase):
    def test_pairs_from_replaced_lists_and_replaced_by(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            write(d / "vcv.xml.gz", [
                archive(100, replaced=[rep(7), rep(8), rep(accession="VCV000000009")],
                        nested=[rep(999, accession="RCV000000999")]),
                archive(8, status="replaced", replaced_by=100),
                archive(200, status="replaced", replaced_by=201),
                archive(201),
            ])
            stats = br.build(d / "vcv.xml.gz", d / "rep.tsv")
            self.assertEqual((d / "rep.tsv").read_text(encoding="utf-8"),
                             "replaced_variation_id\tcurrent_variation_id\n7\t100\n8\t100\n9\t100\n200\t201\n")
            self.assertEqual(stats["records"], 4)
            self.assertEqual(stats["pairs"], 4)
            with contextlib.nullcontext():
                self.assertEqual(idx.read_replacements(d / "rep.tsv"), {"7": "100", "8": "100", "9": "100", "200": "201"})

    def test_conflicts_and_malformed_entries_fail_closed(self):
        cases = {
            "two current records claim one replaced id": [archive(100, replaced=[rep(7)]), archive(300, replaced=[rep(7)])],
            "accession and VariationID disagree": [archive(100, replaced=[rep(7, accession="VCV000000006")])],
            "replaced record without ReplacedBy": [archive(5, status="replaced")],
            "Replaced without any identifier": [archive(100, replaced=[rep(accession="SCV000000001")])],
            "no record": [],
        }
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            for label, archives in cases.items():
                with self.subTest(label=label):
                    write(d / "vcv.xml.gz", archives)
                    with self.assertRaises(idx.InputFormatError):
                        br.build(d / "vcv.xml.gz", d / f"out_{len(label)}.tsv")

    def test_existing_output_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            write(d / "vcv.xml.gz", [archive(100, replaced=[rep(7)])])
            (d / "rep.tsv").write_text("x", encoding="utf-8")
            with self.assertRaises(FileExistsError):
                br.build(d / "vcv.xml.gz", d / "rep.tsv")


if __name__ == "__main__":
    unittest.main()
