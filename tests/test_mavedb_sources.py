"""Tests for the MaveDB part of the development positions, synthetic files."""
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
os.environ["AVI_EVAL_SYNTHETIC"] = "1"

import mavedb_sources as ms  # noqa: E402
import clinvar_index as idx  # noqa: E402

CHROMS = {"NC_000017.11": "17", "NC_000003.12": "3"}


def allele(label, start, end):
    return {"type": "Allele", "location": {"start": start, "end": end, "sequenceReference": {"label": label}}}


def mapped(urn, post, current=True):
    return {"variantUrn": urn, "postMapped": post, "current": current}


class Versions(unittest.TestCase):
    def test_every_version_of_each_experiment(self):
        records = {
            ("experiment", "urn:mavedb:1-a"): {"scoreSetUrns": ["urn:mavedb:1-a-3"]},
            ("score-set", "urn:mavedb:1-a-3"): {"supersededScoreSet": {"urn": "urn:mavedb:1-a-2"}},
            ("score-set", "urn:mavedb:1-a-2"): {"supersededScoreSet": "urn:mavedb:1-a-1"},
            ("score-set", "urn:mavedb:1-a-1"): {"supersededScoreSet": None},
        }
        got = ms.score_set_urns("urn:mavedb:1-a", lambda kind, urn: records[(kind, urn)])
        self.assertEqual(got, ["urn:mavedb:1-a-1", "urn:mavedb:1-a-2", "urn:mavedb:1-a-3"])

    def test_studies_follow_the_source_table(self):
        genes = {s["gene"] for s in ms.STUDIES}
        self.assertEqual(genes, {"BRCA1", "BRCA2", "BAP1", "RAD51C", "VHL", "DDX3X", "RNU4-2", "ATM", "BARD1",
                                 "PALB2", "XRCC2", "SFPQ", "CTCF"})
        self.assertEqual(len(ms.STUDIES), 14)
        region_now = sorted(s["gene"] for s in ms.STUDIES if s["experiments"] is None)
        self.assertEqual(region_now, ["ATM", "BRCA1", "RNU4-2"])            # De Jonghe, Lee, Dace
        findlay = next(s for s in ms.STUDIES if s["study"] == "Findlay 2018")
        self.assertEqual(len(findlay["experiments"]), 27)                  # -0 and -a to -z


class Fetch(unittest.TestCase):
    def test_unavailable_records_are_listed_and_the_fetch_continues(self):
        # A superseded version that MaveDB no longer serves (404) is recorded as unavailable; the
        # builder then applies the gene region to that study. Any other HTTP error stops the fetch.
        import urllib.error
        served = {
            ms.local_name("urn:mavedb:1-a", "experiment"): {"scoreSetUrns": ["urn:mavedb:1-a-2"]},
            ms.local_name("urn:mavedb:1-a-2", "score-set"): {"supersededScoreSet": {"urn": "urn:mavedb:1-a-1"}},
            ms.local_name("urn:mavedb:1-a-2", "scores"): "accession\n",
            ms.local_name("urn:mavedb:1-a-2", "mapped"): [],
            ms.local_name("urn:mavedb:2-a", "experiment"): {"scoreSetUrns": ["urn:mavedb:2-a-1"]},
            ms.local_name("urn:mavedb:2-a-1", "score-set"): {"supersededScoreSet": None},
            ms.local_name("urn:mavedb:2-a-1", "scores"): "accession\n",
            ms.local_name("urn:mavedb:2-a-1", "mapped"): [],
        }
        with tempfile.TemporaryDirectory() as t:
            raw = Path(t)

            def download(url, name):
                if name not in served:
                    raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)
                (raw / name).parent.mkdir(parents=True, exist_ok=True)
                body = served[name]
                (raw / name).write_text(body if isinstance(body, str) else json.dumps(body), encoding="utf-8")

            studies = ({"study": "A", "gene": "G1", "experiments": ("urn:mavedb:1-a",)},
                       {"study": "B", "gene": "G2", "experiments": ("urn:mavedb:2-a",)})
            with mock.patch.multiple(ms, RAW=raw, STUDIES=studies):
                unavailable = ms.fetch(download)
                self.assertEqual(unavailable, [ms.local_name("urn:mavedb:1-a-1", "score-set")])
                self.assertTrue((raw / ms.local_name("urn:mavedb:2-a-1", "mapped")).exists())   # continued
                with self.assertRaises(ms.MissingSource):                                       # gene region
                    ms.study_score_sets(studies[0])

                def failing(url, name):
                    raise urllib.error.HTTPError(url, 500, "Server Error", {}, None)
                with self.assertRaises(urllib.error.HTTPError):
                    with mock.patch.multiple(ms, RAW=Path(t) / "other"):
                        ms.fetch(failing)


class Retry(unittest.TestCase):
    def test_rate_limit_is_retried_with_backoff(self):
        import urllib.error
        calls, slept = [], []

        def flaky(url, name):
            calls.append(name)
            if len(calls) < 3:
                raise urllib.error.HTTPError(url, 429, "Too Many Requests", {"Retry-After": "7"}, None)

        ms.with_retry(flaky, sleep=slept.append, attempts=5)("u", "n")
        self.assertEqual(len(calls), 3)
        self.assertEqual(slept[:2], [7, 7])                       # honours Retry-After

    def test_gateway_errors_and_timeouts_are_retried(self):
        import urllib.error
        for failure in (urllib.error.HTTPError("u", 504, "Gateway Time-out", {}, None),
                        urllib.error.HTTPError("u", 502, "Bad Gateway", {}, None), TimeoutError("read timed out"),
                        urllib.error.URLError("connection reset")):
            calls = []

            def flaky(url, name):
                calls.append(name)
                if len(calls) < 2:
                    raise failure

            ms.with_retry(flaky, sleep=lambda s: None, attempts=4)("u", "n")
            self.assertEqual(len(calls), 2, failure)

    def test_other_errors_and_exhausted_retries_raise(self):
        import urllib.error

        def not_found(url, name):
            raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)

        def always_busy(url, name):
            raise urllib.error.HTTPError(url, 503, "Busy", {}, None)

        with self.assertRaises(urllib.error.HTTPError):
            ms.with_retry(not_found, sleep=lambda s: None, attempts=5)("u", "n")
        with self.assertRaises(urllib.error.HTTPError):
            ms.with_retry(always_busy, sleep=lambda s: None, attempts=3)("u", "n")


class Inventory(unittest.TestCase):
    def test_exact_experiments_of_the_source_table(self):
        # specified here independently from the module.
        def letters(base, last):
            return [f"urn:mavedb:{base}-{c}" for c in "abcdefghijklmnopqrstuvwxyz"[:"abcdefghijklmnopqrstuvwxyz".index(last) + 1]]
        expected = {
            "Findlay 2018": ["urn:mavedb:00000097-0"] + letters("00000097", "z"),
            "Huang 2025": ["urn:mavedb:00001225-a", "urn:mavedb:00001279-a"],
            "Waters 2024": ["urn:mavedb:00000662-0"] + letters("00000662", "r"),
            "Olvera-Leon 2024": ["urn:mavedb:00000673-0"] + letters("00000673", "i"),
            "Buckley 2024": ["urn:mavedb:00000675-a"],
            "Radford 2023": ["urn:mavedb:00000658-0"] + letters("00000658", "q") + ["urn:mavedb:00001228-a"],
            "De Jonghe 2026": None, "Lee 2025": None, "Dace 2025": None,
            "Woo 2025": ["urn:mavedb:00001250-a"],
            "Tejura 2026 PALB2": ["urn:mavedb:00001259-a"], "Tejura 2026 XRCC2": ["urn:mavedb:00001264-a"],
            "Tejura 2026 SFPQ": ["urn:mavedb:00001265-a"], "Tejura 2026 CTCF": ["urn:mavedb:00001262-a"],
        }
        got = {s["study"]: (sorted(s["experiments"]) if s["experiments"] is not None else None) for s in ms.STUDIES}
        self.assertEqual(got, {k: (sorted(v) if v is not None else None) for k, v in expected.items()})


class Positions(unittest.TestCase):
    def test_positions_from_mapped_variants(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            (d / "s.scores.csv").write_text("accession,hgvs_nt,score\nurn:mavedb:1-a-1#1,c.1A>G,0.5\n"
                                            "urn:mavedb:1-a-1#2,c.2A>G,-1\nurn:mavedb:1-a-1#3,c.3_4del,NA\n",
                                            encoding="utf-8")
            (d / "s.mapped.json").write_text(json.dumps([
                mapped("urn:mavedb:1-a-1#1", allele("NC_000017.11", 99, 100)),             # SNV at 100
                mapped("urn:mavedb:1-a-1#1", allele("NC_000017.10", 5, 6), current=False),  # other assembly
                mapped("urn:mavedb:1-a-1#2", {"type": "Haplotype", "members": [allele("NC_000017.11", 199, 200),
                                                                               allele("NC_000017.11", 209, 210)]}),
                mapped("urn:mavedb:1-a-1#3", allele("NC_000017.11", 299, 301)),            # 2-bp deletion
            ]), encoding="utf-8")
            result = ms.score_set_positions(d / "s.scores.csv", d / "s.mapped.json", CHROMS)
            self.assertEqual(result["unpositioned"], 0)
            self.assertEqual(sorted(result["intervals"]),
                             [("17", 100, 100), ("17", 200, 200), ("17", 210, 210), ("17", 300, 301)])
            self.assertEqual(result["forms"], {"Allele": 3, "Haplotype": 1})

    def test_insertion_spans_both_flanking_bases(self):
        self.assertEqual(ms.interval(allele("NC_000003.12", 500, 500), CHROMS), ("3", 500, 501))

    def test_a_scored_variant_without_a_grch38_position_is_counted(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            (d / "s.scores.csv").write_text("accession,score\nurn:x#1,1\nurn:x#2,1\nurn:x#3,1\n", encoding="utf-8")
            (d / "s.mapped.json").write_text(json.dumps([
                mapped("urn:x#1", allele("NC_000017.11", 9, 10)),
                mapped("urn:x#2", allele("ENST000001.1", 9, 10)),                # transcript only
                mapped("urn:x#3", None),                                          # mapping failed
            ]), encoding="utf-8")
            result = ms.score_set_positions(d / "s.scores.csv", d / "s.mapped.json", CHROMS)
            self.assertEqual(result["unpositioned"], 2)
            self.assertEqual(result["intervals"], [("17", 10, 10)])

    def test_rows_without_a_usable_accession_count_as_unpositioned(self):
        # an empty or short row is a scored variant that cannot be located.
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            (d / "s.scores.csv").write_text("accession,hgvs_nt,score\nurn:x#1,c.1A>G,1\n,c.2A>G,2\nurn:x#9\n",
                                            encoding="utf-8")
            (d / "s.mapped.json").write_text(json.dumps([mapped("urn:x#1", allele("NC_000017.11", 9, 10))]),
                                             encoding="utf-8")
            result = ms.score_set_positions(d / "s.scores.csv", d / "s.mapped.json", CHROMS)
            self.assertEqual((result["scored"], result["unpositioned"]), (3, 2))     # empty row and unmapped #9

    def test_partly_mapped_haplotype_is_unpositioned(self):
        # every member of a haplotype must have a GRCh38 location.
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            (d / "s.scores.csv").write_text("accession,score\nurn:x#1,1\n", encoding="utf-8")
            (d / "s.mapped.json").write_text(json.dumps([mapped("urn:x#1", {"type": "Haplotype", "members": [
                allele("NC_000003.12", 99, 100), allele("ENST000001.1", 5, 6)]})]), encoding="utf-8")
            result = ms.score_set_positions(d / "s.scores.csv", d / "s.mapped.json", CHROMS)
            self.assertEqual(result["unpositioned"], 1)
            self.assertIn(("3", 100, 100), result["intervals"])            # still excluded (superset)

    def test_missing_metadata_is_reported_not_fatal(self):
        # an experiment or score-set record that was never obtained routes the study to its
        # gene region; a recorded file that fails its provenance check stays fatal.
        with tempfile.TemporaryDirectory() as t:
            with mock.patch.object(ms, "RAW", Path(t)):
                study = {"study": "X", "gene": "VHL", "experiments": ("urn:mavedb:5-a",)}
                with self.assertRaises(ms.MissingSource):
                    ms.study_score_sets(study)

    def test_scores_file_needs_an_accession_column(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            (d / "s.scores.csv").write_text("hgvs_nt,score\nc.1A>G,1\n", encoding="utf-8")
            (d / "s.mapped.json").write_text("[]", encoding="utf-8")
            with self.assertRaises(idx.InputFormatError):
                ms.score_set_positions(d / "s.scores.csv", d / "s.mapped.json", CHROMS)

    def test_local_names_are_unique_and_safe(self):
        names = {ms.local_name("urn:mavedb:00000097-a-2", k) for k in ("score-set", "scores", "mapped")}
        names.add(ms.local_name("urn:mavedb:00000097-a", "experiment"))
        self.assertEqual(len(names), 4)
        self.assertTrue(all(n.startswith("mavedb/") and ":" not in n and ".." not in n for n in names))


if __name__ == "__main__":
    unittest.main()
