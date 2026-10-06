"""Tests for the download record (synthetic files only)."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import download_public as dp  # noqa: E402


class Record(unittest.TestCase):
    def test_ledger_entry(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            source = root / "src_file.txt"
            source.write_bytes(b"abc\n")
            with mock.patch.object(dp, "STUDY_ROOT", root), mock.patch.object(dp, "RAW_DIR", root / "data" / "raw"), \
                    mock.patch.object(dp, "LEDGER", root / "data" / "downloads.json"):
                self.assertEqual(dp.main([source.as_uri()]), 0)
                ledger = json.loads((root / "data" / "downloads.json").read_text(encoding="utf-8"))
                self.assertEqual(ledger[-1]["file"], "data/raw/src_file.txt")
                self.assertEqual(ledger[-1]["bytes"], 4)
                self.assertNotIn(b"\r", (root / "data" / "downloads.json").read_bytes())   # LF only
                with self.assertRaises(FileExistsError):
                    dp.main([source.as_uri()])                                # never overwrites


class TargetName(unittest.TestCase):
    def test_name_places_the_file_inside_raw(self):
        with tempfile.TemporaryDirectory() as t:
            root = Path(t)
            source = root / "scores"
            source.write_bytes(b"x\n")
            with mock.patch.object(dp, "STUDY_ROOT", root), mock.patch.object(dp, "RAW_DIR", root / "data" / "raw"), \
                    mock.patch.object(dp, "LEDGER", root / "data" / "downloads.json"):
                dp.main(["--name", "mavedb/urn_mavedb_00000675-a-1_scores.csv",
                         source.as_uri()])
                ledger = json.loads((root / "data" / "downloads.json").read_text(encoding="utf-8"))
                self.assertEqual(ledger[-1]["file"], "data/raw/mavedb/urn_mavedb_00000675-a-1_scores.csv")
                for bad in ("../outside.csv", "/abs.csv", "mavedb/../../x.csv", "C:/x.csv"):
                    with self.assertRaises(ValueError, msg=bad):
                        dp.main(["--name", bad, source.as_uri()])
                with self.assertRaises(ValueError):                    # one name, one URL
                    dp.main(["--name", "a.csv", source.as_uri(), source.as_uri()])


if __name__ == "__main__":
    unittest.main()
