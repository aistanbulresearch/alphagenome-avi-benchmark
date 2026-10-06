"""Tests for the MANE Select transcript table, on synthetic rows."""
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

import build_mane_table as bm  # noqa: E402
import development_sources as s1  # noqa: E402
import variant_sets as vsets  # noqa: E402
import clinvar_index as idx  # noqa: E402

HEADER = ("#NCBI_GeneID\tEnsembl_Gene\tHGNC_ID\tsymbol\tname\tRefSeq_nuc\tRefSeq_prot\tEnsembl_nuc\tEnsembl_prot\t"
          "MANE_status\tGRCh38_chr\tchr_start\tchr_end\tchr_strand\n")


def row(gene, transcript, status):
    return (f"GeneID:{gene}\tENSG{gene}.1\tHGNC:{gene}\tG{gene}\tgene {gene}\t{transcript}\tNP_{gene}.1\t"
            f"ENST{gene}.1\tENSP{gene}.1\t{status}\tNC_000003.12\t100\t200\t+\n")


def write(path, rows, header=HEADER):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as h:
        h.write(header + "".join(rows))


class ManeTable(unittest.TestCase):
    def test_select_rows_with_plus_clinical_flag(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            write(d / "mane.gz", [row(1, "NM_000001.4", "MANE Select"), row(2, "NM_000002.1", "MANE Select"),
                                  row(2, "NM_000202.3", "MANE Plus Clinical")])
            counts = bm.build(d / "mane.gz", d / "mane_select.tsv")
            self.assertEqual((d / "mane_select.tsv").read_text(encoding="utf-8"),
                             "refseq_transcript\tgene_id\tplus_clinical\n"
                             "NM_000001.4\t1\tno\nNM_000002.1\t2\tyes\n")
            self.assertEqual(counts, {"MANE Select": 2, "MANE Plus Clinical": 1, "genes with Plus Clinical": 1})
            with contextlib.nullcontext():
                self.assertEqual(vsets.read_mane(d / "mane_select.tsv"), {"NM_000001.4", "NM_000002.1"})

    def test_malformed_input_fails_closed(self):
        cases = {
            "missing column": ([row(1, "NM_000001.4", "MANE Select")], HEADER.replace("\tMANE_status", "")),
            "unknown status": ([row(1, "NM_000001.4", "MANE Other")], HEADER),
            "unversioned transcript": ([row(1, "NM_000001", "MANE Select")], HEADER),
            "two Select transcripts for one gene": ([row(1, "NM_000001.4", "MANE Select"),
                                                     row(1, "NM_000011.1", "MANE Select")], HEADER),
            "Plus Clinical without Select": ([row(3, "NM_000003.1", "MANE Plus Clinical")], HEADER),
            "no rows": ([], HEADER),
        }
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            for label, (rows, header) in cases.items():
                with self.subTest(label=label):
                    write(d / "mane.gz", rows, header)
                    with self.assertRaises(idx.InputFormatError):
                        bm.build(d / "mane.gz", d / "out.tsv")

    def test_existing_output_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            write(d / "mane.gz", [row(1, "NM_000001.4", "MANE Select")])
            (d / "out.tsv").write_text("x", encoding="utf-8")
            with self.assertRaises(FileExistsError):
                bm.build(d / "mane.gz", d / "out.tsv")


if __name__ == "__main__":
    unittest.main()
