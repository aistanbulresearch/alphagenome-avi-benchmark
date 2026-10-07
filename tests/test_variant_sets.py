"""End-to-end build of the ClinVar variant list on synthetic gzipped files in a temporary folder."""
import gzip
import json
import os
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import clinvar_records as cr  # noqa: E402
import hgvs_strata as hs  # noqa: E402
import variant_sets as vsets  # noqa: E402
import clinvar_index as idx  # noqa: E402

HEADER = (Path(__file__).resolve().parent.parent / "data" / "derived"
          / "header_variant_summary_2026-10.txt").read_text(encoding="utf-8").strip().split("\t")
REV = {1: "criteria provided, single submitter", 2: "criteria provided, multiple submitters, no conflicts"}
DEEP = "NM_000001.1(G1):c.100+60A>G"


def vs(vid, aid, clnsig, chrom="3", pos=1000, ref="A", alt="G", stars=2, assembly="GRCh38",
       name=DEEP, gene="101", review=None):
    row = {h: "-" for h in HEADER}
    row.update({"#AlleleID": str(aid), "VariationID": str(vid), "ClinicalSignificance": clnsig,
                "ReviewStatus": review or REV[stars], "Assembly": assembly, "Chromosome": chrom,
                "PositionVCF": str(pos), "ReferenceAlleleVCF": ref, "AlternateAlleleVCF": alt,
                "Name": name, "GeneID": gene})
    return "\t".join(row[h] for h in HEADER)


def write_vs(path, rows, header=HEADER):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as h:
        h.write("\t".join(header) + "\n" + "".join(r + "\n" for r in rows))


def write_vcf(path, lines, header=True):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as h:
        h.write("##fileformat=VCFv4.1\n")
        if header:
            h.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        h.write("".join(l + "\n" for l in lines))


T1_ROWS = [
    vs(1, 11, "Pathogenic", pos=1001),                                        # new label
    vs(1, 11, "Pathogenic", assembly="GRCh37", pos=901),                      # GRCh37 twin
    vs(2, 12, "Likely benign", pos=1002),                                     # seen in T0
    vs(3, 13, "Pathogenic", pos=1003),                                        # re-entry (E-hist only)
    vs(4, 14, "Pathogenic", pos=1004), vs(4, 14, "Benign", pos=1004),         # ambiguous duplicate
    vs(5, 15, "Pathogenic", assembly="GRCh37"),                               # no GRCh38 row
    vs(6, 16, "Pathogenic", chrom="Y", pos=1006),                             # excluded chromosome
    vs(7, 17, "Benign", pos=1007, review="something new"),                    # unknown star string
    vs(8, 18, "Benign", pos=1008),                                            # collision with E-dev
    vs(9, 19, "Uncertain significance", pos=1009),                            # not strict
    vs(10, 20, "Benign", pos=1010, stars=1, chrom="2"),                       # validation chromosome, 1 star
    vs(11, 21, "Pathogenic", pos=1011, ref="CA", alt="C", name="NM_000001.1(G1):c.100+70del"),
    vs(12, 22, "Benign", pos=1012, name="GRCh38/hg38 3p21(chr3:1-2)x3"),     # unparsed
    vs(13, 23, "Pathogenic", pos=1013, review="no assertion criteria provided"),  # 0 stars
    vs(14, 24, "Pathogenic", pos=1014, name="NM_000001.1(G1):c.100+7A>G"),   # offset +7
    vs(15, 25, "Pathogenic", pos=1015, name="NM_000001.1(G1):c.100+9A>G"),   # offset +9
    vs(16, 26, "Benign", pos=1016),                                           # development position
    vs(17, 27, "Pathogenic", pos=1017),                                       # development submitter
    vs(18, 28, "Pathogenic", pos=1018), vs(18, 28, "Pathogenic", assembly="GRCh37", pos=918),
    vs(19, 29, "Pathogenic", pos=1019),                                       # E-dev line without ALT
    vs(20, 30, "Pathogenic", pos=1020),                                       # seen via non-strict GRCh37 twin
    vs(20, 30, "Uncertain significance", assembly="GRCh37", pos=920),
    vs(21, 31, "Pathogenic", pos=1021),                                       # blocked: unknown record at position
    vs(22, 32, "Pathogenic", pos=1022),                                       # seen via replacement history
    vs(23, -1, "Uncertain significance", pos=1023),                           # missing AlleleID in T1
    vs(24, 34, "Pathogenic", pos=1024),                                       # blocked: unknown E-dev line, no ALT
    vs(25, 35, "Pathogenic", pos=1025, name="NM_000002.1(G2):c.100+60A>G"),   # not a MANE transcript
    vs(26, 36, "Pathogenic", pos=1026, name="NM_000001.2(G1):c.100+60A>G"),   # MANE accession, other version
]
T0_ROWS = [vs(2, 12, "Benign/Likely benign", pos=1002), vs(3, 13, "Uncertain significance", pos=1003),
           vs(918, 918, "Pathogenic", assembly="GRCh37", pos=918),            # merged: GRCh37 allele only
           vs(920, 920, "Likely pathogenic", assembly="GRCh37", pos=920),     # matches a non-strict twin
           vs(7022, 7022, "Pathogenic", pos="-1", ref="na", alt="na"),        # replaced by VID 22
           vs(9021, 9021, "Pathogenic", pos=1021, alt="T"),                   # unknown at T1, placed
           vs(999, 9999, "Pathogenic", pos="-1", ref="na", alt="na"),         # unknown, unplaced
           vs(9998, -1, "Pathogenic", pos="-1", ref="na", alt="na")]          # unknown, no AlleleID
EHIST_ROWS = [vs(3, 13, "Likely pathogenic", pos=1003)]
EDEV_LINES = ["3\t1999\t8\tA\tT\t.\t.\tALLELEID=18;CLNSIG=Uncertain_significance;"
              "CLNREVSTAT=criteria_provided,_single_submitter",
              "3\t1019\t19\tA\t.\t.\t.\tALLELEID=29;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter",
              "3\t1024\t9024\tA\t.\t.\t.\tALLELEID=99024;CLNSIG=Pathogenic"]


def write_inputs(d):
    write_vs(d / "t1.gz", T1_ROWS)
    write_vs(d / "t0.gz", T0_ROWS)
    write_vs(d / "ehist.gz", EHIST_ROWS)
    write_vcf(d / "edev.gz", EDEV_LINES)
    (d / "pos.tsv").write_text("chrom\tstart\tend\tsource\nchr3\t1016\t1016\tbenchmark X\n", encoding="utf-8")
    (d / "rep.tsv").write_text("replaced_variation_id\tcurrent_variation_id\n7022\t22\n", encoding="utf-8")
    (d / "mane.tsv").write_text("refseq_transcript\tgene_id\tplus_clinical\nNM_000001.1\t101\tyes\n",
                                encoding="utf-8")
    (d / "submitters.txt").write_text("organization_id\tsubmitter_name\n506627\t"
                                      "Broad Center for Mendelian Genomics, Broad Institute of MIT and Harvard\n",
                                      encoding="utf-8")
    with gzip.open(d / "subs.gz", "wt", encoding="utf-8", newline="") as h:
        h.write("##Overview\n#VariationID\tClinicalSignificance\tSubmitter\tSCV\n"
                "17\tPathogenic\tBroad Center for Mendelian Genomics, Broad Institute of MIT and Harvard\tSCV1\n"
                "1\tPathogenic\tOther Lab\tSCV2\n")


def derived(d):
    return {"positions": d / "pos.tsv", "replacements": d / "rep.tsv", "mane": d / "mane.tsv"}


def exposures(d):
    return {"E-dev": d / "edev.gz", "T0": d / "t0.gz", "E-hist": d / "ehist.gz"}


def read_rows(path):
    with gzip.open(path, "rt", encoding="utf-8") as h:
        header = h.readline().rstrip("\n").split("\t")
        return {r[0]: dict(zip(header, r)) for r in (line.rstrip("\n").split("\t") for line in h)}


class Build(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        d = cls.d = Path(cls.tmp.name)
        write_inputs(d)
        files = dict(vsets.FILES, positions=d / "pos.tsv", submitters=d / "submitters.txt",
                     replacements=d / "rep.tsv", mane=d / "mane.tsv", chain=d / "pos.tsv")
        with mock.patch.dict(vsets.FILES, files), mock.patch.object(vsets.liftover, "load", lambda p: None):
            cls.summary = vsets.build(d / "t1.gz", exposures(d), d / "subs.gz", d / "out")
        cls.rows = read_rows(d / "out" / "clinvar_variants.tsv.gz")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_one_row_per_candidate(self):
        # not candidates: 4 (ambiguous duplicate), 5 (no GRCh38 row), 9 and 23 (not P/LP or B/LB)
        self.assertEqual(set(self.rows), {"1", "2", "3", "6", "7", "8", "10", "11", "12", "13", "14", "15",
                                          "16", "17", "18", "19", "20", "21", "22", "24", "25", "26"})

    def test_new_label_set(self):
        new = {v for v, r in self.rows.items() if r["new_label"] == "yes"}
        self.assertEqual(new, {"1", "6", "7", "10", "11", "12", "13", "14", "15", "25", "26"})

    def test_exclusion_reasons(self):
        r = self.rows
        for vid in ("2", "3", "18", "19", "20", "22"):
            self.assertEqual(r[vid]["exclusion"], "labelled before", vid)
        self.assertEqual(r["8"]["exclusion"], "identity collision with E-dev")
        self.assertEqual(r["16"]["development_material"], "benchmark X")
        self.assertEqual(r["17"]["development_material"], "development submitters")
        for vid in ("21", "24"):
            self.assertEqual(r[vid]["exclusion"], "unresolved historical identity (same position)", vid)
        self.assertEqual(r["2"]["labelled_before_in"], "T0")

    def test_regions_and_fields(self):
        r = self.rows
        self.assertEqual((r["1"]["region"], r["1"]["family"], r["1"]["intron_offset"]), ("deep_intronic", "intronic", "60"))
        self.assertEqual(r["14"]["intron_offset"], "7")
        self.assertEqual(r["1"]["mane_select"], "yes")
        self.assertEqual(r["25"]["mane_select"], "no")
        self.assertEqual(r["26"]["mane_select"], "no")
        self.assertEqual(r["11"]["variant_type"], "small_indel")
        self.assertEqual(r["6"]["chromosome_group"], "excluded")
        self.assertEqual((r["1"]["chrom"], r["1"]["pos"], r["1"]["label"]), ("3", "1001", "P/LP"))
        self.assertEqual(r["1"]["inside_long_unknown_record"], "no")

    def test_summary(self):
        s = self.summary
        self.assertEqual(s["development_material"], {"benchmark X": 1, "development submitters": 1})
        self.assertEqual(s["unknown_record_blocking"]["candidates_blocked"], 2)
        n = sum(c["n"] for c in s["counts"] if c["set"] == "new labels")
        self.assertEqual(n, 11)
        self.assertEqual(sum(c["n"] for c in s["counts"] if c["set"] == "all labels"), 22)


class UnknownRecordSpan(unittest.TestCase):
    def test_only_small_unknown_records_block(self):
        rec = lambda vid, pos: cr.Record("T1", vid, vid, "GRCh38", "3", f"3:{pos}:A:G", "n", "Pathogenic",  # noqa: E731
                                         "r", 2, True, "1")
        cands = {"1": idx.Candidate(rec("1", 100), "P/LP"), "2": idx.Candidate(rec("2", 5000), "P/LP")}
        index = idx.ExposureIndex(cands)
        unknown = {"vid:a": {"pos38": {("3", 100, 100)}, "pos37": set(), "sources": {"T0"}},
                   "vid:b": {"pos38": {("3", 4000, 9000)}, "pos37": set(), "sources": {"T0"}}}
        with mock.patch.object(index, "unmatched", lambda: unknown):
            out = index.block_unmatched()
        self.assertEqual(cands["1"].unresolved, "unresolved historical identity (same position)")
        self.assertIsNone(cands["2"].unresolved)                       # inside a 5,001-base record only
        self.assertTrue(cands["2"].near_long_unknown)                  # flagged for the sensitivity analysis
        self.assertFalse(cands["1"].near_long_unknown)
        self.assertEqual(out["records_with_a_span_over_limit"], 1)
        self.assertEqual(out["candidates_inside_a_longer_record"], 1)


if __name__ == "__main__":
    unittest.main()
