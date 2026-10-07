import gzip
import io
import math
import os
import struct
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

os.environ["AVI_EVAL_SYNTHETIC"] = "1"
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import part2_data as p2  # noqa: E402

# A toy genome: chr1 forward sequence, 1-based positions.
GENOME = {"1": "ACGTTTTGCAAGGCATCGATCGGA"}


def fetch(chrom, start, end):
    seq = GENOME.get(str(chrom), "")
    if start < 1 or end > len(seq) or end < start:
        return ""
    return seq[start - 1:end]


# ---------------------------------------------------------------------------------------------------
# Spreadsheet readers
# ---------------------------------------------------------------------------------------------------
def make_xlsx() -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("xl/workbook.xml",
                   '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
                   'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
                   '<sheets><sheet name="Data" sheetId="1" r:id="rId1"/></sheets></workbook>')
        z.writestr("xl/_rels/workbook.xml.rels",
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                   '<Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
        z.writestr("xl/sharedStrings.xml",
                   '<sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
                   '<si><t>name</t></si><si><r><t>val</t></r><r><t>ue</t></r></si><si><t>x</t></si></sst>')
        z.writestr("xl/worksheets/sheet1.xml",
                   '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>'
                   '<row r="1"><c r="A1" t="s"><v>0</v></c><c r="C1" t="s"><v>1</v></c></row>'
                   '<row r="2"><c r="A2" t="s"><v>2</v></c><c r="B2" t="b"><v>1</v></c><c r="C2"><v>2.5</v></c></row>'
                   '<row r="3"/></sheetData></worksheet>')
    return buf.getvalue()


def rec(rid: int, body: bytes) -> bytes:
    return struct.pack("<HH", rid, len(body)) + body


def biff_string(s: str, len_bytes: int = 2) -> bytes:
    return struct.pack("<H" if len_bytes == 2 else "<B", len(s)) + b"\x00" + s.encode("latin-1")


def make_biff() -> bytes:
    """Workbook globals (BOUNDSHEET, SST split by CONTINUE) and one sheet with every cell record."""
    sheet_name = biff_string("S1", 1)
    sst_head = struct.pack("<II", 3, 3) + biff_string("alpha") + struct.pack("<H", 6) + b"\x00" + b"bet"
    cont = b"\x01" + "ahZ".encode("utf-16-le")             # rest of 'betahZ' in 16-bit characters
    third = biff_string("gamma")
    glob = rec(0x809, b"\x00" * 16)
    bound_placeholder = rec(0x85, struct.pack("<I", 0) + b"\x00\x00" + sheet_name)
    glob_rest = rec(0xFC, sst_head) + rec(0x3C, cont + third) + rec(0x0A, b"")
    sheet_offset = len(glob) + len(bound_placeholder) + len(glob_rest)
    bound = rec(0x85, struct.pack("<I", sheet_offset) + b"\x00\x00" + sheet_name)
    rk_int = (7 << 2) | 2                                  # integer 7
    rk_div = (1234 << 2) | 2 | 1                           # 12.34
    sheet = (rec(0x809, b"\x00" * 16)
             + rec(0xFD, struct.pack("<HHHI", 0, 0, 0, 1))           # A1 = 'betahZ'
             + rec(0x203, struct.pack("<HHHd", 0, 1, 0, 3.5))        # B1 = 3.5
             + rec(0x27E, struct.pack("<HHHI", 1, 0, 0, rk_int))     # A2 = 7
             + rec(0xBD, struct.pack("<HHHIHIH", 2, 0, 0, rk_div, 0, rk_int, 1))  # A3, B3
             + rec(0x204, struct.pack("<HHH", 3, 0, 0) + biff_string("lbl"))      # A4
             + rec(0x205, struct.pack("<HHHBB", 3, 1, 0, 1, 0))      # B4 = True
             + rec(0x06, struct.pack("<HHH", 4, 0, 0) + b"\x00\x00\x00\x00\x00\x00\xff\xff" + b"\x00" * 6)
             + rec(0x207, biff_string("fs"))                         # A5 = 'fs' (formula string)
             + rec(0x06, struct.pack("<HHH", 4, 1, 0) + struct.pack("<d", 2.0) + b"\x00" * 6)  # B5 = 2.0
             + rec(0x0A, b""))
    return glob + bound + glob_rest + sheet


def make_cfb(stream: bytes, name: str = "Workbook") -> bytes:
    """Minimal compound file: sector 0 = FAT, sector 1 = directory, sectors 2.. = the stream."""
    sec = 512
    stream = stream + b"\x00" * (-len(stream) % sec)
    n_stream = len(stream) // sec
    end, fatsect, free = 0xFFFFFFFE, 0xFFFFFFFD, 0xFFFFFFFF
    fat = [fatsect, end] + [i + 1 for i in range(2, 2 + n_stream - 1)] + [end]
    fat += [free] * (sec // 4 - len(fat))
    header = bytearray(512)
    header[0:8] = bytes.fromhex("D0CF11E0A1B11AE1")
    struct.pack_into("<HH", header, 30, 9, 6)
    struct.pack_into("<II", header, 44, 1, 1)
    struct.pack_into("<IIIII", header, 56, 4096, end, 0, end, 0)
    struct.pack_into("<109I", header, 76, 0, *([free] * 108))

    def entry(nm: str, typ: int, start: int, size: int) -> bytes:
        e = bytearray(128)
        raw = (nm + "\x00").encode("utf-16-le")
        e[0:len(raw)] = raw
        struct.pack_into("<H", e, 64, len(raw))
        e[66] = typ
        struct.pack_into("<III", e, 68, free, free, free)
        struct.pack_into("<II", e, 116, start, size)
        return bytes(e)

    directory = entry("Root Entry", 5, end, 0) + entry(name, 2, 2, len(stream)) + b"\x00" * 256
    return bytes(header) + struct.pack(f"<{sec // 4}I", *fat) + directory + stream


class Readers(unittest.TestCase):
    def test_xlsx(self):
        sheets = p2.read_xlsx(make_xlsx())
        self.assertEqual(sheets["Data"], [["name", None, "value"], ["x", True, 2.5]])

    def test_biff_cells_and_continued_shared_string(self):
        rows = p2.parse_biff(make_biff())["S1"]
        self.assertEqual(rows[0][:2], ["betahZ", 3.5])
        self.assertEqual(rows[1][0], 7.0)
        self.assertAlmostEqual(rows[2][0], 12.34)
        self.assertEqual(rows[2][1], 7.0)
        self.assertEqual(rows[3][:2], ["lbl", True])
        self.assertEqual(rows[4][:2], ["fs", 2.0])

    def test_xls_container(self):
        stream = make_biff()
        stream = stream + b"\x00" * (4200 - len(stream))      # at least the mini-stream cutoff
        sheets = p2.read_xls(make_cfb(stream))
        self.assertEqual(sheets["S1"][0][0], "betahZ")

    def test_sheet_frame_keeps_only_named_columns(self):
        rows = [["title", None, None], ["a", "Score X", "b "], [1, 9, 2], [None, None, None]]
        df = p2.sheet_frame(rows, 1, ["a", "b"])
        self.assertEqual(list(df.columns), ["a", "b"])
        self.assertEqual(df.values.tolist(), [[1, 2]])
        self.assertEqual(list(p2.sheet_frame(rows, 1).columns), ["a", "Score X", "b"])


# ---------------------------------------------------------------------------------------------------
# Genome access and alleles
# ---------------------------------------------------------------------------------------------------
class Genome(unittest.TestCase):
    def test_fai_and_fetch(self):
        with tempfile.TemporaryDirectory() as d:
            fa = Path(d) / "g.fa"
            fa.write_bytes(b">chr1 desc\nACGT\nTTTG\nCA\n>chr2\nGGGG\nCC\n")
            f = p2.Fasta(fa)
            self.assertEqual(f.index["chr1"], (10, 11, 4, 5))
            self.assertEqual(f.fetch("1", 3, 6), "GTTT")
            self.assertEqual(f.fetch("chr2", 4, 6), "GCC")
            self.assertEqual(f.fetch("2", 6, 7), "")
            self.assertEqual(f.fetch("9", 1, 1), "")
            f.handle.close()

    def test_normalize(self):
        # deletion of one T inside TTTT (positions 4-7) written at the right end -> left-aligned at G3
        self.assertEqual(p2.normalize(fetch, "1", 6, "TT", "T"), (3, "GT", "G"))
        # reverse-complement form with the shared base at the end (AG -> G) gets an anchor on the left
        self.assertEqual(p2.normalize(fetch, "1", 10, "AG", "G"), (9, "CA", "C"))
        self.assertEqual(p2.normalize(fetch, "1", 2, "C", "T"), (2, "C", "T"))
        self.assertEqual(p2.normalize(fetch, "1", 4, "T", "TT"), (3, "G", "GT"))

    def test_variant_type_and_revcomp(self):
        self.assertEqual(p2.variant_type("A", "G"), "SNV")
        self.assertEqual(p2.variant_type("AT", "A"), "deletion")
        self.assertEqual(p2.variant_type("A", "AT"), "insertion")
        self.assertEqual(p2.variant_type("AT", "GC"), "delins")
        self.assertEqual(p2.revcomp("AACG"), "CGTT")

    def test_bh(self):
        q = p2.bh(np.array([0.01, 0.04, 0.03, 0.2]))
        np.testing.assert_allclose(q, [0.04, 0.16 / 3, 0.16 / 3, 0.2])


# ---------------------------------------------------------------------------------------------------
# Datasets
# ---------------------------------------------------------------------------------------------------
def ron_rep(values: dict, mean=50.0, sd=1.0) -> pd.DataFrame:
    iso = p2.RON_ISOFORMS
    rows = {"wt mean": [mean] * 6, "wt SD": [sd] * 6}
    rows.update({k: v for k, v in values.items()})
    return pd.DataFrame.from_dict(rows, orient="index", columns=iso)


class Datasets(unittest.TestCase):
    def test_pou1f1_complements_transcript_alleles(self):
        df = pd.DataFrame([{"Chromosome": 3.0, "Position (HG38)": 100.0, "Position (HGVS)": "c.1-5", "Ref": "G",
                            "Alt": "A", "Variant Consequence": "Intronic", "Splice Disruptive": False,
                            "Intermediate": True, "Beta fold-change": 2.5, "Skip fold-change": 1.0,
                            "Other fold-change": 0.5}])
        r = p2.harmonise_pou1f1(df)[0]
        self.assertEqual((r["chrom"], r["pos"], r["ref"], r["alt"]), ("3", 100, "C", "T"))
        self.assertEqual((r["splice_altering"], r["effect_size"]), ("intermediate", 2.5))

    def test_wt1_hgvs_and_liftover(self):
        df = pd.DataFrame([{"HGVS Variant": "c.1447+40T>A", "HG19 Position": 500.0, "Variant Consequence": "intronic",
                            "Log2 KTS+/KTS- Ratio": -1.2, "Splice Disruptive": True, "Intermediate": False},
                           {"HGVS Variant": "c.1355-3C>G", "HG19 Position": 9.0, "Variant Consequence": "intronic",
                            "Log2 KTS+/KTS- Ratio": 0.1, "Splice Disruptive": False, "Intermediate": False}])
        lift = lambda span: [("11", span[1] + 1000, span[2] + 1000)] if span[1] == 500 else []  # noqa: E731
        rows, unlifted = p2.harmonise_wt1(df, lift)
        self.assertEqual(unlifted, 1)
        self.assertEqual((rows[0]["pos"], rows[0]["ref"], rows[0]["alt"], rows[0]["splice_altering"]),
                         (1500, "A", "T", "yes"))

    def test_ron_calls(self):
        big = [50, 50, 50, 50, 50, 70]          # +20 points, z = 20 in every replicate
        small = [50, 50, 50, 50, 50, 53]        # +3 points: significant z but below 5 points
        reps = [ron_rep({"M1": big, "M2": small, "M3": [50] * 6}) for _ in range(3)]
        calls = p2.ron_calls(reps)
        self.assertTrue(calls.loc["M1", "splice_altering"])
        self.assertFalse(calls.loc["M2", "splice_altering"])
        self.assertFalse(calls.loc["M3", "splice_altering"])
        self.assertAlmostEqual(calls.loc["M1", "max_abs_change"], 20.0)

    def test_ron_minus_strand_and_indel_anchor(self):
        const = 1000
        reps = [ron_rep({"A26T": [50, 50, 50, 50, 50, 70], "CT38C": [50] * 6}) for _ in range(3)]
        average = pd.DataFrame({"mutation type": ["SNV", "DEL"],
                                "genomic position": [f"chr3:{const - 26}", f"chr3:{const - 37}"]},
                               index=["A26T", "CT38C"])
        median = pd.DataFrame({"Mutation": ["G40A"], "Mutation type": ["SNV"], "Genomic position": [f"chr3:{const - 40}"],
                               "Median delta AE inclusion (%)": [-8.0], "Median delta other (%)": [2.0]})
        lift = lambda span: [("3", span[1] + 5, span[2] + 5)]  # noqa: E731
        rows, info = p2.harmonise_ron(reps, average, median, lift)
        self.assertEqual(info["hg19_offset"], const)
        by = {r["variant_id_in_source"]: r for r in rows}
        self.assertEqual((by["A26T"]["pos"], by["A26T"]["ref"], by["A26T"]["alt"]), (const - 26 + 5, "T", "A"))
        self.assertEqual(by["A26T"]["splice_altering"], "yes")
        # anchor C at minigene 37 (from the genomic position), deleted T at 38: forward AG -> G at C-38
        self.assertEqual((by["CT38C"]["pos"], by["CT38C"]["ref"], by["CT38C"]["alt"]), (const - 38 + 5, "AG", "G"))
        self.assertEqual((by["G40A"]["subset"], by["G40A"]["splice_altering"]), ("median_estimate", "yes"))

    def test_cd19(self):
        validated = pd.DataFrame({"coordinate (hg38)": [None, "chr16:1465"], "corresponding minigene mutation":
                                  [None, "T465A"]})
        off = p2.cd19_offset(validated)
        self.assertEqual(off, 1000)
        effects = pd.DataFrame({"mutation": ["T465A"] * 2 + ["CT10C"] * 2, "type": ["SNV"] * 2 + ["INDEL"] * 2,
                                "region": ["intron 1"] * 2 + ["upstream region"] * 2,
                                "average delta frequency": [0.1, -0.3, 0.01, 0.0],
                                "splicing-effective mutation": ["No", "Yes", "No", "No"], "isoform": ["a", "b"] * 2})
        rows = {r["variant_id_in_source"]: r for r in p2.harmonise_cd19(effects, off)}
        self.assertEqual((rows["T465A"]["pos"], rows["T465A"]["splice_altering"]), (1465, "yes"))
        self.assertAlmostEqual(rows["T465A"]["effect_size"], 0.3)
        self.assertEqual((rows["CT10C"]["pos"], rows["CT10C"]["ref"], rows["CT10C"]["alt"]), (1010, "CT", "C"))

    def test_compass_subsets_and_effect(self):
        def line(ref_id, source, var, dl, refpsi=(0.5,) * 5, seq_type="single"):
            r = {"Reference": ref_id, "snp": "x", "source": source, "seq_type": seq_type, "variant_hg38": var,
                 "gene_exon": "GENE1 exon 2", "transcript_id": "ENST1", "exon_start_hg38": 100, "exon_end_hg38": 150}
            for c, v, w in zip(p2.COMPASS_CELLS, dl, refpsi):
                r[f"{c}_delta_logit_pooled"] = v
                r[f"{c}_wt_pooled_psi_raw"] = w
            return r

        nan = float("nan")
        df = pd.DataFrame([
            line(1, "clinvar_single", "chr1:120:A>G", [1.5, 0.7, nan, nan, nan]),           # random, mean 1.1
            line(2, "clinvar_single", "chr1:121:C>T", [0.2, 0.1, 0, 0, 0]),                 # in ClinVar
            line(3, "exac", "chr1:90:G>A", [-0.5, nan, nan, nan, nan], refpsi=(0.0,) * 5),  # reference PSI 0
            line(4, "geuvadis", "chr1:122:G>A", [0.0] * 5),
            line(5, "clinvar_single", "chr1:123:G>A", [nan] * 5),                           # not measured
            line(6, "clinvar_double", "chr1:1:A>G;chr1:2:A>G", [3] * 5, seq_type="double"),
            {**line(7, "clinvar_single", None, [nan] * 5, seq_type="none"), "snp": "none"},
        ])
        rows, info = p2.harmonise_compass(df, {("1", 121, "C", "T")})
        by = {r["variant_id_in_source"].split(":")[0]: r for r in rows}
        self.assertEqual((info["doubles"], info["references"], info["single_not_measured"]), (1, 1, 1))
        self.assertEqual((by["1"]["subset"], by["1"]["use"], by["1"]["splice_altering"]), ("random_snv", "main", "yes"))
        self.assertAlmostEqual(by["1"]["effect_size"], 1.1)
        self.assertEqual((by["2"]["subset"], by["2"]["use"]), ("clinvar", "not_used"))
        self.assertEqual((by["3"]["subset"], by["3"]["use"]), ("exac;reference_psi_0_or_1", "not_used"))
        self.assertTrue(by["3"]["region_note"].startswith("intron"))
        self.assertEqual(by["4"]["use"], "not_used")

    def test_opensplice(self):
        df = pd.DataFrame([
            {"CHROM": "chr1", "POS": "5", "REF": "T", "ALT": "C", "exon_id": "G1_e2", "variant_id": "v1", "lib_id": "L",
             "mut_type": "sub", "region": "Exon", "psi": 40.0, "delta_logit": -1.0, "significant": "yes"},
            {"CHROM": "chr1", "POS": "5_6del", "REF": "TTT", "ALT": "T", "exon_id": "G1_e2", "variant_id": "v2",
             "lib_id": "L", "mut_type": "d2", "region": "Exon", "psi": 40.0, "delta_logit": 0.1, "significant": "no"},
            {"CHROM": "chr1", "POS": "9", "REF": "C", "ALT": "A", "exon_id": "G1_e2", "variant_id": "v3", "lib_id": "L",
             "mut_type": "sub", "region": "Exon", "psi": float("nan"), "delta_logit": float("nan"), "significant": "no"},
            {"CHROM": None, "POS": None, "REF": None, "ALT": None, "exon_id": "G1_e2", "variant_id": "wt",
             "lib_id": "L", "mut_type": "wt", "region": "Exon", "psi": 50.0, "delta_logit": 0.0, "significant": None},
        ])
        rows, info = p2.harmonise_opensplice(df)
        self.assertEqual(info, {"wt_rows": 1, "not_measured": 1})
        self.assertEqual([(r["pos"], r["ref"], r["alt"], r["splice_altering"], r["gene"]) for r in rows],
                         [(5, "T", "C", "yes", "G1"), (4, "TTT", "T", "no", "G1")])

    def test_mfass(self):
        df = pd.DataFrame({"id": ["e_000", "e_001", "e_002"], "ensembl_id": ["E1"] * 3, "chr": ["chr1"] * 3,
                           "ref_allele": [None, "C", "G"], "alt_allele": [None, "T", "A"],
                           "snp_position_hg38_1based": [None, 2.0, 3.0], "label": ["exon"] * 3,
                           "category": ["natural", "mutant", "mutant"], "v2_dpsi": [0.0, -0.7, -0.1],
                           "strong_lof": [False, True, None]})
        rows, info = p2.harmonise_mfass(df)
        self.assertEqual(info, {"mutant_rows": 2, "no_call": 1})
        self.assertEqual((rows[0]["use"], rows[0]["splice_altering"], rows[0]["pos"]), ("comparison", "yes", 2))


class Assembly(unittest.TestCase):
    def rows(self):
        base = dict(effect_size=1.0, effect_measure="m", splice_altering="no", threshold_rule="t", gene="G",
                    region_note="", subset="s", use="main")
        return [
            {**base, "dataset": "d", "variant_id_in_source": "b", "chrom": "1", "pos": 2, "ref": "C", "alt": "T"},
            {**base, "dataset": "d", "variant_id_in_source": "a", "chrom": "1", "pos": 2, "ref": "C", "alt": "T"},
            {**base, "dataset": "d", "variant_id_in_source": "c", "chrom": "1", "pos": 3, "ref": "A", "alt": "T"},
            {**base, "dataset": "d", "variant_id_in_source": "d", "chrom": "1", "pos": 6, "ref": "TT", "alt": "T"},
            {**base, "dataset": "e", "variant_id_in_source": "a", "chrom": "1", "pos": 2, "ref": "C", "alt": "T"},
        ]

    def test_finalize(self):
        df, checks = p2.finalize(self.rows(), fetch)
        self.assertEqual(checks["ref_mismatch"], {"d": 1})                  # pos 3 is G, not A
        self.assertEqual(checks["duplicates"], {"d": 1})
        d = df[df["dataset"] == "d"].set_index("variant_id_in_source")
        self.assertEqual(d.loc["a", "use"], "main")
        self.assertEqual((d.loc["b", "use"], d.loc["b", "subset"]), ("not_used", "s;duplicate"))
        self.assertEqual((d.loc["d", "pos"], d.loc["d", "ref"], d.loc["d", "alt"], d.loc["d", "variant_type"]),
                         (3, "GT", "G", "deletion"))
        self.assertEqual(df[df["dataset"] == "e"]["use"].tolist(), ["main"])   # other dataset untouched

    def test_write_and_clinvar_keys(self):
        df, _ = p2.finalize(self.rows(), fetch)
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "x" / "t.tsv.gz"
            p2.write(df, out)
            back = pd.read_csv(out, sep="\t", dtype=str)
            self.assertEqual(list(back.columns), p2.COLUMNS)
            self.assertEqual(len(back), 4)
            vcf = Path(d) / "c.vcf.gz"
            with gzip.open(vcf, "wt") as h:
                h.write("##x\n#CHROM\tPOS\tID\tREF\tALT\n1\t5\t.\tA\tG,T\t.\n1\t6\t.\tC\tT\t.\n")
            self.assertEqual(p2.clinvar_keys(vcf, {("1", 5)}), {("1", 5, "A", "G"), ("1", 5, "A", "T")})


if __name__ == "__main__":
    unittest.main()
