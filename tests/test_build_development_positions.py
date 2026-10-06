"""Tests for the development positions (source table sections 3 and 4), on synthetic files."""
import contextlib
import gzip
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
os.environ["AVI_EVAL_SYNTHETIC"] = "1"

import build_development_positions as bp  # noqa: E402
import clinvar_index as idx  # noqa: E402

GFF = "\n".join([
    "##gff-version 3",
    "NC_000017.11\tRefSeq\tregion\t1\t83257441\t.\t+\t.\tID=NC_000017.11:1..83257441;chromosome=17;genome=chromosome",
    "NC_000023.11\tRefSeq\tregion\t1\t156040895\t.\t+\t.\tID=NC_000023.11:1..156040895;chromosome=X;genome=chromosome",
    "NC_000024.10\tRefSeq\tregion\t1\t57227415\t.\t+\t.\tID=NC_000024.10:1..57227415;chromosome=Y;genome=chromosome",
    "NT_187361.1\tRefSeq\tregion\t1\t175055\t.\t+\t.\tID=NT_187361.1:1..175055;genome=genomic",
    "NC_000017.11\tBestRefSeq\tgene\t43044295\t43170327\t.\t-\t.\tID=gene-BRCA1;Dbxref=GeneID:672,HGNC:HGNC:1100;Name=BRCA1;gene=BRCA1",
    "NT_187361.1\tBestRefSeq\tgene\t100\t200\t.\t-\t.\tID=gene-BRCA1-2;Dbxref=GeneID:672;Name=BRCA1;gene=BRCA1",
    "NC_000023.11\tBestRefSeq\tgene\t1000\t2000\t.\t+\t.\tID=gene-PAR1;Dbxref=GeneID:999;Name=PARG;gene=PARG",
    "NC_000024.10\tBestRefSeq\tgene\t5000\t6000\t.\t+\t.\tID=gene-PAR1-2;Dbxref=GeneID:999;Name=PARG;gene=PARG",
    "NC_000017.11\tBestRefSeq\tgene\t10\t20\t.\t+\t.\tID=gene-TWIN;Dbxref=GeneID:1;Name=TWIN;gene=TWIN",
    "NC_000017.11\tBestRefSeq\tgene\t30\t40\t.\t+\t.\tID=gene-TWIN-2;Dbxref=GeneID:2;Name=TWIN;gene=TWIN",
]) + "\n"


def write_gff(path):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as h:
        h.write(GFF)


class GeneRegions(unittest.TestCase):
    def test_primary_chromosomes_only(self):
        with tempfile.TemporaryDirectory() as t:
            write_gff(Path(t) / "g.gff.gz")
            chroms = bp.chromosome_map(Path(t) / "g.gff.gz")
            self.assertEqual(chroms, {"NC_000017.11": "17", "NC_000023.11": "X", "NC_000024.10": "Y"})

    def test_region_per_chromosome_padded(self):
        with tempfile.TemporaryDirectory() as t:
            write_gff(Path(t) / "g.gff.gz")
            regions = bp.gene_regions(Path(t) / "g.gff.gz", ["BRCA1", "PARG"])
            self.assertEqual(regions["BRCA1"], {"gene_id": "672", "intervals": [("17", 43043295, 43171327)]})
            # A gene on X and Y: one span per chromosome, never across chromosomes.
            self.assertEqual(regions["PARG"]["intervals"], [("X", 1, 3000), ("Y", 4000, 7000)])

    def test_a_symbol_must_name_one_gene(self):
        with tempfile.TemporaryDirectory() as t:
            write_gff(Path(t) / "g.gff.gz")
            for symbol in ("TWIN", "ABSENT"):
                with self.assertRaises(idx.InputFormatError, msg=symbol):
                    bp.gene_regions(Path(t) / "g.gff.gz", [symbol])


class Benchmarks(unittest.TestCase):
    def test_traitgym_positions_of_every_row(self):
        import pyarrow as pa
        import pyarrow.parquet as pq
        with tempfile.TemporaryDirectory() as t:
            p = Path(t) / "test.parquet"
            pq.write_table(pa.table({"chrom": ["1", "X", "2"], "pos": [100, 200, 300], "ref": ["A", "CTG", "T"],
                                     "alt": ["G", "C", "TAA"], "label": [True, False, True]}), p)
            # SNV; deletion (anchor and deleted bases); VCF-style insertion (anchor and the next base).
            self.assertEqual(bp.traitgym_intervals(p), [("1", 100, 100), ("X", 200, 202), ("2", 300, 301)])

    def test_ncvardb_rules(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            (d / "p.tsv").write_text("ID\tgenome\tchr\tpos\tref\talt\tmut_type\n"
                                     "a\thg38\t1\t100\tA\tG\tsubstitution\n"
                                     "b\thg38\t2\t200\tAAGT\t-\tdeletion\n"
                                     "c\thg38\t3\t300\t-\tT\tinsertion\n"
                                     "i\thg38\t6\t600\tA\tAGG\tinsertion\n"
                                     "x\thg38\t5\t500\tCC\tA\tdeletion\n", encoding="utf-8")
            (d / "p.vcf").write_text("##fileformat=VCFv4.3\n#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                                     "5\t500\trs1\tCC\tA\t50\tPASS\tNC=x\n"
                                     "6\t600\trs2\tA\tAGG\t50\tPASS\tNC=i\n", encoding="utf-8")
            # '-' insertion: p-1..p+1 (right under either reading of pos); VCF-style insertion: anchor and next base.
            self.assertEqual(bp.ncvardb_intervals([d / "p.tsv"], [d / "p.vcf"]),
                             [("1", 100, 100), ("2", 200, 203), ("3", 299, 301), ("5", 500, 501), ("6", 600, 601)])
            # a VCF row is accepted only at the hg38 position of its ncVarDB row.
            bad = {"declared GRCh37": "##reference=GRCh37\n#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                                      "5\t500\trs1\tCC\tA\t50\tPASS\tNC=x\n",
                   "other position": "#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n5\t9500\trs1\tCC\tA\t50\tPASS\tNC=x\n",
                   "unknown identifier": "#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n5\t500\trs1\tCC\tA\t50\tPASS\tNC=zz\n",
                   "no identifier": "#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n5\t500\trs1\tCC\tA\t50\tPASS\t.\n"}
            for label, text in bad.items():
                with self.subTest(label=label):
                    (d / "b.vcf").write_text(text, encoding="utf-8")
                    with self.assertRaises(idx.InputFormatError):
                        bp.ncvardb_intervals([d / "p.tsv"], [d / "b.vcf"])

    def test_comma_separated_identifiers_on_both_sides(self):
        # identifiers are indexed one by one; conflicting assignments are refused.
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            (d / "p.tsv").write_text("ID\tgenome\tchr\tpos\tref\talt\na,b\thg38\t3\t100\tA\tG\n", encoding="utf-8")
            (d / "p.vcf").write_text("#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n3\t100\t.\tA\tG\t50\tPASS\tNC=a,b\n",
                                     encoding="utf-8")
            self.assertEqual(bp.ncvardb_intervals([d / "p.tsv"], [d / "p.vcf"]), [("3", 100, 100)])
            (d / "c.tsv").write_text("ID\tgenome\tchr\tpos\tref\talt\na\thg38\t3\t100\tA\tG\na\thg38\t4\t100\tA\tG\n",
                                     encoding="utf-8")
            with self.assertRaises(idx.InputFormatError):
                bp.ncvardb_intervals([d / "c.tsv"], [])

    def test_file_inventories(self):
        # specified here independently from the module.
        self.assertEqual(sorted(f for files in bp.TRAITGYM.values() for f in files),
                         sorted(["mendelian_traits_all", "mendelian_traits_matched_9", "mendelian_traits_v21_matched_9",
                                 "mendelian_traits_v22_matched_9", "complex_traits_all", "complex_traits_matched_9",
                                 "complex_traits_v22_matched_9"]))
        self.assertTrue(all(f.startswith("mendelian") for f in bp.TRAITGYM["TraitGym Mendelian"]))
        self.assertTrue(all(f.startswith("complex") for f in bp.TRAITGYM["TraitGym Complex"]))
        self.assertEqual(sorted(bp.NCVAR_TSV + bp.NCVAR_VCF), ["ncVar_benign.tsv", "ncVar_benign.vcf",
                                                              "ncVar_pathogenic.tsv", "ncVar_pathogenic.vcf"])

    def test_ncvardb_rows_must_be_hg38(self):
        with tempfile.TemporaryDirectory() as t:
            d = Path(t)
            (d / "q.tsv").write_text("ID\tgenome\tchr\tpos\tref\talt\tmut_type\na\thg19\t1\t100\tA\tG\ts\n",
                                     encoding="utf-8")
            with self.assertRaises(idx.InputFormatError):
                bp.ncvardb_intervals([d / "q.tsv"], [])


class Assemble(unittest.TestCase):
    def test_merged_intervals_by_source(self):
        rows = bp.merged_rows({"saturation genome editing": [("17", 10, 20), ("17", 15, 30), ("17", 32, 40)],
                               "ncVarDB": [("chr1", 5, 5)]})
        self.assertEqual(rows, [("17", 10, 30, "saturation genome editing"), ("17", 32, 40, "saturation genome editing"),
                                ("1", 5, 5, "ncVarDB")])

    def test_written_file_is_read_by_the_counting_reader(self):
        with tempfile.TemporaryDirectory() as t:
            out = Path(t) / "positions.tsv"
            bp.write_positions([("17", 10, 30, "saturation genome editing"), ("1", 5, 5, "ncVarDB")], out)
            with contextlib.nullcontext():
                got = idx.read_positions(out)
            self.assertEqual(dict(got), {"saturation genome editing": [("17", 10, 30)], "ncVarDB": [("1", 5, 5)]})

    def test_missing_file_keeps_the_observed_diagnostics(self):
        # a missing file selects the gene region but the parsed score sets stay reported.
        region = {"gene_id": "7428", "intervals": [("3", 1, 9)]}
        parsed = {"scored": 7, "unpositioned": 5, "intervals": [("3", 5, 5)]}
        self.assertEqual(bp.sge_choice([parsed], region, missing=["mavedb/x.mapped.json"]),
                         ("gene region", [("3", 1, 9)]))
        self.assertEqual(bp.sge_choice([dict(parsed, unpositioned=0)], region, missing=[]),
                         ("positions", [("3", 5, 5)]))

    def test_sge_study_uses_gene_region_when_any_variant_is_unpositioned(self):
        region = {"gene_id": "7428", "intervals": [("3", 1, 9)]}
        positioned = {"scored": 2, "unpositioned": 0, "intervals": [("3", 5, 5)]}
        partly = {"scored": 2, "unpositioned": 1, "intervals": [("3", 5, 5)]}
        self.assertEqual(bp.sge_choice([positioned], region), ("positions", [("3", 5, 5)]))
        self.assertEqual(bp.sge_choice([positioned, partly], region), ("gene region", [("3", 1, 9)]))
        self.assertEqual(bp.sge_choice(None, region), ("gene region", [("3", 1, 9)]))       # fixed now
        self.assertEqual(bp.sge_choice([], region), ("gene region", [("3", 1, 9)]))         # no file


class EndToEnd(unittest.TestCase):
    def test_build_writes_positions_evidence_and_sources(self):
        import pyarrow as pa
        import pyarrow.parquet as pq
        import mavedb_sources as ms
        import development_sources as s1
        with tempfile.TemporaryDirectory() as t:
            raw, exposure = Path(t) / "raw", Path(t) / "derived" / "exposure"
            raw.mkdir()
            write_gff(raw / "g.gff.gz")
            sentinel, expected_tg = 1000, {}
            for source, dirs in bp.TRAITGYM.items():
                for name in dirs:
                    sentinel += 1000                                  # a distinct position per file
                    (raw / "traitgym" / name).mkdir(parents=True)
                    pq.write_table(pa.table({"chrom": ["17"], "pos": [sentinel], "ref": ["A"], "alt": ["G"]}),
                                   raw / "traitgym" / name / "test.parquet")
                    expected_tg.setdefault(source, []).append(("17", sentinel, sentinel))
            (raw / "ncvardb").mkdir()
            expected_nc = []
            for i, name in enumerate(bp.NCVAR_TSV):
                pos = 600 + 10 * i
                (raw / "ncvardb" / name).write_text(f"ID\tgenome\tchr\tpos\tref\talt\nv{i}\thg38\t17\t{pos}\tA\tG\n"
                                                    f"w{i}\thg38\t17\t{pos + 100}\tC\tT\n", encoding="utf-8")
                expected_nc += [("17", pos, pos), ("17", pos + 100, pos + 100)]
            for i, name in enumerate(bp.NCVAR_VCF):
                (raw / "ncvardb" / name).write_text(f"#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                                                    f"17\t{700 + 10 * i}\t.\tC\tT\t50\tPASS\tNC=w{i}\n", encoding="utf-8")
            (raw / "organization_summary.txt").write_text("#organization\torganization ID\n", encoding="utf-8")
            (raw / "mavedb").mkdir()
            urn = "urn:mavedb:9-a-1"
            records = {("experiment", "urn:mavedb:9-a"): {"scoreSetUrns": [urn]},
                       ("score-set", urn): {"supersededScoreSet": None}}
            for (kind, u), body in records.items():
                (raw / ms.local_name(u, kind)).write_text(json.dumps(body), encoding="utf-8")
            (raw / ms.local_name(urn, "scores")).write_text(f"accession,score\n{urn}#1,1\n", encoding="utf-8")
            (raw / ms.local_name(urn, "mapped")).write_text(json.dumps([{"variantUrn": f"{urn}#1", "postMapped": {
                "type": "Allele", "location": {"start": 43050000, "end": 43050001,
                                               "sequenceReference": {"label": "NC_000017.11"}}}}]), encoding="utf-8")
            studies = ({"study": "S1", "gene": "BRCA1", "experiments": ("urn:mavedb:9-a",)},
                       {"study": "S2", "gene": "BRCA1", "experiments": None},
                       {"study": "S3", "gene": "PARG", "experiments": ("urn:mavedb:8-a",)})   # never obtained
            with mock.patch.multiple(bp, RAW=raw, EXPOSURE=exposure, GFF=raw / "g.gff.gz", ROOT=Path(t)), \
                    mock.patch.multiple(ms, RAW=raw, STUDIES=studies):
                result = bp.build()
            self.assertEqual([s["used"] for s in result["studies"]], ["positions", "gene region", "gene region"])
            self.assertEqual(result["studies"][2]["missing"], ["mavedb/urn_mavedb_8-a.experiment.json"])
            self.assertEqual((result["studies"][0]["score_sets"], result["studies"][0]["unpositioned"]), (1, 0))
            with contextlib.nullcontext(), \
                    contextlib.nullcontext():
                positions = idx.read_positions(exposure / "development_positions.tsv")
                sources = s1.read_sources(exposure / "development_sources.tsv")
            self.assertEqual(positions["saturation genome editing"],
                             [("17", 43043295, 43171327), ("X", 1, 3000), ("Y", 4000, 7000)])
            for source, intervals in expected_tg.items():
                self.assertEqual(positions[source], intervals)                     # every file contributes
            self.assertEqual(sorted(positions["ncVarDB"]), sorted(expected_nc))
            self.assertEqual(s1.coverage_problems(sources, exposure), [])
            evidence = json.loads((exposure / sources["saturation genome editing"][1]).read_text(encoding="utf-8"))
            self.assertEqual(evidence["studies"][0]["score_sets"], 1)
            self.assertEqual(evidence["studies"][0]["score_set_urns"], ["urn:mavedb:9-a-1"])
            self.assertEqual(evidence["studies"][0]["post_mapped_forms"], {"Allele": 1})
            # a second build from the same inputs is byte-identical.
            first = {f: (exposure / f).read_bytes() for f in ("development_positions.tsv", "development_sources.tsv")}
            again = Path(t) / "again" / "exposure"
            with mock.patch.multiple(bp, RAW=raw, EXPOSURE=again, GFF=raw / "g.gff.gz", ROOT=Path(t)), \
                    mock.patch.multiple(ms, RAW=raw, STUDIES=studies):
                bp.build()
            self.assertEqual({f: (again / f).read_bytes() for f in first}, first)
            for name in sorted(p.name for p in (exposure / "evidence").iterdir()):
                self.assertEqual((again / "evidence" / name).read_bytes(), (exposure / "evidence" / name).read_bytes())
            with self.assertRaises(FileExistsError):                                   # never overwritten
                with mock.patch.multiple(bp, RAW=raw, EXPOSURE=exposure, GFF=raw / "g.gff.gz", ROOT=Path(t)), \
                        mock.patch.multiple(ms, RAW=raw, STUDIES=studies):
                    bp.build()


if __name__ == "__main__":
    unittest.main()
