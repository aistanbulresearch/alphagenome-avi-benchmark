"""Distance of each Part 2 variant to the exon its assay measured (0 inside that exon), from each dataset's own
description of the construct, and the position class built on it. Replaces the nearest GENCODE exon, which can
be another transcript's exon or a neighbouring gene's.

    POU1F1, WT1    HGVS name in the source (c.N+k / c.N-k: k bases; no offset: exonic)
    CD19           exon intervals spanned by the variants the authors label exon 1-3
    MST1R          exons 10-12 of the MANE Select transcript (the minigene holds exons 10-12)
    COMPASS        the assayed exon's hg38 start and end, per sequence (Reference)
    OpenSplice     exon interval spanned by the variants the authors place in the exon region, per exon_id
    MFASS          exon inside the 170-nt construct: left flank intron1_len, then exon_len (genomic orientation)
"""
from __future__ import annotations

import csv
import gzip
import re
import zipfile
from collections import defaultdict
from pathlib import Path

import pandas as pd

from paths import RAW

PART2 = RAW / "part2"
_HGVS = re.compile(r"c\.\*?-?\d+(?:([+-])(\d+))?")


def hgvs_distance(note: str) -> int | None:
    """Bases into the intron from an HGVS c. name in a region note; 0 for an exonic name; None if absent."""
    m = _HGVS.search(note or "")
    if not m:
        return None
    return int(m.group(2)) if m.group(1) else 0


def interval_distance(pos: int, intervals: list) -> int | None:
    if not intervals:
        return None
    return min(0 if s <= pos <= e else min(abs(pos - s), abs(pos - e)) for s, e in intervals)


def mfass_exon(row) -> tuple:
    start = int(row["start_hg38_0based"]) + int(row["intron1_len"]) + 1
    return start, start + int(row["exon_len"]) - 1


def intervals_from_labels(d: pd.DataFrame, key: str, pos: str, is_exon: str) -> dict:
    """{key: [(first, last)]}: the span of the positions labelled exonic, per construct."""
    ex = d[d[is_exon]]
    return {k: [(int(g[pos].min()), int(g[pos].max()))] for k, g in ex.groupby(key)}


def distance_class(dist: int | None) -> str:
    if dist is None:
        return "unmapped"
    if dist == 0:
        return "exonic"
    return "intronic 1-2" if dist <= 2 else "intronic 3-50" if dist <= 50 else "intronic >50"


def compass_exons() -> dict:
    out = {}
    with gzip.open(PART2 / "compass" / "GSE307247_Processed_PSIs_All_Cells.csv.gz", "rt", encoding="utf-8") as h:
        for r in csv.DictReader(h):
            if r["exon_start_hg38"] and r["exon_end_hg38"]:
                out[r["Reference"]] = (int(r["exon_start_hg38"]), int(r["exon_end_hg38"]))
    return out


def opensplice_exons() -> dict:
    z = zipfile.ZipFile(PART2 / "opensplice" / "Supplementary_Table4_media-6.zip")
    with z.open("Supplementary_Table4.tsv") as h:
        d = pd.read_csv(h, sep="\t", usecols=["POS", "exon_id", "region", "mut_type"])
    d = d[d["mut_type"] == "sub"]
    return intervals_from_labels(d.assign(is_exon=d["region"] == "Exon"), "exon_id", "POS", "is_exon")


def mfass_exons() -> dict:
    files = list((PART2 / "mfass").glob("github_*/snv_data_clean.txt"))
    d = pd.read_csv(files[0], sep="\t", usecols=["ensembl_id", "start_hg38_0based", "intron1_len", "exon_len"])
    return {k: mfass_exon(g.iloc[0]) for k, g in d.groupby("ensembl_id")}


def mst1r_exons(gtf: Path) -> list:
    exons = []
    with gzip.open(gtf, "rt", encoding="utf-8") as h:
        for line in h:
            if 'gene_name "MST1R"' in line and 'tag "MANE_Select"' in line:
                f = line.split("\t")
                if f[2] == "exon" and re.search(r'exon_number (\d+);', f[8]).group(1) in ("10", "11", "12"):
                    exons.append((int(f[3]), int(f[4])))
    return sorted(exons)


def distances(d: pd.DataFrame, gtf: Path) -> pd.Series:
    """Distance to the assayed exon for every row of the Part 2 table (None where it cannot be placed)."""
    note, pos = d["region_note"].fillna(""), d["pos"].astype(int)
    out = pd.Series([None] * len(d), index=d.index, dtype=object)
    kz = d["dataset"].isin(["kitzman_pou1f1_ex2", "kitzman_wt1_ex9"])
    out[kz] = [hgvs_distance(n) for n in note[kz]]
    cd = d["dataset"] == "cortes_lopez2022_cd19"
    lab = note[cd].str.split(";").str[0].str.strip()
    spans = defaultdict(list)
    for (_, g) in pd.DataFrame({"lab": lab, "pos": pos[cd]}).groupby("lab"):
        if g["lab"].iloc[0].startswith("exon"):
            spans["cd19"].append((int(g["pos"].min()), int(g["pos"].max())))
    out[cd] = [interval_distance(p, spans["cd19"]) for p in pos[cd]]
    ms = d["dataset"] == "braun2018_mst1r_ex11"
    mexons = mst1r_exons(gtf)
    out[ms] = [interval_distance(p, mexons) for p in pos[ms]]
    co = d["dataset"] == "compass"
    cex = compass_exons()
    out[co] = [interval_distance(p, [cex[r]] if r in cex else [])
               for p, r in zip(pos[co], d.loc[co, "variant_id_in_source"].str.split(":").str[0])]
    op = d["dataset"] == "opensplice"
    oex = opensplice_exons()
    out[op] = [interval_distance(p, oex.get(e, []))
               for p, e in zip(pos[op], note[op].str.split(";").str[1].str.strip())]
    mf = d["dataset"] == "mfass"
    mex = mfass_exons()
    out[mf] = [interval_distance(p, [mex[e]] if e in mex else [])
               for p, e in zip(pos[mf], d.loc[mf, "variant_id_in_source"].str.rsplit("_", n=1).str[0])]
    return out


def classes(d: pd.DataFrame, gtf: Path) -> pd.Series:
    return distances(d, gtf).map(distance_class)


__all__ = ["hgvs_distance", "interval_distance", "mfass_exon", "intervals_from_labels", "distance_class",
           "distances", "classes"]
