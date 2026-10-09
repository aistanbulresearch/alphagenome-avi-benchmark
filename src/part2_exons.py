"""Position of each Part 2 variant relative to the exon its assay measured, from each dataset's own description
of the construct, and the position class built on it. Replaces the nearest GENCODE exon, which can be another
transcript's exon or a neighbouring gene's.

    POU1F1, WT1    HGVS name in the source (c.N+k / c.N-k: k bases; no offset: exonic)
    CD19           exon intervals spanned by the variants the authors label exon 1-3; exon 2 assayed
    MST1R          exons 10-12 of the MANE Select transcript (the minigene holds them); exon 11 assayed
    COMPASS        the assayed exon's hg38 start and end, per sequence (Reference); kept only where both ends
                   are GENCODE v38 exon ends on that sequence
    OpenSplice     per exon_id, the substitutions the authors place in the exon region plus the exonic bases of
                   the splice-site regions (5'SS: last 3 exonic bases; 3'SS: first exonic base), kept only where
                   the span's length equals the table's exon_length
    MFASS          exon inside the 170-nt construct: left flank intron1_len, then exon_len (genomic orientation)

In the two multi-exon minigenes (CD19, MST1R) a variant inside an exon other than the assayed one is classed
"other exon", and an intronic variant's distance is to the nearest exon of the construct (as in Part 1).
Outside the outer end of a gene's first or last exon (MANE Select) a variant is classed "flank", not intronic.
"""
from __future__ import annotations

import csv
import gzip
import re
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

from paths import RAW

PART2 = RAW / "part2"
_HGVS = re.compile(r"c\.\*?-?\d+(?:([+-])(\d+))?")
SPLICE_SITE_EXONIC = {"5'SS": 3, "3'SS": 1}      # exonic bases inside each OpenSplice splice-site region


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


def nearest(pos: int, intervals: list) -> Optional[tuple]:
    if not intervals:
        return None
    return min(intervals, key=lambda iv: interval_distance(pos, [iv]))


def mfass_exon(row) -> tuple:
    start = int(row["start_hg38_0based"]) + int(row["intron1_len"]) + 1
    return start, start + int(row["exon_len"]) - 1


def intervals_from_labels(d: pd.DataFrame, key: str, pos: str, is_exon: str) -> dict:
    """{key: [(first, last)]}: the span of the positions labelled exonic, per construct."""
    ex = d[d[is_exon]]
    return {k: [(int(g[pos].min()), int(g[pos].max()))] for k, g in ex.groupby(key)}


def distance_class(dist: int | None) -> str:
    if dist is None or pd.isna(dist):
        return "unmapped"
    if dist == 0:
        return "exonic"
    return "intronic 1-2" if dist <= 2 else "intronic 3-50" if dist <= 50 else "intronic >50"


def position_class(dist: int | None, override: str | None) -> str:
    """The class: an override ("flank", "other exon") where set, else the class of the distance."""
    return override if isinstance(override, str) else distance_class(dist)


def opensplice_spans(t: pd.DataFrame) -> dict:
    """{exon_id: [(start, end)]} from OpenSplice substitutions (columns exon_id, region, POS, exon_length):
    the exon-region span extended by the exonic bases of each splice-site region on the side where that region
    lies (so on either strand). Exons without exon-region rows or whose span length differs from exon_length
    are left out."""
    out = {}
    for e, g in t.groupby("exon_id"):
        ex = g.loc[g["region"] == "Exon", "POS"]
        length = g["exon_length"].dropna()
        if ex.empty or length.empty:
            continue
        s, t_end = int(ex.min()), int(ex.max())
        for region, n in SPLICE_SITE_EXONIC.items():
            ss = g.loc[g["region"] == region, "POS"]
            if ss.empty:
                continue
            if int(ss.max()) < s:
                s -= n
            elif int(ss.min()) > t_end:
                t_end += n
        if t_end - s + 1 == int(length.iloc[0]):
            out[e] = [(s, t_end)]
    return out


def construct_position(pos: int, exons: list, assayed: int) -> tuple:
    """(distance, override) in a multi-exon construct; `exons` in genomic order, `assayed` the index of the
    measured exon. Inside another exon: "other exon"; otherwise the distance to the nearest exon."""
    if any(s <= pos <= e for i, (s, e) in enumerate(exons) if i != assayed):
        return None, "other exon"
    return interval_distance(pos, exons), None


def outer_sides(lines: Iterable[str]) -> tuple[dict, dict]:
    """From GTF exon lines of MANE Select transcripts: per sequence, the ends of exons whose low side lies outside
    the gene (first exon on +, last on -) and the starts of exons whose high side does (last on +, first on -)."""
    tx = defaultdict(list)
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon" or 'tag "MANE_Select"' not in f[8]:
            continue
        tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        tx[tid].append((f[0], int(f[3]), int(f[4])))
    low, high = defaultdict(set), defaultdict(set)
    for exons in tx.values():
        exons.sort(key=lambda x: x[1])
        low[exons[0][0]].add(exons[0][2])
        high[exons[-1][0]].add(exons[-1][1])
    return low, high


def is_flank(pos: int, span: tuple, low_ends: set, high_starts: set) -> bool:
    """True if `pos` lies beyond the outer end of an exon that is a gene's first or last exon."""
    s, e = span
    return (pos < s and e in low_ends) or (pos > e and s in high_starts)


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
        d = pd.read_csv(h, sep="\t", usecols=["POS", "exon_id", "region", "mut_type", "exon_length"])
    return opensplice_spans(d[d["mut_type"] == "sub"])


def mfass_exons() -> dict:
    files = list((PART2 / "mfass").glob("github_*/snv_data_clean.txt"))
    d = pd.read_csv(files[0], sep="\t", usecols=["ensembl_id", "start_hg38_0based", "intron1_len", "exon_len"])
    return {k: mfass_exon(g.iloc[0]) for k, g in d.groupby("ensembl_id")}


def gtf_tables(gtf: Path) -> tuple:
    """(MST1R exons 10-12, outer sides (low, high), exon starts and ends per sequence) from one GTF pass."""
    mst1r, mane, starts, ends = [], [], defaultdict(set), defaultdict(set)
    with gzip.open(gtf, "rt", encoding="utf-8") as h:
        for line in h:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "exon":
                continue
            starts[f[0]].add(int(f[3]))
            ends[f[0]].add(int(f[4]))
            if 'tag "MANE_Select"' in f[8]:
                mane.append(line)
                if 'gene_name "MST1R"' in f[8] and re.search(r'exon_number (\d+);', f[8]).group(1) in ("10", "11", "12"):
                    mst1r.append((int(f[3]), int(f[4])))
    return sorted(mst1r), outer_sides(mane), starts, ends


def positions(d: pd.DataFrame, gtf: Path) -> pd.DataFrame:
    """Per row of the Part 2 table: distance to the assayed exon ('dist', None where it cannot be placed) and an
    override class ('override': "flank", "other exon" or None)."""
    n = len(d)
    note, pos = d["region_note"].fillna("").tolist(), d["pos"].astype(int).tolist()
    ds, src = d["dataset"].tolist(), d["variant_id_in_source"].fillna("").astype(str).tolist()
    chrom = [c if c.startswith("chr") else "chr" + c for c in d["chrom"].astype(str)]
    mexons, (low, high), starts, ends = gtf_tables(gtf)
    dist, override, spans = [None] * n, [None] * n, [[] for _ in range(n)]   # spans: assayed exon interval(s)
    cd_rows = [i for i in range(n) if ds[i] == "cortes_lopez2022_cd19"]
    cd_spans = defaultdict(list)
    for i in cd_rows:
        label = note[i].split(";")[0].strip()
        if label.startswith("exon"):
            cd_spans[label].append(pos[i])
    cd19 = sorted((min(v), max(v)) for v in cd_spans.values())
    cex, oex, mex = compass_exons(), opensplice_exons(), mfass_exons()
    for i in range(n):
        if ds[i] in ("kitzman_pou1f1_ex2", "kitzman_wt1_ex9"):
            dist[i] = hgvs_distance(note[i])
            continue
        if ds[i] == "braun2018_mst1r_ex11":
            dist[i], override[i] = construct_position(pos[i], mexons, 1)
            spans[i] = mexons
        elif ds[i] == "cortes_lopez2022_cd19":
            dist[i], override[i] = construct_position(pos[i], cd19, 1)
            spans[i] = cd19
        elif ds[i] == "compass":
            sp = cex.get(src[i].split(":")[0])
            if sp and sp[0] in starts[chrom[i]] and sp[1] in ends[chrom[i]]:
                spans[i] = [sp]
        elif ds[i] == "opensplice":
            parts = note[i].split(";")
            spans[i] = oex.get(parts[1].strip(), []) if len(parts) > 1 else []
        elif ds[i] == "mfass":
            sp = mex.get(src[i].rsplit("_", 1)[0])
            spans[i] = [sp] if sp else []
        else:
            continue
        if override[i] is None:
            dist[i] = interval_distance(pos[i], spans[i])
            iv = nearest(pos[i], spans[i])
            if iv is not None and is_flank(pos[i], iv, low[chrom[i]], high[chrom[i]]):
                override[i] = "flank"
    return pd.DataFrame({"dist": pd.Series(dist, index=d.index, dtype=object),
                         "override": pd.Series(override, index=d.index, dtype=object)})


def classes(d: pd.DataFrame, gtf: Path) -> pd.Series:
    p = positions(d, gtf)
    return pd.Series([position_class(x, o) for x, o in zip(p["dist"], p["override"])], index=d.index)


__all__ = ["hgvs_distance", "interval_distance", "nearest", "mfass_exon", "intervals_from_labels", "distance_class",
           "position_class", "opensplice_spans", "construct_position", "outer_sides", "is_flank", "gtf_tables",
           "positions", "classes"]
