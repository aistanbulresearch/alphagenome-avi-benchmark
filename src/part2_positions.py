"""Where the Part 2 variants lie relative to exons (a description of the data, not an analysis): exonic,
intronic 1-2 bases from an exon (canonical splice site), 3-50, or more than 50 bases. Exons are those of the
protein-coding transcripts of the variant's gene in GENCODE v38 (any protein-coding gene if the gene name is
not found). Writes results/part2/variant_positions.tsv (counts per dataset and class).

Usage:
    python src/part2_positions.py [--gtf <gencode.v38.annotation.gtf.gz>]
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

import part2_analysis as p2
from paths import RAW, ROOT

CLASSES = ("exonic", "intronic 1-2", "intronic 3-50", "intronic >50")


def read_exons(gtf: Path) -> dict:
    """{"by_gene": {(chrom, gene_name): [(start, end)]}, "all": {chrom: sorted [(start, end)]}} for exons of
    protein-coding transcripts; chromosomes without the 'chr' prefix."""
    by_gene, all_ = defaultdict(set), defaultdict(set)
    with gzip.open(gtf, "rt", encoding="utf-8") as h:
        for line in h:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] != "exon" or 'transcript_type "protein_coding"' not in f[8]:
                continue
            chrom = f[0][3:] if f[0].startswith("chr") else f[0]
            iv = (int(f[3]), int(f[4]))
            name = re.search(r'gene_name "([^"]+)"', f[8])
            if name:
                by_gene[(chrom, name.group(1))].add(iv)
            all_[chrom].add(iv)
    return {"by_gene": {k: sorted(v) for k, v in by_gene.items()}, "all": {k: sorted(v) for k, v in all_.items()}}


def _distance(intervals: list, pos: int) -> int | None:
    """0 inside an exon, else bases to the nearest exon end; None without exons."""
    if not intervals:
        return None
    best = None
    i = bisect.bisect_right(intervals, (pos, float("inf")))
    near = intervals if len(intervals) <= 5000 else intervals[max(0, i - 200):i + 200]   # one gene: all exons
    for s, e in near:
        dist = 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e))
        best = dist if best is None else min(best, dist)
    return best


def classify(exons: dict, chrom: str, pos: int, gene: str) -> str:
    intervals = exons["by_gene"].get((chrom, gene)) or exons["all"].get(chrom, [])
    dist = _distance(intervals, pos)
    if dist is None:
        return "no exon"
    if dist == 0:
        return "exonic"
    return "intronic 1-2" if dist <= 2 else "intronic 3-50" if dist <= 50 else "intronic >50"


def table(d: pd.DataFrame, exons: dict) -> pd.DataFrame:
    d = d.assign(**{"class": [classify(exons, str(c), int(p), str(g)) for c, p, g in
                              zip(d["chrom"], d["pos"], d["gene"])]})
    return (d.groupby(["dataset", "class"]).agg(variants=("y", "size"), splice_altering=("y", "sum"))
            .reset_index())


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gtf", type=Path, default=RAW / "gencode.v38.annotation.gtf.gz")
    a = ap.parse_args(argv)
    d = p2.load()
    rows = []
    exons = read_exons(a.gtf)
    for name, mask in p2.analysis_sets(d).items():
        t = table(d[mask], exons)
        t["dataset"] = name
        rows.append(t)
    out = pd.concat(rows, ignore_index=True)
    path = ROOT / "results" / "part2" / "variant_positions.tsv"
    out.to_csv(path, sep="\t", index=False)
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
