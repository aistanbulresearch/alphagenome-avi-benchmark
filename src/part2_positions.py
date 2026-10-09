"""Where the Part 2 variants lie relative to the exon their assay measured (a description of the data, not an
analysis): exonic, intronic 1-2 bases from that exon, 3-50, more than 50, flank (beyond the outer end of a
gene's first or last exon), other exon (MST1R exons 10 and 12), or unmapped. Distances come from each
dataset's own construct description (src/part2_exons.py). Writes results/part2/variant_positions.tsv.

Usage:
    python src/part2_positions.py [--gtf <gencode.v38.annotation.gtf.gz>]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

import part2_analysis as p2
import part2_exons as px
from paths import RAW, ROOT

CLASSES = ("exonic", "intronic 1-2", "intronic 3-50", "intronic >50", "flank", "other exon", "unmapped")


def table(d: pd.DataFrame, cls: pd.Series) -> pd.DataFrame:
    """Variants and splice-altering variants per dataset and class."""
    d = d.assign(**{"class": cls})
    return (d.groupby(["dataset", "class"]).agg(variants=("y", "size"), splice_altering=("y", "sum"))
            .reset_index())


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gtf", type=Path, default=RAW / "gencode.v38.annotation.gtf.gz")
    a = ap.parse_args(argv)
    d = p2.load()
    cls = px.classes(d, a.gtf)
    rows = []
    for name, mask in p2.analysis_sets(d).items():
        t = table(d[mask], cls[mask])
        t["dataset"] = name
        rows.append(t)
    out = pd.concat(rows, ignore_index=True)
    path = ROOT / "results" / "part2" / "variant_positions.tsv"
    out.to_csv(path, sep="\t", index=False)
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
