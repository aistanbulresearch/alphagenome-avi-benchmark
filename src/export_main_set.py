"""The Part 1 main analysis set as a list of identities and labels (no tool scores), so that readers can use the
same variants without rebuilding the ClinVar variant list. The set is the 'main' mask of part1_data.sets.

Output: data/derived/variants/main_set.tsv (variation_id, chrom, pos, ref, alt, label, region, chromosome_group).

Usage:
    python src/export_main_set.py
"""
from __future__ import annotations

import sys

import pandas as pd

import part1_data as p1

OUT = p1.VARIANTS / "main_set.tsv"
COLUMNS = ("variation_id", "chrom", "pos", "ref", "alt", "label", "region", "chromosome_group")


def main_rows(v: pd.DataFrame, mask: pd.Series) -> pd.DataFrame:
    """Rows of the main set, identity and label columns only."""
    return v.loc[mask.values, list(COLUMNS)].reset_index(drop=True)


def main(argv: list[str] | None = None) -> int:
    v = p1.load()
    out = main_rows(v, p1.sets(v)["main"])
    out = out.sort_values("variation_id", key=lambda s: s.astype(int))
    out.to_csv(OUT, sep="\t", index=False, lineterminator="\n")
    counts = out["label"].value_counts().to_dict()
    print(f"{len(out)} variants written to {OUT}: {counts}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
