"""Additional Part 2 analysis, added after a code-consistency review of the manuscript (2026-10-11), defined
before it was computed: the analysed COMPASS set split by the authors' source subsets.

The planned COMPASS set (STUDY_DESIGN.md: random-SNV and ExAC subsets) holds the authors' random SNVs and their
ExAC SNVs, the latter chosen by presence in a population database. This analysis
- takes the subset labels as given in the harmonised assay table (`subset`: random_snv, exac);
- splits the COMPASS rows that enter the planned analysis (same exclusions) into the two subsets;
- computes, per subset, AUROC and AUPRC for AVI, SpliceAI raw and masked, CADD and phyloP 100-way, and the paired
  AUROC differences AVI minus each tool, with the gene bootstrap used for COMPASS (2,000 replicates, same seed);
- recomputes the random-effects pooled AVI minus SpliceAI raw and masked differences (DerSimonian-Laird and
  Hartung-Knapp) over the planned main datasets, with COMPASS replaced by its random-SNV subset;
- reports the ExAC subset per dataset only, descriptively; it never enters a pooled estimate.
The planned results (results/part2/*.tsv, pooled.json) are not changed.

Output: results/part2/compass_subsets_discrimination.tsv, compass_subsets_paired_differences.tsv,
compass_subsets_pooled.json.

Usage:
    python src/part2_compass_subsets.py [--replicates 2000]
"""
from __future__ import annotations

import argparse
import json
import sys

import pandas as pd

import metrics as m
import part2_analysis as p2

RANDOM = "compass, random SNVs"
EXAC = "compass, ExAC SNVs"
LABELS = {RANDOM: "random_snv", EXAC: "exac"}
TOOLS = ["AVI", "SpliceAI", "SpliceAI masked", "CADD", "phyloP 100-way"]
POOLED_PAIRS = ["AVI - SpliceAI", "AVI - SpliceAI masked"]


def subset_masks(d: pd.DataFrame, compass_mask: pd.Series) -> dict:
    """{name: mask} of the analysed COMPASS rows in each source subset."""
    return {name: compass_mask & (d["subset"] == label) for name, label in LABELS.items()}


def planned_rows(table: pd.DataFrame) -> list:
    """Planned per-dataset paired differences (results/part2/paired_differences.tsv) as numeric records."""
    return [{"dataset": r["dataset"], "pair": r["pair"], "difference": float(r["difference"]), "se": float(r["se"])}
            for r in table[["dataset", "pair", "difference", "se"]].to_dict("records")]


def pooled_replacing(planned: list, subsets: list, main_sets: list, pairs: list) -> dict:
    """Pooled AVI minus tool differences over the planned main datasets with COMPASS replaced by its random-SNV
    subset; the ExAC subset is not pooled."""
    out = {}
    for pair in pairs:
        rows = [r for r in planned if r["pair"] == pair and r["dataset"] in main_sets and r["dataset"] != "compass"
                and r["se"] > 0]
        rows += [r for r in subsets if r["pair"] == pair and r["dataset"] == RANDOM and r["se"] > 0]
        est, se = [r["difference"] for r in rows], [r["se"] for r in rows]
        out[pair] = {**p2.dersimonian_laird(est, se), "datasets": [r["dataset"] for r in rows]}
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--replicates", type=int, default=m.REPLICATES)
    a = ap.parse_args(argv)
    d = p2.load()
    sets = p2.analysis_sets(d)
    rows, diffs = [], []
    for name, mask in subset_masks(d, sets["compass"]).items():
        r, df = p2.dataset_tables(d, mask, name, TOOLS, a.replicates)
        rows += r
        diffs += df
    planned = planned_rows(pd.read_csv(p2.OUT / "paired_differences.tsv", sep="\t"))
    pooled_planned = json.loads((p2.OUT / "pooled.json").read_text(encoding="utf-8"))["pooled"]
    main_sets = pooled_planned["AVI - SpliceAI"]["datasets"]
    pooled = pooled_replacing(planned, diffs, main_sets, POOLED_PAIRS)
    counts = {n: {"variants": int(mk.sum()), "splice_altering": int(d.loc[mk, "y"].sum())}
              for n, mk in subset_masks(d, sets["compass"]).items()}
    pd.DataFrame(rows).to_csv(p2.OUT / "compass_subsets_discrimination.tsv", sep="\t", index=False)
    pd.DataFrame(diffs).to_csv(p2.OUT / "compass_subsets_paired_differences.tsv", sep="\t", index=False)
    (p2.OUT / "compass_subsets_pooled.json").write_text(
        json.dumps({"pooled_with_random_compass": pooled, "counts": counts}, indent=1) + "\n", encoding="utf-8")
    print("tables written to", p2.OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
