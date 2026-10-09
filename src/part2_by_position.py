"""Part 2 discrimination by position relative to the exon each assay measured (exonic, intronic 3-50, intronic
>50; the same 50-base cut as Part 1's deep intronic region), per dataset and pooled over the main datasets. A dataset's
class is analysed when it holds at least 10 splice-altering and 10 other variants. Added after the main
Part 2 results. Writes results/part2/by_position_*.

Usage:
    python src/part2_by_position.py [--replicates 2000] [--gtf <gencode.v38.annotation.gtf.gz>]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

import metrics as m
import part2_analysis as p2
import part2_exons as px
from paths import RAW

CLASSES = ("exonic", "intronic 3-50", "intronic >50")
MIN_EACH = 10


def by_position(d: pd.DataFrame, sets: dict, tools: list, reps: int, main_sets: list) -> tuple[list, list, dict]:
    """Rows and AVI − tool differences per (dataset, class), and pooled differences per class."""
    rows, diffs, pooled = [], [], {}
    for cls in CLASSES:
        class_diffs = []
        for ds, mask in sets.items():
            mk = mask & (d["class"] == cls)
            y = d.loc[mk, "y"]
            if (y == 1).sum() < MIN_EACH or (y == 0).sum() < MIN_EACH:
                continue
            r, df = p2.dataset_tables(d, mk, ds, tools, reps)
            rows += [{**x, "class": cls} for x in r]
            class_diffs += [{**x, "class": cls} for x in df]
        diffs += class_diffs
        pooled[cls] = p2.pooled_differences(class_diffs, main_sets, tools)
    return rows, diffs, pooled


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--replicates", type=int, default=m.REPLICATES)
    ap.add_argument("--gtf", type=Path, default=RAW / "gencode.v38.annotation.gtf.gz")
    a = ap.parse_args(argv)
    d = p2.load()
    d["class"] = px.classes(d, a.gtf)                 # distance to the exon each assay measured
    tools = [t for t in p2.TOOLS if t in d.columns]
    sets = p2.analysis_sets(d)
    main_sets = [n for n in sets if n != "mfass" and not n.startswith("opensplice, ")]
    rows, diffs, pooled = by_position(d, sets, tools, a.replicates, main_sets)
    p2.OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(p2.OUT / "by_position_discrimination.tsv", sep="\t", index=False)
    pd.DataFrame(diffs).to_csv(p2.OUT / "by_position_differences.tsv", sep="\t", index=False)
    (p2.OUT / "by_position_pooled.json").write_text(json.dumps(pooled, indent=1) + "\n", encoding="utf-8")
    print("tables written to", p2.OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
