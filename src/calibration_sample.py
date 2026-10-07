"""The threshold-derivation set for the calibration: rows of the all-labels set outside the new-label set,
non-canonical intronic and UTR regions (canonical +/-1, +/-2 sites are evidence for a loss-of-function
criterion, not for a predictor criterion), at least one star, chromosomes 1-22 and X.

Every P/LP row is kept. B/LB rows are sampled at random, stratified by region and chromosome group with
proportional allocation, to a fixed total; the calibration weights the two classes to the prior
probability, so the B/LB sampling fraction does not change the posterior. The sample is drawn once,
with a fixed seed, before any score is retrieved, and its list is committed.

Output: data/derived/variants/calibration_sample.tsv (variation_id, label, region, chromosome_group).

Usage:
    python src/calibration_sample.py [--benign 40000] [--seed 20261007]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import random
import sys
from collections import defaultdict
from pathlib import Path

import gnomad_frequencies as gf
from paths import DERIVED

OUT = DERIVED / "variants" / "calibration_sample.tsv"


def eligible(row: dict) -> bool:
    return (row["chrom"] in gf.CHROMS and row["stars"].isdigit() and int(row["stars"]) >= 1
            and row["family"] in ("intronic", "utr") and row["region"] != "canonical"
            and row["new_label"] == "no")


def draw(rows: list, benign: int, seed: int) -> list:
    """All P/LP rows and a proportional stratified random sample of `benign` B/LB rows."""
    pathogenic = [r for r in rows if r["label"] == "P/LP"]
    strata = defaultdict(list)
    for r in rows:
        if r["label"] == "B/LB":
            strata[(r["region"], r["chromosome_group"])].append(r)
    total = sum(len(v) for v in strata.values())
    rng = random.Random(seed)
    chosen = []
    keys = sorted(strata)
    quotas = {k: len(strata[k]) * benign / total for k in keys}
    counts = {k: int(quotas[k]) for k in keys}
    # largest remainders get the rest, ties by key order
    for k in sorted(keys, key=lambda k: (-(quotas[k] - counts[k]), k))[: benign - sum(counts.values())]:
        counts[k] += 1
    for k in keys:
        pool = sorted(strata[k], key=lambda r: int(r["variation_id"]))
        chosen += rng.sample(pool, min(counts[k], len(pool)))
    return sorted(pathogenic + chosen, key=lambda r: int(r["variation_id"]))


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--benign", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=20261007)
    a = ap.parse_args(argv)
    with gzip.open(gf.VARIANTS, "rt", encoding="utf-8", newline="") as h:
        rows = [r for r in csv.DictReader(h, delimiter="\t") if eligible(r)]
    sample = draw(rows, a.benign, a.seed)
    with open(OUT, "w", encoding="utf-8", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(("variation_id", "label", "region", "chromosome_group"))
        for r in sample:
            w.writerow((r["variation_id"], r["label"], r["region"], r["chromosome_group"]))
    n_p = sum(r["label"] == "P/LP" for r in sample)
    print(f"eligible {len(rows)}; sample {len(sample)} = {n_p} P/LP + {len(sample) - n_p} B/LB; seed {a.seed}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
