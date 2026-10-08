"""Part 1 supplementary tables: paired AUPRC differences and sensitivity at fixed specificity (95% and
99%), with gene-bootstrap intervals, AVI against each comparator; and published thresholds (CADD for
missense, SpliceAI) evaluated as score intervals on the main set. Writes
results/part1/supplementary_discrimination.tsv and results/part1/published_thresholds.tsv.

Usage:
    python src/part1_supplement.py [--replicates 2000]
"""
from __future__ import annotations

import argparse
import sys

import pandas as pd

import calibration as cal
import metrics as m
import part1_data as pdata
from part1_analysis import OUT, base_weights

COMPARATORS = ("SpliceAI", "SpliceAI 4999", "Pangolin", "SpliceAI2", "CADD")
SETS = ("main", "predictor not mentioned")
REGIONS = ("intronic >=3 bp (main)", "donor region", "acceptor region", "near intronic", "deep intronic", "UTR")
SPECIFICITIES = (0.95, 0.99)


def _metrics(s, y, w, prep) -> dict:
    out = {"AUPRC": m.auprc(s, y, w, prep)}
    for spec in SPECIFICITIES:
        out[f"sensitivity at {spec:.0%} specificity"] = m.sensitivity_at_specificity(s, y, spec, w, prep)
    return out


def paired(v: pd.DataFrame, mask: pd.Series, a: str, b: str, reps: int) -> list:
    """Rows (metric, value_a, value_b, difference, ci) on the variants both tools score. Design weights
    of sampled tools (w_<tool>) enter every metric, so the values refer to the full set."""
    d = v[mask & v[a].notna() & v[b].notna()]
    y = d["y"].to_numpy()
    if y.sum() == 0 or (1 - y).sum() == 0:
        return []
    sa, sb = d[a].to_numpy(), d[b].to_numpy()
    pa, pb, base = m._order(sa), m._order(sb), base_weights(d, (a, b))
    va, vb = _metrics(sa, y, base, pa), _metrics(sb, y, base, pb)
    boots = {k: [] for k in va}
    for w in m.cluster_weights(m.clusters_of(d["gene_id"].tolist()), reps):
        ba, bb = _metrics(sa, y, w * base, pa), _metrics(sb, y, w * base, pb)
        for k in boots:
            boots[k].append(ba[k] - bb[k])
    return [{"metric": k, "n_pos": int(y.sum()), "n_neg": int((1 - y).sum()), "value_a": va[k], "value_b": vb[k],
             "difference": va[k] - vb[k], "ci": m.interval(boots[k])} for k in va]


INF = float("inf")
# Published single-tool thresholds, as score intervals: (level, side, lo, hi); "ge" is lo <= s < hi and
# "le" is lo < s <= hi. CADD: Pejaver et al. 2022 Table 2 (derived for missense variants). SpliceAI:
# Walker et al. 2023 (derived on raw scores with a 4,999-base window).
PUBLISHED = {
    "CADD (Pejaver 2022, missense)": [("PP3 supporting", "ge", 25.3, 28.1), ("PP3 moderate", "ge", 28.1, INF),
                                      ("BP4 supporting", "le", 17.3, 22.7), ("BP4 moderate", "le", 0.15, 17.3),
                                      ("BP4 strong", "le", -INF, 0.15)],
    "SpliceAI (Walker 2023)": [("PP3", "ge", 0.2, INF), ("BP4", "le", -INF, 0.1)],
}


def published_intervals(v: pd.DataFrame, reps: int) -> list:
    """Interval likelihood ratios of published thresholds on the calibration test set of the main set
    (non-canonical intronic and UTR), with gene-bootstrap intervals."""
    test = pdata.sets(v)["main"] & v["region"].isin(list(pdata.MAIN_REGION) + ["utr5", "utr3"])
    rows = []
    for label, intervals in PUBLISHED.items():
        tool = label.split(" (")[0]
        if tool not in v.columns:
            continue
        d = v[test & v[tool].notna()]
        s, y = d[tool].to_numpy(), d["y"].to_numpy()
        clusters = m.clusters_of(d["gene_id"].tolist())
        for level, side, lo, hi in intervals:
            hit = (s >= lo) & (s < hi) if side == "ge" else (s > lo) & (s <= hi)
            boots = [cal.interval_lr(s, y, lo, hi, side, w) for w in m.cluster_weights(clusters, reps)]
            rows.append({"tool": label, "level": level, "side": side, "lo": lo, "hi": hi,
                         "lr": cal.interval_lr(s, y, lo, hi, side), "ci": m.interval(boots),
                         "n_pos": int(hit[y == 1].sum()), "n_neg": int(hit[y == 0].sum()),
                         "fraction_pos": float(hit[y == 1].mean()), "fraction_neg": float(hit[y == 0].mean())})
    return rows


def tables(v: pd.DataFrame, tools: list, reps: int) -> pd.DataFrame:
    set_masks, region_masks = pdata.sets(v), pdata.regions(v)
    rows = []
    for sname in SETS:
        for rname in REGIONS:
            for b in [t for t in COMPARATORS if t in tools]:
                for r in paired(v, set_masks[sname] & region_masks[rname], "AVI", b, reps):
                    rows.append({"set": sname, "region": rname, "pair": f"AVI - {b}", **r})
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--replicates", type=int, default=m.REPLICATES)
    a = ap.parse_args(argv)
    v = pdata.load()
    out = tables(v, pdata.tools_present(v), a.replicates)
    OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT / "supplementary_discrimination.tsv", sep="\t", index=False)
    print("written", OUT / "supplementary_discrimination.tsv", len(out), "rows")
    pd.DataFrame(published_intervals(v, a.replicates)).to_csv(OUT / "published_thresholds.tsv", sep="\t", index=False)
    print("written", OUT / "published_thresholds.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
