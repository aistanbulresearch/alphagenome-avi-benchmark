"""Part 1 supplementary discrimination: paired AUPRC differences and sensitivity at fixed specificity
(95% and 99%), with gene-bootstrap intervals, AVI against each comparator. Writes
results/part1/supplementary_discrimination.tsv.

Usage:
    python src/part1_supplement.py [--replicates 2000]
"""
from __future__ import annotations

import argparse
import sys

import pandas as pd

import metrics as m
import part1_data as pdata
from part1_analysis import OUT, base_weights

COMPARATORS = ("SpliceAI", "SpliceAI 4999", "Pangolin", "CADD")
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
    return 0


if __name__ == "__main__":
    sys.exit(main())
