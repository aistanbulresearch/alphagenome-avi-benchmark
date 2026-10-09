"""Part 1 analysis (STUDY_DESIGN.md section 4.4): coverage, discrimination with gene-bootstrap intervals,
paired differences, and the ClinGen calibration. Writes tables to results/part1/.

Usage:
    python src/part1_analysis.py [--replicates 2000] [--calibration-replicates 10000]
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd

import calibration as cal
import metrics as m
import part1_data as pdata
from paths import ROOT

OUT = ROOT / "results" / "part1"
PAIRS = [("AVI", "SpliceAI"), ("AVI", "CADD"), ("AVI", "SpliceAI 4999"), ("AVI", "Pangolin"), ("AVI", "SpliceAI2"),
         ("AVI", "phyloP 100-way"), ("AVI", "phastCons 470-way"), ("AVI", "phyloP 241-way")]


def coverage(v: pd.DataFrame, set_masks: dict, tools: list) -> pd.DataFrame:
    rows = []
    for name, mask in set_masks.items():
        sub = v[mask]
        for t in tools:
            rows.append({"set": name, "tool": t, "variants": len(sub), "scored": int(sub[t].notna().sum())})
    return pd.DataFrame(rows)


def base_weights(d: pd.DataFrame, tools: tuple) -> np.ndarray:
    """Design weights of a variant: the inverse B/LB sampling fraction for tools scored on a sample
    (column w_<tool>), else 1."""
    w = np.ones(len(d))
    for t in tools:
        if f"w_{t}" in d.columns:
            w = w * d[f"w_{t}"].fillna(1.0).to_numpy()
    return w


def discrimination(v: pd.DataFrame, mask: pd.Series, tools: list, reps: int) -> tuple[list, list]:
    """Per tool: AUROC, AUPRC with intervals; per pair: AUROC difference with interval (common variants).
    Tools scored on a B/LB sample carry inverse-sampling weights, so AUPRC refers to the full set."""
    sub = v[mask]
    rows, diffs = [], []
    for t in tools:
        d = sub[sub[t].notna()]
        y, s = d["y"].to_numpy(), d[t].to_numpy()
        if y.sum() == 0 or (1 - y).sum() == 0:
            rows.append({"tool": t, "n_pos": int(y.sum()), "n_neg": int((1 - y).sum())})
            continue
        prep, base = m._order(s), base_weights(d, (t,))
        boots_roc, boots_pr = [], []
        for w in m.cluster_weights(m.clusters_of(d["gene_id"].tolist()), reps):
            boots_roc.append(m.auroc(s, y, w * base, prep))
            boots_pr.append(m.auprc(s, y, w * base, prep))
        rows.append({"tool": t, "n_pos": int(y.sum()), "n_neg": int((1 - y).sum()),
                     "genes_pos": int(d.loc[d["y"] == 1, "gene_id"].nunique()),
                     "pos_fraction": float((base * y).sum() / base.sum()),
                     "auroc": m.auroc(s, y, base, prep), "auroc_ci": m.interval(boots_roc),
                     "auprc": m.auprc(s, y, base, prep), "auprc_ci": m.interval(boots_pr)})
    for a, b in PAIRS:
        if a not in tools or b not in tools:
            continue
        d = sub[sub[a].notna() & sub[b].notna()]
        y = d["y"].to_numpy()
        if y.sum() == 0 or (1 - y).sum() == 0:
            continue
        sa, sb = d[a].to_numpy(), d[b].to_numpy()
        pa, pb, base = m._order(sa), m._order(sb), base_weights(d, (a, b))
        boots = [m.auroc(sa, y, w * base, pa) - m.auroc(sb, y, w * base, pb)
                 for w in m.cluster_weights(m.clusters_of(d["gene_id"].tolist()), reps)]
        diffs.append({"pair": f"{a} - {b}", "n_pos": int(y.sum()), "n_neg": int((1 - y).sum()),
                      "difference": m.auroc(sa, y, base, pa) - m.auroc(sb, y, base, pb),
                      "ci": m.interval(boots)})
    return rows, diffs


def lr_interval(boots: list, hit: np.ndarray, y: np.ndarray) -> tuple[list, str]:
    """95% interval of an interval or threshold likelihood ratio: the gene bootstrap, or exact Clopper-Pearson
    bounds when the interval holds no P/LP or no B/LB (the bootstrap cannot vary an empty class)."""
    k_pos, k_neg = int(hit[y == 1].sum()), int(hit[y == 0].sum())
    if k_pos == 0 or k_neg == 0:
        return list(cal.exact_lr_interval(k_pos, int((y == 1).sum()), k_neg, int((y == 0).sum()))),             "exact (a zero count in the interval)"
    return list(m.interval(boots)), "gene bootstrap"


def calibration_tables(v: pd.DataFrame, set_masks: dict, tools: list, reps: int, boot_reps: int) -> dict:
    """Thresholds derived on the calibration sample (non-canonical intronic and UTR, rare) and their
    likelihood ratios on the main set (same regions), plus SpliceAI's published thresholds."""
    calib_regions = v["region"].isin(list(pdata.MAIN_REGION) + ["utr5", "utr3"])
    derive = set_masks["calibration sample"] & calib_regions
    test = set_masks["main"] & calib_regions
    out = {}
    for t in [x for x in tools if x not in ("SpliceAI 4999", "Distance")]:
        d = v[derive & v[t].notna()]
        if d.empty:                     # scored on the main set only (Pangolin): nothing to derive from
            continue
        res = cal.calibrate(d[t].to_numpy(), d["y"].to_numpy(), boot_reps)
        tv = v[test & v[t].notna()]
        s, y = tv[t].to_numpy(), tv["y"].to_numpy()
        clusters = m.clusters_of(tv["gene_id"].tolist())
        evaluated = {"interval": {}, "cumulative": {}}
        # main test (ClinGen logic): each level within its own score interval
        for lvl, side, lo, hi in cal.level_intervals(res["thresholds"]):
            boots = [cal.interval_lr(s, y, lo, hi, side, w) for w in m.cluster_weights(clusters, reps)]
            hit = (s >= lo) & (s < hi) if side == "ge" else (s > lo) & (s <= hi)
            ci, how = lr_interval(boots, hit, y)
            evaluated["interval"][f"{lvl} {'pathogenic' if side == 'ge' else 'benign'}"] = {
                "interval": [lo, hi], "lr": cal.interval_lr(s, y, lo, hi, side), "ci": ci, "ci_method": how,
                "n_pos": int(hit[y == 1].sum()), "n_neg": int(hit[y == 0].sum()),
                "fraction_pos": float(hit[y == 1].mean()), "fraction_neg": float(hit[y == 0].mean())}
        # supplementary: everything beyond each threshold; the published single-level SpliceAI rules
        checks = [(lvl, "ge", p) for lvl, (p, _) in res["thresholds"].items() if p is not None]
        checks += [(lvl, "le", b) for lvl, (_, b) in res["thresholds"].items() if b is not None]
        if t == "SpliceAI":
            checks += [("published PP3 (>=0.2)", "ge", 0.2), ("published BP4 (<=0.1)", "le", 0.1)]
        for lvl, side, thr in checks:
            boots = [cal.likelihood_ratio(s, y, thr, side, w) for w in m.cluster_weights(clusters, reps)]
            ci, how = lr_interval(boots, (s >= thr) if side == "ge" else (s <= thr), y)
            evaluated["cumulative"][f"{lvl} {'pathogenic' if side == 'ge' else 'benign'}"] = {
                "threshold": thr, "lr": cal.likelihood_ratio(s, y, thr, side), "ci": ci, "ci_method": how,
                "fraction_pos_meeting": float(((s >= thr) if side == "ge" else (s <= thr))[y == 1].mean()),
                "fraction_neg_meeting": float(((s >= thr) if side == "ge" else (s <= thr))[y == 0].mean())}
        # supplementary: the same intervals within each region of the main set
        by_region = {}
        for rname, rmask in pdata.regions(v).items():
            if rname in ("intronic >=3 bp (main)", "canonical", "synonymous", "5' UTR", "3' UTR"):
                continue
            rv = v[test & v[t].notna() & rmask]
            rs, ry = rv[t].to_numpy(), rv["y"].to_numpy()
            rc = m.clusters_of(rv["gene_id"].tolist())
            rows = {}
            for lvl, side, lo, hi in cal.level_intervals(res["thresholds"]):
                hit = (rs >= lo) & (rs < hi) if side == "ge" else (rs > lo) & (rs <= hi)
                boots = [cal.interval_lr(rs, ry, lo, hi, side, w) for w in m.cluster_weights(rc, reps)]
                ci, how = lr_interval(boots, hit, ry)
                rows[f"{lvl} {'pathogenic' if side == 'ge' else 'benign'}"] = {
                    "lr": cal.interval_lr(rs, ry, lo, hi, side), "ci": ci, "ci_method": how,
                    "n_pos": int(hit[ry == 1].sum()), "n_neg": int(hit[ry == 0].sum())}
            by_region[rname] = {"n_pos": int(ry.sum()), "n_neg": int((1 - ry).sum()), "levels": rows}
        evaluated["interval_by_region"] = by_region
        # supplementary: thresholds for other priors (targets by the rule of Pejaver et al.)
        prior_sens = {}
        for prior in (0.01, 0.0441, 0.10):
            lv = cal.levels_for_prior(prior)
            th = cal.thresholds_from_curve(res["curve"], lv)
            prior_sens[str(prior)] = {"targets": lv, "thresholds": th,
                                      "fraction_pos_at_least_supporting": float((s >= th["supporting"][0])[y == 1].mean())
                                      if th["supporting"][0] is not None else None,
                                      "fraction_neg_benign_supporting": float((s <= th["supporting"][1])[y == 0].mean())
                                      if th["supporting"][1] is not None else None}
        evaluated["prior_sensitivity"] = prior_sens
        out[t] = {"derivation": {k: res[k] for k in ("n_pos", "n_neg", "thresholds")},
                  "curve": res["curve"], "test_n_pos": int(y.sum()), "test_n_neg": int((1 - y).sum()),
                  "evaluated_on_main": evaluated}
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--replicates", type=int, default=m.REPLICATES)
    ap.add_argument("--calibration-replicates", type=int, default=cal.BOOT)
    a = ap.parse_args(argv)
    v = pdata.load()
    tools = pdata.tools_present(v)
    set_masks, region_masks = pdata.sets(v), pdata.regions(v)
    OUT.mkdir(parents=True, exist_ok=True)
    coverage(v, set_masks, tools).to_csv(OUT / "coverage.tsv", sep="\t", index=False)
    disc, pairs = [], []
    for sname, smask in set_masks.items():
        for rname, rmask in region_masks.items():
            rows, diffs = discrimination(v, smask & rmask, tools, a.replicates)
            disc += [{"set": sname, "region": rname, **r} for r in rows]
            pairs += [{"set": sname, "region": rname, **d} for d in diffs]
    pd.DataFrame(disc).to_csv(OUT / "discrimination.tsv", sep="\t", index=False)
    pd.DataFrame(pairs).to_csv(OUT / "paired_differences.tsv", sep="\t", index=False)
    calib = calibration_tables(v, set_masks, tools, a.replicates, a.calibration_replicates)
    (OUT / "calibration.json").write_text(json.dumps(calib, indent=1, default=float) + "\n", encoding="utf-8")
    print("tables written to", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
