"""Part 2 analysis (STUDY_DESIGN.md section 5.3): per-dataset discrimination of the authors' splice-altering
calls by each tool, rank correlation with the measured effect, paired AUROC differences and their
random-effects pooling. Writes tables to results/part2/.

Usage:
    python src/part2_analysis.py [--replicates 2000]
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import json
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

import metrics as m
from paths import DERIVED, EXPOSURE, ROOT

TABLE = DERIVED / "part2" / "assay_variants.tsv.gz"
SCORES = DERIVED / "scores"
OUT = ROOT / "results" / "part2"
TOOLS = {"AVI": ("part2_avi.tsv.gz", "avi_raw"), "SpliceAI": ("part2_spliceai.tsv.gz", "spliceai_max"),
         "SpliceAI masked": ("part2_spliceai_masked.tsv.gz", "spliceai_max"),
         "CADD": ("part2_cadd.tsv.gz", "cadd_phred"), "phyloP 100-way": ("part2_conservation.tsv.gz", "phyloP100"),
         "Pangolin": ("part2_pangolin.tsv.gz", "pangolin_max"),
         "SpliceAI2": ("part2_spliceai2.tsv.gz", "spliceai2_summary")}
SINGLE_GENE = {"kitzman_pou1f1_ex2", "kitzman_wt1_ex9", "braun2018_mst1r_ex11", "cortes_lopez2022_cd19"}
# Datasets a tool's developers evaluated it on: reported per dataset, left out of that tool's pooled estimate.
DEVELOPER_EVALUATED = {"SpliceAI2": {"opensplice"}}


def _read(path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as h:
        return pd.read_csv(h, sep="\t", dtype=str, keep_default_na=False)


def development_hits(alleles: pd.DataFrame, path=EXPOSURE / "development_positions.tsv") -> set:
    """Allele keys whose position lies in a developers' benchmark interval (streamed, low memory)."""
    by_chrom = defaultdict(list)
    for k, c, p in zip(alleles["allele"], alleles["chrom"], alleles["pos"].astype(int)):
        by_chrom[c].append((p, k))
    for c in by_chrom:
        by_chrom[c].sort()
    starts = {c: [p for p, _ in v] for c, v in by_chrom.items()}
    hit = set()
    with open(path, encoding="utf-8") as h:
        next(h)
        for line in h:
            chrom, s, e, _ = line.rstrip("\n").split("\t")
            if chrom not in starts:
                continue
            i, j = bisect.bisect_left(starts[chrom], int(s)), bisect.bisect_right(starts[chrom], int(e))
            hit.update(k for _, k in by_chrom[chrom][i:j])
    return hit


def part1_new_labels() -> set:
    v = _read(DERIVED / "variants" / "clinvar_variants.tsv.gz")
    v = v[v["new_label"] == "yes"]
    return set(v["chrom"] + ":" + v["pos"] + ":" + v["ref"] + ":" + v["alt"])


def load(dev_path=EXPOSURE / "development_positions.tsv") -> pd.DataFrame:
    d = _read(TABLE)
    d = d[d["use"].isin(["main", "comparison"]) & (d["variant_type"] == "SNV")
          & d["splice_altering"].isin(["yes", "no"])].copy()
    d["allele"] = d["chrom"] + ":" + d["pos"] + ":" + d["ref"] + ":" + d["alt"]
    d["y"] = (d["splice_altering"] == "yes").astype(int)
    d["effect"] = pd.to_numeric(d["effect_size"], errors="coerce").abs()
    d["in_part1"] = d["allele"].isin(part1_new_labels())
    d["development"] = d["allele"].isin(development_hits(d.drop_duplicates("allele"), dev_path))
    for tool, (name, col) in TOOLS.items():
        if (SCORES / name).exists():
            s = _read(SCORES / name)[["allele", col]].rename(columns={col: tool})
            s[tool] = pd.to_numeric(s[tool], errors="coerce")
            d = d.merge(s, on="allele", how="left")
    return d.reset_index(drop=True)


def analysis_sets(d: pd.DataFrame) -> dict:
    """{name: mask}: each main dataset after exclusions; the developers'-gene rows; MFASS (comparison)."""
    clean = ~d["in_part1"] & ~d["development"]
    out = {}
    for ds in sorted(d["dataset"].unique()):
        use = d["use"] == ("comparison" if ds == "mfass" else "main")
        out[ds] = (d["dataset"] == ds) & use & clean
    out["opensplice, developers' benchmark positions"] = (d["dataset"] == "opensplice") & d["development"]
    return out


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    ok = ~np.isnan(x) & ~np.isnan(y)
    if ok.sum() < 3:
        return float("nan")
    return float(pd.Series(x[ok]).rank().corr(pd.Series(y[ok]).rank()))


def dataset_tables(d: pd.DataFrame, mask: pd.Series, name: str, tools: list, reps: int) -> tuple[list, list]:
    sub = d[mask]
    clusters_by = (lambda x: np.arange(len(x)).astype(str)) if name in SINGLE_GENE or name.startswith("kitzman") \
        else (lambda x: m.clusters_of(x["gene"].tolist()))
    rows, diffs = [], {}
    for t in tools:
        x = sub[sub[t].notna()]
        y, s = x["y"].to_numpy(), x[t].to_numpy()
        if y.sum() == 0 or (1 - y).sum() == 0:
            rows.append({"dataset": name, "tool": t, "n_pos": int(y.sum()), "n_neg": int((1 - y).sum())})
            continue
        prep = m._order(s)
        b_roc, b_pr = [], []
        for w in m.cluster_weights(clusters_by(x), reps):
            b_roc.append(m.auroc(s, y, w, prep))
            b_pr.append(m.auprc(s, y, w, prep))
        rows.append({"dataset": name, "tool": t, "n_pos": int(y.sum()), "n_neg": int((1 - y).sum()),
                     "pos_fraction": float(y.mean()), "auroc": m.auroc(s, y, prepared=prep),
                     "auroc_ci": m.interval(b_roc), "auprc": m.auprc(s, y, prepared=prep),
                     "auprc_ci": m.interval(b_pr), "spearman_abs_effect": spearman(s, x["effect"].to_numpy())})
    for t in [t for t in tools if t != "AVI"]:
        x = sub[sub["AVI"].notna() & sub[t].notna()]
        y = x["y"].to_numpy()
        if y.sum() == 0 or (1 - y).sum() == 0:
            continue
        sa, sb = x["AVI"].to_numpy(), x[t].to_numpy()
        pa, pb = m._order(sa), m._order(sb)
        boots = [m.auroc(sa, y, w, pa) - m.auroc(sb, y, w, pb) for w in m.cluster_weights(clusters_by(x), reps)]
        diffs[t] = {"dataset": name, "pair": f"AVI - {t}", "n_pos": int(y.sum()), "n_neg": int((1 - y).sum()),
                    "difference": m.auroc(sa, y, prepared=pa) - m.auroc(sb, y, prepared=pb),
                    "ci": m.interval(boots), "se": float(np.nanstd(boots, ddof=1))}
    return rows, list(diffs.values())


def dersimonian_laird(est: list, se: list) -> dict:
    """Random-effects pooled estimate with 95% interval, tau^2, Q and I^2; plus the Hartung-Knapp interval
    (t distribution with k-1 degrees of freedom and the weighted residual variance), which is wider and
    better behaved when the number of datasets is small."""
    est, se = np.asarray(est, float), np.asarray(se, float)
    k = len(est)
    w = 1 / se ** 2
    fixed = (w * est).sum() / w.sum()
    q = (w * (est - fixed) ** 2).sum()
    tau2 = max(0.0, (q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0.0
    wr = 1 / (se ** 2 + tau2)
    pooled = (wr * est).sum() / wr.sum()
    se_p = (1 / wr.sum()) ** 0.5
    out = {"estimate": float(pooled), "ci": [float(pooled - 1.96 * se_p), float(pooled + 1.96 * se_p)],
           "tau2": float(tau2), "q": float(q), "i2": float(max(0.0, (q - (k - 1)) / q)) if q > 0 else 0.0,
           "k": int(k)}
    if k > 1:
        half = stats.t.ppf(0.975, k - 1) * ((wr * (est - pooled) ** 2).sum() / (k - 1) / wr.sum()) ** 0.5
        out["ci_hartung_knapp"] = [float(pooled - half), float(pooled + half)]
    return out


def pooled_differences(diffs: list, main_sets: list, tools: list) -> dict:
    """Random-effects pooled AVI − tool AUROC differences over the main datasets, leaving out the datasets
    the tool's developers evaluated it on."""
    pooled = {}
    for t in [t for t in tools if t != "AVI"]:
        keep = [s for s in main_sets if s not in DEVELOPER_EVALUATED.get(t, set())]
        x = [r for r in diffs if r["pair"] == f"AVI - {t}" and r["dataset"] in keep and r["se"] > 0]
        if x:
            pooled[f"AVI - {t}"] = {**dersimonian_laird([r["difference"] for r in x], [r["se"] for r in x]),
                                    "datasets": [r["dataset"] for r in x]}
    return pooled


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--replicates", type=int, default=m.REPLICATES)
    a = ap.parse_args(argv)
    d = load()
    tools = [t for t in TOOLS if t in d.columns]
    OUT.mkdir(parents=True, exist_ok=True)
    rows, diffs = [], []
    sets = analysis_sets(d)
    for name, mask in sets.items():
        r, df = dataset_tables(d, mask, name, tools, a.replicates)
        rows += r
        diffs += df
    pd.DataFrame(rows).to_csv(OUT / "discrimination.tsv", sep="\t", index=False)
    pd.DataFrame(diffs).to_csv(OUT / "paired_differences.tsv", sep="\t", index=False)
    main_sets = [n for n in sets if n != "mfass" and not n.startswith("opensplice, ")]
    pooled = pooled_differences(diffs, main_sets, tools)
    counts = {n: {"variants": int(mk.sum()), "splice_altering": int(d.loc[mk, "y"].sum())} for n, mk in sets.items()}
    (OUT / "pooled.json").write_text(json.dumps({"pooled": pooled, "counts": counts}, indent=1) + "\n",
                                     encoding="utf-8")
    print("tables written to", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
