"""Numbers sheet for the manuscript: every number the text uses, read from the result tables of this repository at its
current commit. Usage (from the repository root): python manuscript/tools/numbers_sheet.py > NUMBERS.md"""
import ast
import json
import subprocess
import sys

import pandas as pd

R1, R2 = "results/part1/", "results/part2/"


def ci(x):
    if isinstance(x, str):
        x = ast.literal_eval(x.replace("inf", "1e999").replace("nan", "None"))
    if not isinstance(x, (list, tuple)) or len(x) != 2:
        return "[NA]"
    lo, hi = x
    return f"[{f3(lo)}, {f3(hi)}]"


def f3(v):
    if v is None:
        return "NA"
    v = float(v)
    if v != v:
        return "NA"
    if v in (float("inf"), float("-inf")):
        return "∞" if v > 0 else "−∞"
    a = abs(v)
    s = f"{v:.3f}" if a < 10 else f"{v:.1f}" if a < 100 else f"{v:.0f}"
    return s.replace("-", "−")


def pct(v):
    return f"{100 * float(v):.1f}%"


head = subprocess.run(["git", "describe", "--tags", "--always"], capture_output=True, text=True).stdout.strip()
sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
print(f"# Numbers sheet (generated from the result tables at {head}, {sha})\n")

print("## Part 1 flow (results/part1/flow.tsv)\n")
print(pd.read_csv(R1 + "flow.tsv", sep="\t").fillna("").to_string(index=False))

print("\n## Part 1 coverage, main set\n")
cov = pd.read_csv(R1 + "coverage.tsv", sep="\t")
print(cov[cov["set"] == "main"].to_string(index=False))

print("\n## Part 1 discrimination, main set (AUROC / AUPRC with gene-bootstrap 95% CI)\n")
d = pd.read_csv(R1 + "discrimination.tsv", sep="\t")
m = d[d["set"] == "main"]
for reg, g in m.groupby("region", sort=False):
    print(f"### {reg}")
    for _, r in g.iterrows():
        print(f"- {r['tool']}: n {int(r['n_pos'])} P/LP / {int(r['n_neg'])} B/LB; AUROC {f3(r['auroc'])} {ci(r['auroc_ci'])}; "
              f"AUPRC {f3(r['auprc'])} {ci(r['auprc_ci'])} (baseline {f3(r['pos_fraction'])})")
print("\nOther sets present:", sorted(set(d["set"])))

print("\n## Part 1 paired AUROC differences (all sets)\n")
p = pd.read_csv(R1 + "paired_differences.tsv", sep="\t")
cols = list(p.columns)
for (s, reg), g in p.groupby([cols[0], cols[1]], sort=False):
    items = "; ".join(f"{r[cols[2]]} {f3(r['difference'])} {ci(r['ci'])} (n {int(r[cols[3]])}/{int(r[cols[4]])})"
                      for _, r in g.iterrows())
    print(f"- [{s}] {reg}: {items}")

print("\n## Part 1 supplementary (AUPRC, sensitivity at fixed specificity), all rows\n")
sd = pd.read_csv(R1 + "supplementary_discrimination.tsv", sep="\t")
for _, r in sd.iterrows():
    print(f"- [{r['set']}] {r['region']} | {r['pair']} | {r['metric']}: {f3(r['value_a'])} vs {f3(r['value_b'])}, "
          f"diff {f3(r['difference'])} {ci(r['ci'])} (n {int(r['n_pos'])}/{int(r['n_neg'])})")

print("\n## Part 1 calibration (results/part1/calibration.json)\n")
c = json.load(open(R1 + "calibration.json", encoding="utf-8"))
for tool, t in c.items():
    der = t["derivation"]
    print(f"### {tool}: derivation {der['n_pos']} P/LP / {der['n_neg']} B/LB; test {t.get('test_n_pos')} / {t.get('test_n_neg')}")
    print("- thresholds (pathogenic, benign): " + "; ".join(f"{k} {f3(v[0])}, {f3(v[1])}" for k, v in der["thresholds"].items()))
    ev = t.get("evaluated_on_main") or {}
    for kind in ("interval", "cumulative"):
        for lvl, v in (ev.get(kind) or {}).items():
            extra = (f"n {v['n_pos']}/{v['n_neg']}, share P/LP {pct(v['fraction_pos'])}, share B/LB {pct(v['fraction_neg'])}"
                     if kind == "interval" else
                     f"P/LP meeting {pct(v['fraction_pos_meeting'])}, B/LB meeting {pct(v['fraction_neg_meeting'])}")
            print(f"  - {kind} {lvl}: LR {f3(v['lr'])} {ci(v['ci'])} ({v.get('ci_method')}); {extra}")
    for reg, rv in (ev.get("interval_by_region") or {}).items():
        lv = "; ".join(f"{k} {f3(x['lr'])} {ci(x['ci'])} ({x['n_pos']}/{x['n_neg']})" for k, x in rv["levels"].items())
        print(f"  - by region {reg} (n {rv['n_pos']}/{rv['n_neg']}): {lv}")
    ps = ev.get("prior_sensitivity")
    if ps:
        print(f"  - prior sensitivity: {json.dumps(ps, default=str)[:1500]}")

print("\n## Published thresholds on the main set\n")
pt = pd.read_csv(R1 + "published_thresholds.tsv", sep="\t")
for _, r in pt.iterrows():
    print(f"- {r['tool']} {r['level']} ({r['side']} {r['lo']}..{r['hi']}): LR {f3(r['lr'])} {ci(r['ci'])}; "
          f"n {r['n_pos']}/{r['n_neg']}; share P/LP {pct(r['fraction_pos'])}, B/LB {pct(r['fraction_neg'])}")

print("\n## Part 2 per-dataset discrimination\n")
d2 = pd.read_csv(R2 + "discrimination.tsv", sep="\t")
for ds, g in d2.groupby("dataset", sort=False):
    print(f"### {ds} (n {int(g['n_pos'].iloc[0])} splice-altering / {int(g['n_neg'].iloc[0])} other)")
    for _, r in g.iterrows():
        print(f"- {r['tool']}: n {int(r['n_pos'])}/{int(r['n_neg'])}; AUROC {f3(r['auroc'])} {ci(r['auroc_ci'])}; "
              f"AUPRC {f3(r['auprc'])}; Spearman |effect| {f3(r['spearman_abs_effect'])}")

print("\n## Part 2 paired differences per dataset\n")
p2 = pd.read_csv(R2 + "paired_differences.tsv", sep="\t")
for ds, g in p2.groupby("dataset", sort=False):
    print(f"- {ds}: " + "; ".join(f"{r['pair']} {f3(r['difference'])} {ci(r['ci'])}" for _, r in g.iterrows()))

print("\n## Part 2 pooled (DerSimonian-Laird; Hartung-Knapp; leave-one-out)\n")
pj = json.load(open(R2 + "pooled.json", encoding="utf-8"))
for pair, v in pj["pooled"].items():
    hk = v.get("ci_hartung_knapp")
    print(f"- {pair}: {f3(v['estimate'])} {ci(v['ci'])}; HK {ci(hk) if hk else 'NA'}; I² {pct(v['i2']) if v.get('i2') is not None else 'NA'}; "
          f"tau² {f3(v['tau2'])}; Q {f3(v['q'])}; k {v['k']}; datasets {v['datasets']}")
    for name, lo in (v.get("leave_one_out") or {}).items():
        print(f"    - leave out {name}: {f3(lo.get('estimate'))} {ci(lo['ci']) if lo.get('ci') else ''} "
              f"HK {ci(lo['ci_hartung_knapp']) if lo.get('ci_hartung_knapp') else 'NA'}")
print("\n- counts: " + json.dumps(pj["counts"]))

print("\n## Part 2 by position (additional analysis)\n")
bp = json.load(open(R2 + "by_position_pooled.json", encoding="utf-8"))
for cls, pairs in bp.items():
    for pair, v in pairs.items():
        if isinstance(v, dict) and "estimate" in v:
            hk = v.get("ci_hartung_knapp")
            print(f"- {cls} | {pair}: {f3(v['estimate'])} {ci(v['ci'])}; HK {ci(hk) if hk else 'NA'}; k {v['k']}; {v['datasets']}")
bd = pd.read_csv(R2 + "by_position_differences.tsv", sep="\t")
print("\nPer dataset, AVI − SpliceAI by class:")
for _, r in bd[bd["pair"] == "AVI - SpliceAI"].iterrows():
    print(f"- {r['class']} | {r['dataset']}: {f3(r['difference'])} {ci(r['ci'])} (n {r['n_pos']}/{r['n_neg']})")

print("\n## Part 2 variant positions (analysis sets)\n")
vp = pd.read_csv(R2 + "variant_positions.tsv", sep="\t")
print(vp.pivot_table(index="dataset", columns="class", values="variants", aggfunc="sum").fillna(0).astype(int).to_string())
print()
print(vp.pivot_table(index="dataset", columns="class", values="splice_altering", aggfunc="sum").fillna(0).astype(int).to_string())

print("\n## Derived totals and shares (arithmetic on the tables above)\n")
counts = pj["counts"]
main_ds = ["kitzman_pou1f1_ex2", "kitzman_wt1_ex9", "braun2018_mst1r_ex11", "cortes_lopez2022_cd19", "compass", "opensplice"]
print(f"- Part 2 main datasets: {sum(counts[k]['variants'] for k in main_ds):,} variants, "
      f"{sum(counts[k]['splice_altering'] for k in main_ds):,} splice-altering")
fl = pd.read_csv(R1 + "flow.tsv", sep="\t")
for step in ("new-label set", "main set"):
    r = fl[fl["step"] == step].iloc[0]
    print(f"- {step}: {int(r['total']):,} ({int(r['P/LP']):,} P/LP, {int(r['B/LB']):,} B/LB)")
cv = cov[cov["set"] == "main"].set_index("tool")
print(f"- SpliceAI coverage of the main set: {cv.loc['SpliceAI', 'scored']:,}/{cv.loc['SpliceAI', 'variants']:,} = "
      f"{100 * cv.loc['SpliceAI', 'scored'] / cv.loc['SpliceAI', 'variants']:.1f}%")
for tool in ("AVI", "SpliceAI", "SpliceAI2", "CADD"):
    iv = c[tool]["evaluated_on_main"]["interval"]
    n = c[tool]["test_n_pos"]
    path = sum(v["n_pos"] for k, v in iv.items() if k.endswith("pathogenic"))
    strong = iv.get("strong pathogenic", {}).get("n_pos", 0)
    benign = sum(v["fraction_neg"] for k, v in iv.items() if k.endswith("benign"))
    print(f"- {tool}: any pathogenic evidence {path}/{n} = {100 * path / n:.2f}%; strong {strong}/{n} = "
          f"{100 * strong / n:.2f}%; B/LB with benign evidence {100 * benign:.2f}%")

print("\n## Calibration sample composition (variant list exclusion reason; no scores)\n")
import csv as _csv
import gzip as _gzip
from collections import Counter as _Counter
_lab = {}
with open("data/derived/variants/calibration_sample.tsv", encoding="utf-8", newline="") as _h:
    for _r in _csv.DictReader(_h, delimiter="\t"):
        _lab[_r["variation_id"]] = _r["label"]
_cnt = _Counter()
with _gzip.open("data/derived/variants/clinvar_variants.tsv.gz", "rt", encoding="utf-8", newline="") as _h:
    for _r in _csv.DictReader(_h, delimiter="\t"):
        if _r["variation_id"] in _lab:
            _cnt[_r["exclusion"]] += 1
for _k, _n in sorted(_cnt.items()):
    print(f"- {_k or '(none)'}: {_n:,}")
print(f"- total: {sum(_cnt.values()):,} of {len(_lab):,} sample rows")

print("\n## Part 2 COMPASS by source subset (additional analysis added after a code-consistency review; defined in 582cb74)\n")
import os as _os
if _os.path.exists(R2 + "compass_subsets_pooled.json"):
    _cd = pd.read_csv(R2 + "compass_subsets_discrimination.tsv", sep="\t")
    for _, r in _cd.iterrows():
        print(f"- [{r['dataset']}] {r['tool']}: n {int(r['n_pos'])}/{int(r['n_neg'])}; AUROC {f3(r['auroc'])} {ci(r['auroc_ci'])}; "
              f"AUPRC {f3(r['auprc'])} {ci(r['auprc_ci'])}")
    _cp = pd.read_csv(R2 + "compass_subsets_paired_differences.tsv", sep="\t")
    for _, r in _cp.iterrows():
        print(f"- [{r['dataset']}] {r['pair']}: {f3(r['difference'])} {ci(r['ci'])} (n {int(r['n_pos'])}/{int(r['n_neg'])})")
    _cj = json.load(open(R2 + "compass_subsets_pooled.json", encoding="utf-8"))
    for pair, v in _cj["pooled_with_random_compass"].items():
        hk = v.get("ci_hartung_knapp")
        print(f"- pooled with random-SNV COMPASS | {pair}: {f3(v['estimate'])} {ci(v['ci'])}; HK {ci(hk) if hk else 'NA'}; "
              f"I² {pct(v['i2'])}; k {v['k']}; {v['datasets']}")
    print("- counts: " + json.dumps(_cj["counts"]))

print("\n## COMPASS random-variant proxy vs archived ClinVar releases (labels only; 4b83fb8)\n")
if _os.path.exists(R2 + "compass_clinvar_check.json"):
    _cc = json.load(open(R2 + "compass_clinvar_check.json", encoding="utf-8"))
    print("- COMPASS main SNVs: " + json.dumps(_cc["compass_main_snvs"]))
    for _rel, _cnt in _cc["present_in"].items():
        print(f"- present in {_rel}: " + json.dumps(_cnt))
    if "analysed_snvs" in _cc:
        print("- analysed COMPASS SNVs: " + json.dumps(_cc["analysed_snvs"]))
        for _rel, _cnt in _cc["analysed_present_in"].items():
            print(f"- analysed, present in {_rel}: " + json.dumps(_cnt))
    if "analysed_random_significance_2026_10" in _cc:
        print("- analysed random alleles in October 2026, ClinVar classification: "
              + json.dumps(_cc["analysed_random_significance_2026_10"]))
