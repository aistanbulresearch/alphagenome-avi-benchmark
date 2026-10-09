"""Manuscript figures from results/part1 and results/part2 (and, for the ROC curves, the Part 1 score table).
Static print figures, PDF and PNG, in manuscript/figures/.

    Figure 1  Part 1 discrimination: ROC curves (main region) and AVI − tool AUROC differences by region
    Figure 2  Part 1 calibration: interval likelihood ratios against the ClinGen targets; share of variants
              receiving evidence
    Figure 3  Part 2 forest plots: AVI − tool AUROC differences per dataset, pooled (DerSimonian-Laird and
              Hartung-Knapp)
    Figure 4  Part 2 by position relative to the exon: pooled AVI − tool AUROC differences
    Figure S1 SpliceAI2 (added comparator): Part 1 by region and Part 2 per dataset
    Figure S2 Published missense CADD thresholds evaluated on intronic variants

Every tool keeps one colour in every figure (palette checked for colour-vision deficiency separation).

Usage:
    python src/figures.py [--only 1 2 ...]
"""
from __future__ import annotations

import argparse
import ast
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from paths import ROOT  # noqa: E402

R1, R2 = ROOT / "results" / "part1", ROOT / "results" / "part2"
OUT = ROOT / "manuscript" / "figures"
MAIN = "intronic >=3 bp (main)"
COLOR = {"AVI": "#2a78d6", "SpliceAI": "#eb6834", "SpliceAI masked": "#eb6834", "Pangolin": "#1baf7a",
         "CADD": "#eda100", "phyloP 100-way": "#e87ba4", "SpliceAI2": "#4a3aa7"}
INK, INK2, MUTED, GRID, BASE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
TARGETS = {"supporting": 2.406, "moderate": 5.790, "strong": 33.53, "very_strong": 1124.0}
REGIONS = [("donor region", "Donor region"), ("acceptor region", "Acceptor region"),
           ("near intronic", "Near intronic"), ("deep intronic", "Deep intronic"), (MAIN, "All intronic ≥3 bp")]
DATASETS = [("braun2018_mst1r_ex11", "MST1R exon 11"), ("compass", "COMPASS"), ("cortes_lopez2022_cd19", "CD19"),
            ("kitzman_pou1f1_ex2", "POU1F1 exon 2"), ("kitzman_wt1_ex9", "WT1 exon 9"), ("opensplice", "OpenSplice")]
COMPARISON = [("mfass", "MFASS"), ("opensplice, developers' benchmark positions", "OpenSplice, developers' positions")]
MM = 1 / 25.4


def style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7,
        "axes.edgecolor": BASE, "axes.linewidth": 0.6, "axes.labelcolor": INK2, "axes.titlesize": 7.5,
        "axes.titleweight": "bold", "axes.titlecolor": INK, "axes.titlelocation": "left",
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.grid": False, "grid.color": GRID, "grid.linewidth": 0.5, "legend.frameon": False,
        "legend.fontsize": 6.5, "axes.spines.top": False, "axes.spines.right": False,
        "savefig.dpi": 300, "pdf.fonttype": 42, "figure.facecolor": "white", "axes.facecolor": "white"})


def ci(text) -> tuple:
    if isinstance(text, (list, tuple)):
        return tuple(float(x) for x in text)
    return tuple(float(x) for x in ast.literal_eval(str(text).replace("nan", "None").replace("inf", "1e308"))
                 if x is not None) or (np.nan, np.nan)


def save(fig, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("written", OUT / f"{name}.pdf")


def roc_points(s: np.ndarray, y: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Weighted ROC curve: false and true positive rates at every distinct threshold, descending."""
    order = np.argsort(-s, kind="mergesort")
    s, y, w = s[order], y[order], w[order]
    last = np.r_[np.nonzero(np.diff(s))[0], len(s) - 1]
    tp, fp = np.cumsum(w * (y == 1))[last], np.cumsum(w * (y == 0))[last]
    return np.r_[0, fp / fp[-1]], np.r_[0, tp / tp[-1]]


def diff_rows(ax, rows: list, tools: list, offsets=None) -> None:
    """Dot-and-interval rows: rows = [(label, {tool: (estimate, lo, hi)})], top to bottom."""
    offsets = offsets or {t: (len(tools) - 1) / 2 * 0.22 - i * 0.22 for i, t in enumerate(tools)}
    for r, (_, vals) in enumerate(rows):
        yy = len(rows) - 1 - r
        for t in tools:
            if t not in vals:
                continue
            est, lo, hi = vals[t]
            hollow = t == "SpliceAI masked"
            ax.plot([lo, hi], [yy + offsets[t]] * 2, color=COLOR[t], lw=1.0, solid_capstyle="round")
            ax.plot(est, yy + offsets[t], "o", ms=3.6, mfc="white" if hollow else COLOR[t], mec=COLOR[t],
                    mew=1.0 if hollow else 0.6 if not hollow else 1.0)
    ax.axvline(0, color=BASE, lw=0.7, zorder=0)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([lab for lab, _ in rows][::-1])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_axisbelow(True)


def figure1() -> None:
    import part1_analysis as pa
    import part1_data as pdata
    v = pdata.load()
    mask = pdata.sets(v)["main"] & pdata.regions(v)[MAIN]
    disc = pd.read_csv(R1 / "discrimination.tsv", sep="\t")
    pairs = pd.read_csv(R1 / "paired_differences.tsv", sep="\t")
    fig, (a, b) = plt.subplots(1, 2, figsize=(180 * MM, 72 * MM), gridspec_kw={"width_ratios": [1, 1.25]})
    tools = ["AVI", "SpliceAI", "Pangolin", "CADD", "phyloP 100-way"]
    for t in tools:
        d = v[mask & v[t].notna()]
        fpr, tpr = roc_points(d[t].to_numpy(), d["y"].to_numpy(), pa.base_weights(d, (t,)))
        auc = disc[(disc.set == "main") & (disc.region == MAIN) & (disc.tool == t)]["auroc"].iloc[0]
        a.plot(fpr, tpr, color=COLOR[t], lw=1.1, label=f"{t.replace(' 100-way', '')}  {auc:.3f}",
               zorder=3 if t == "AVI" else 2)
    a.set_xscale("log")
    a.set_xlim(1e-3, 1)
    a.set_ylim(0, 1.005)
    a.set_xlabel("False positive rate (log scale)")
    a.set_ylabel("True positive rate")
    a.grid(color=GRID, lw=0.5)
    a.set_axisbelow(True)
    a.legend(title="AUROC", title_fontsize=6.5, loc="lower right", handlelength=1.4)
    a.set_title("a  ROC, intronic ≥3 bp (449 P/LP, 28,188 B/LB)")
    rows = []
    for key, label in REGIONS:
        vals = {}
        for t in ("SpliceAI", "Pangolin", "CADD"):
            r = pairs[(pairs.set == "main") & (pairs.region == key) & (pairs.pair == f"AVI - {t}")]
            if len(r):
                lo, hi = ci(r["ci"].iloc[0])
                vals[t] = (r["difference"].iloc[0], lo, hi)
        npos = int(disc[(disc.set == "main") & (disc.region == key) & (disc.tool == "AVI")]["n_pos"].iloc[0])
        rows.append((f"{label} ({npos} P/LP)", vals))
    diff_rows(b, rows, ["SpliceAI", "Pangolin", "CADD"])
    b.set_xlabel("AUROC difference, AVI − tool (95% gene-bootstrap interval)")
    b.set_title("b  AVI − tool, by region")
    handles = [plt.Line2D([], [], color=COLOR[t], marker="o", ms=3.6, lw=1.0, label=f"AVI − {t}")
               for t in ("SpliceAI", "Pangolin", "CADD")]
    b.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -0.17), ncol=3, columnspacing=1.2,
             handlelength=1.4)
    fig.tight_layout(w_pad=2.5)
    save(fig, "figure1_part1_discrimination")


LEVELS = [("supporting pathogenic", "Supporting", "supporting", False),
          ("moderate pathogenic", "Moderate", "moderate", False),
          ("strong pathogenic", "Strong", "strong", False),
          ("very_strong pathogenic", "Very strong", "very_strong", False),
          ("supporting benign", "Supporting (benign)", "supporting", True),
          ("moderate benign", "Moderate (benign)", "moderate", True)]


def calibration_panel(ax, cal: dict, tools: list) -> None:
    n = len(LEVELS)
    offs = {t: (len(tools) - 1) / 2 * 0.22 - i * 0.22 for i, t in enumerate(tools)}
    lo_lim, hi_lim = 0.01, 3000
    for r, (key, label, lvl, benign) in enumerate(LEVELS):
        yy = n - 1 - r
        target = 1 / TARGETS[lvl] if benign else TARGETS[lvl]
        ax.plot([target, target], [yy - 0.38, yy + 0.38], color=INK2, lw=1.2, solid_capstyle="butt", zorder=1)
        for t in tools:
            e = cal[t]["evaluated_on_main"]["interval"].get(key)
            if not e or e["lr"] is None or np.isnan(e["lr"]) or e["lr"] <= 0:
                continue
            if np.isinf(e["lr"]):                       # P/LP in the interval, no B/LB: drawn at the edge
                ax.annotate("", xy=(hi_lim * 1.25, yy + offs[t]), xytext=(hi_lim / 3, yy + offs[t]),
                            arrowprops={"arrowstyle": "->", "color": COLOR[t], "lw": 1.0}, annotation_clip=False)
                ax.text(hi_lim * 1.25, yy + offs[t] - 0.07, f"{e['n_pos']} P/LP, 0 B/LB", ha="right", va="top",
                        fontsize=5.5, color=INK2, zorder=4,
                        bbox={"boxstyle": "square,pad=0.15", "fc": "white", "ec": "none"})
                continue
            lo, hi = e["ci"]
            lo = max(lo, lo_lim) if np.isfinite(lo) and lo > 0 else lo_lim
            hi = min(hi, hi_lim) if np.isfinite(hi) else hi_lim
            ax.plot([lo, hi], [yy + offs[t]] * 2, color=COLOR[t], lw=1.0, solid_capstyle="round", zorder=2)
            ax.plot(e["lr"], yy + offs[t], "o", ms=3.6, color=COLOR[t], mec="white", mew=0.4, zorder=3)
            if not np.isfinite(e["ci"][1]) or e["ci"][1] > hi_lim:
                ax.annotate("", xy=(hi_lim * 1.25, yy + offs[t]), xytext=(hi_lim, yy + offs[t]),
                            arrowprops={"arrowstyle": "->", "color": COLOR[t], "lw": 0.8}, annotation_clip=False)
    ax.set_xscale("log")
    ax.set_xlim(lo_lim, hi_lim)
    ax.axvline(1, color=BASE, lw=0.7, zorder=0)
    ax.axhline(1.5, color=GRID, lw=0.6, zorder=0)
    ax.set_yticks(range(n))
    ax.set_yticklabels([lab for _, lab, _, _ in LEVELS][::-1])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.set_xlabel("Likelihood ratio of the score interval on the main set (log scale)")


def evidence_share(cal: dict, tool: str) -> tuple[float, float]:
    iv = cal[tool]["evaluated_on_main"]["interval"]
    pos = sum(e["fraction_pos"] for k, e in iv.items() if k.endswith("pathogenic"))
    neg = sum(e["fraction_neg"] for k, e in iv.items() if k.endswith("benign"))
    return pos, neg


def figure2() -> None:
    cal = json.load(open(R1 / "calibration.json", encoding="utf-8"))
    tools = ["AVI", "SpliceAI", "CADD"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(180 * MM, 66 * MM), gridspec_kw={"width_ratios": [1.6, 1]})
    calibration_panel(a, cal, tools)
    handles = [plt.Line2D([], [], color=COLOR[t], marker="o", ms=3.6, lw=1.0, label=t) for t in tools]
    handles.append(plt.Line2D([], [], color=INK2, lw=1.2, marker="|", ms=7, ls="none", label="ClinGen target"))
    a.legend(handles=handles, loc="upper left", ncol=1, handlelength=1.4, borderaxespad=0.2)
    a.set_title("a  Evidence levels, each tested within its own score interval")
    rows = [t for t in tools]
    for i, t in enumerate(rows):
        pos, neg = evidence_share(cal, t)
        yy = len(rows) - 1 - i
        b.plot([0, pos * 100], [yy + 0.13] * 2, color=COLOR[t], lw=1.0)
        b.plot(pos * 100, yy + 0.13, "o", ms=4, color=COLOR[t])
        b.plot([0, neg * 100], [yy - 0.13] * 2, color=COLOR[t], lw=1.0)
        b.plot(neg * 100, yy - 0.13, "o", ms=4, mfc="white", mec=COLOR[t], mew=1.0)
        b.text(pos * 100 + 3, yy + 0.13, f"{pos:.0%}", va="center", fontsize=6, color=INK2)
        b.text(neg * 100 + 3, yy - 0.13, f"{neg:.0%}" if neg > 0 else "no benign level reached", va="center",
               fontsize=6, color=INK2)
    b.set_yticks(range(len(rows)))
    b.set_yticklabels(rows[::-1])
    b.tick_params(axis="y", length=0)
    b.spines["left"].set_visible(False)
    b.set_xlim(0, 112)
    b.set_xticks([0, 25, 50, 75, 100])
    b.grid(axis="x", color=GRID, lw=0.5)
    b.set_axisbelow(True)
    b.set_xlabel("Share of variants receiving evidence (%)")
    hb = [plt.Line2D([], [], color=INK2, marker="o", ms=4, ls="none", label="P/LP with pathogenic evidence"),
          plt.Line2D([], [], color=INK2, marker="o", ms=4, mfc="white", ls="none",
                     label="B/LB with benign evidence")]
    b.legend(handles=hb, loc="upper left", bbox_to_anchor=(0, -0.2), ncol=1, handlelength=1.0)
    b.set_title("b  Evidence coverage (prior 0.0441)")
    fig.tight_layout(w_pad=2.5)
    save(fig, "figure2_part1_calibration")


def forest(ax, pairs: pd.DataFrame, pooled: dict, tool: str, datasets: list, comparison: list,
           title: str) -> None:
    rows = []
    for key, label in datasets:
        r = pairs[(pairs.dataset == key) & (pairs.pair == f"AVI - {tool}")]
        rows.append((label, None if r.empty else (r["difference"].iloc[0], *ci(r["ci"].iloc[0])), "main"))
    p = pooled.get(f"AVI - {tool}")
    rows.append(("Pooled", p, "pooled"))
    for key, label in comparison:
        r = pairs[(pairs.dataset == key) & (pairs.pair == f"AVI - {tool}")]
        rows.append((label, None if r.empty else (r["difference"].iloc[0], *ci(r["ci"].iloc[0])), "comparison"))
    n = len(rows)
    col = COLOR[tool]
    hollow = tool == "SpliceAI masked"
    for i, (label, val, kind) in enumerate(rows):
        yy = n - 1 - i - (0.6 if kind == "comparison" else 0)
        if val is None:
            ax.text(0, yy, "not scored" if kind != "pooled" else "", ha="center", va="center", fontsize=5.5,
                    color=MUTED)
            continue
        if kind == "pooled":
            est, lo, hi = val["estimate"], *val["ci"]
            hk = val.get("ci_hartung_knapp")
            if hk:
                ax.plot(hk, [yy] * 2, color=INK, lw=0.7)
            ax.fill([lo, est, hi, est], [yy, yy + 0.32, yy, yy - 0.32], color=INK, lw=0)
            continue
        est, lo, hi = val
        c = MUTED if kind == "comparison" else col
        ax.plot([lo, hi], [yy] * 2, color=c, lw=1.0, solid_capstyle="round")
        ax.plot(est, yy, "s" if kind == "main" else "o", ms=3.4, mfc="white" if hollow and kind == "main" else c,
                mec=c, mew=1.0 if hollow else 0.5)
    ax.axvline(0, color=BASE, lw=0.7, zorder=0)
    ax.axhline(n - 1 - len(datasets) + 0.5, color=GRID, lw=0.6)
    ax.axhline(n - 1 - len(datasets) - 0.5 - 0.3, color=GRID, lw=0.6)
    ax.set_yticks([n - 1 - i - (0.6 if k == "comparison" else 0) for i, (_, _, k) in enumerate(rows)])
    ax.set_yticklabels([lab for lab, _, _ in rows])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.set_title(title)


def figure3() -> None:
    pairs = pd.read_csv(R2 / "paired_differences.tsv", sep="\t")
    pooled = json.load(open(R2 / "pooled.json", encoding="utf-8"))["pooled"]
    tools = [("SpliceAI", "a  AVI − SpliceAI (raw)"), ("SpliceAI masked", "b  AVI − SpliceAI (masked)"),
             ("Pangolin", "c  AVI − Pangolin"), ("CADD", "d  AVI − CADD")]
    fig, axes = plt.subplots(1, 4, figsize=(180 * MM, 70 * MM), sharey=True, sharex=True)
    for ax, (t, title) in zip(axes, tools):
        forest(ax, pairs, pooled, t, DATASETS, COMPARISON, title)
    for ax in axes[1:]:
        ax.tick_params(axis="y", labelleft=False)
    axes[0].set_xlim(-0.4, 0.25)
    fig.supxlabel("AUROC difference, AVI − tool (95% interval; pooled: DerSimonian-Laird diamond, Hartung-Knapp line)",
                  fontsize=6.5, color=INK2, y=0.02)
    fig.tight_layout(w_pad=0.8, rect=(0, 0.04, 1, 1))
    save(fig, "figure3_part2_forest")


def figure4() -> None:
    bp = json.load(open(R2 / "by_position_pooled.json", encoding="utf-8"))
    classes = [("exonic", "Exonic"), ("intronic 3-50", "Intronic 3–50 bp"), ("intronic >50", "Intronic >50 bp")]
    tools = ["SpliceAI", "SpliceAI masked", "Pangolin", "CADD"]
    fig, ax = plt.subplots(figsize=(120 * MM, 58 * MM))
    offs = {t: 0.33 - i * 0.22 for i, t in enumerate(tools)}
    for r, (key, label) in enumerate(classes):
        yy = len(classes) - 1 - r
        for t in tools:
            p = bp.get(key, {}).get(f"AVI - {t}")
            if not p:
                continue
            hk = p.get("ci_hartung_knapp")
            if hk and np.all(np.isfinite(hk)):
                ax.plot(hk, [yy + offs[t]] * 2, color=COLOR[t], lw=0.5)
            ax.plot(p["ci"], [yy + offs[t]] * 2, color=COLOR[t], lw=1.6, solid_capstyle="round")
            hollow = t == "SpliceAI masked"
            ax.plot(p["estimate"], yy + offs[t], "o", ms=3.8, mfc="white" if hollow else COLOR[t], mec=COLOR[t],
                    mew=1.0)
    ax.axvline(0, color=BASE, lw=0.7, zorder=0)
    ax.set_yticks(range(len(classes)))
    ax.set_yticklabels([lab for _, lab in classes][::-1])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.set_xlabel("Pooled AUROC difference, AVI − tool (thick: DerSimonian-Laird 95%; thin: Hartung-Knapp)")
    handles = [plt.Line2D([], [], color=COLOR[t], marker="o", ms=3.8, lw=1.4,
                          mfc="white" if t == "SpliceAI masked" else COLOR[t], label=f"AVI − {t}") for t in tools]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, handlelength=1.4)
    fig.tight_layout()
    save(fig, "figure4_part2_by_position")


def figure_s1() -> None:
    p1 = pd.read_csv(R1 / "paired_differences.tsv", sep="\t")
    pairs = pd.read_csv(R2 / "paired_differences.tsv", sep="\t")
    pooled = json.load(open(R2 / "pooled.json", encoding="utf-8"))["pooled"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(180 * MM, 64 * MM))
    rows = []
    for key, label in REGIONS:
        r = p1[(p1.set == "main") & (p1.region == key) & (p1.pair == "AVI - SpliceAI2")]
        rows.append((f"{label} ({int(r['n_pos'].iloc[0])} P/LP)",
                     {"SpliceAI2": (r["difference"].iloc[0], *ci(r["ci"].iloc[0]))}))
    diff_rows(a, rows, ["SpliceAI2"], {"SpliceAI2": 0})
    a.set_xlabel("AUROC difference, AVI − SpliceAI2 (95% interval)")
    a.set_title("a  Part 1, clinical labels")
    datasets = [d for d in DATASETS if d[0] != "opensplice"]
    comparison = [("opensplice", "OpenSplice (developer-evaluated)")] + COMPARISON
    forest(b, pairs, pooled, "SpliceAI2", datasets, comparison, "b  Part 2, assays")
    b.set_xlabel("AUROC difference, AVI − SpliceAI2 (95% interval)")
    fig.tight_layout(w_pad=2.5)
    save(fig, "figureS1_spliceai2")


def figure_s2() -> None:
    pub = pd.read_csv(R1 / "published_thresholds.tsv", sep="\t")
    cadd = pub[pub["tool"].str.startswith("CADD")]
    fig, ax = plt.subplots(figsize=(120 * MM, 56 * MM))
    target = {"PP3 supporting": 2.406, "PP3 moderate": 5.790, "BP4 supporting": 1 / 2.406,
              "BP4 moderate": 1 / 5.790, "BP4 strong": 1 / 33.53}
    rows = list(cadd.itertuples())
    for i, r in enumerate(rows):
        yy = len(rows) - 1 - i
        lo, hi = ci(r.ci)
        ax.plot([target[r.level]] * 2, [yy - 0.35, yy + 0.35], color=INK2, lw=1.2)
        ax.plot([max(lo, 0.01), hi], [yy] * 2, color=COLOR["CADD"], lw=1.0)
        ax.plot(r.lr, yy, "o", ms=4, color=COLOR["CADD"])
        rng = f"({r.lo:g}, {r.hi:g}]" if r.side == "le" else f"[{r.lo:g}, {r.hi:g})"
        ax.text(1.02, yy, f"{r.fraction_pos:.0%} of P/LP", transform=ax.get_yaxis_transform(), va="center",
                fontsize=6, color=INK2)
        rows[i] = (f"{r.level}  {rng.replace('inf', '∞').replace('-∞', '−∞')}", r)
    ax.set_xscale("log")
    ax.set_xlim(0.01, 1000)
    ax.axvline(1, color=BASE, lw=0.7, zorder=0)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([lab for lab, _ in rows][::-1])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.set_xlabel("Likelihood ratio on intronic and UTR variants of the main set (log scale)")
    ax.set_title("Missense CADD intervals (Pejaver et al. 2022) on intronic variants")
    handles = [plt.Line2D([], [], color=COLOR["CADD"], marker="o", ms=4, lw=1.0, label="CADD, interval LR (95%)"),
               plt.Line2D([], [], color=INK2, lw=1.2, marker="|", ms=7, ls="none", label="ClinGen target")]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -0.3), ncol=2, handlelength=1.4)
    fig.tight_layout()
    save(fig, "figureS2_cadd_missense_thresholds")


FIGURES = {"1": figure1, "2": figure2, "3": figure3, "4": figure4, "S1": figure_s1, "S2": figure_s2}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--only", nargs="*", choices=list(FIGURES))
    a = ap.parse_args(argv)
    style()
    for key in a.only or FIGURES:
        FIGURES[key]()
    return 0


if __name__ == "__main__":
    sys.exit(main())
