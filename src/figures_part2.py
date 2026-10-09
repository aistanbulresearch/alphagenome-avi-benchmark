"""Part 2 figures: forest plots (Figure 4), position relative to the exon (Figure 5) and the SpliceAI2
comparison over both parts (Figure S1)."""
from __future__ import annotations

import json

import pandas as pd

from figures_common import (BASE, COLOR, GRID, INK, INK2, MM, MUTED, R1, R2, ci, clean_rows_axis, diff_rows,
                            np, plt, save, tool_handles)

DATASETS = [("braun2018_mst1r_ex11", "MST1R exon 11"), ("compass", "COMPASS"), ("cortes_lopez2022_cd19", "CD19"),
            ("kitzman_pou1f1_ex2", "POU1F1 exon 2"), ("kitzman_wt1_ex9", "WT1 exon 9"), ("opensplice", "OpenSplice")]
AVI_BENCH = ("opensplice, developers' benchmark positions", "OpenSplice at AVI benchmark positions")
COMPARISON = [("mfass", "MFASS"), AVI_BENCH]
CLASSES = [("exonic", "Exonic"), ("intronic 3-50", "Intronic 3–50 bp"), ("intronic >50", "Intronic >50 bp")]


def load2() -> tuple:
    pairs = pd.read_csv(R2 / "paired_differences.tsv", sep="\t")
    pooled = json.load(open(R2 / "pooled.json", encoding="utf-8"))
    return pairs, pooled["pooled"], pooled["counts"]


def forest(ax, pairs, pooled: dict, tool: str, datasets: list, comparison: list, title: str) -> None:
    """Per-dataset differences (squares sized by random-effects weight), the pooled diamond (DerSimonian-Laird
    95%) with the Hartung-Knapp interval as a line, and comparison datasets in grey (not pooled)."""
    p = pooled.get(f"AVI - {tool}")
    rows = []
    for key, label in datasets:
        r = pairs[(pairs.dataset == key) & (pairs.pair == f"AVI - {tool}")]
        rows.append((label, None if r.empty else (r["difference"].iloc[0], *ci(r["ci"].iloc[0]), r["se"].iloc[0]),
                     "main"))
    rows.append(("Pooled", p, "pooled"))
    for key, label in comparison:
        r = pairs[(pairs.dataset == key) & (pairs.pair == f"AVI - {tool}")]
        rows.append((label, None if r.empty else (r["difference"].iloc[0], *ci(r["ci"].iloc[0]), r["se"].iloc[0]),
                     "comparison"))
    weights = {}
    if p:
        w = {lab: 1 / (v[3] ** 2 + p["tau2"]) for lab, v, k in rows if k == "main" and v is not None}
        weights = {lab: x / max(w.values()) for lab, x in w.items()}
    n = len(rows)
    col, hollow = COLOR[tool], tool == "SpliceAI masked"
    ys = []
    for i, (label, val, kind) in enumerate(rows):
        yy = n - 1 - i - (0.6 if kind == "comparison" else 0)
        ys.append(yy)
        if val is None:
            ax.text(0, yy, "not scored" if kind != "pooled" else "", ha="center", va="center", fontsize=5.5,
                    color=MUTED)
            continue
        if kind == "pooled":
            est, (lo, hi) = val["estimate"], val["ci"]
            hk = val.get("ci_hartung_knapp")
            if hk:
                ax.plot(hk, [yy] * 2, color=INK, lw=0.7)
            ax.fill([lo, est, hi, est], [yy, yy + 0.3, yy, yy - 0.3], color=INK, lw=0)
            ax.text(1.0, yy - 0.42, f"k = {val['k']}, I² = {val['i2']:.0%}", transform=ax.get_yaxis_transform(),
                    ha="right", va="top", fontsize=5.5, color=INK2)
            continue
        est, lo, hi = val[:3]
        c = MUTED if kind == "comparison" else col
        ax.plot([lo, hi], [yy] * 2, color=c, lw=1.0, solid_capstyle="round")
        ms = 2.2 + 3.6 * np.sqrt(weights.get(label, 0.4)) if kind == "main" else 3.2
        ax.plot(est, yy, "s" if kind == "main" else "o", ms=ms, mfc="white" if hollow and kind == "main" else c,
                mec=c, mew=1.0 if hollow else 0.5)
    ax.axvline(0, color=BASE, lw=0.7, zorder=0)
    ax.axhline(n - 1 - len(datasets) + 0.5, color=GRID, lw=0.6)
    ax.axhline(n - 1 - len(datasets) - 0.5 - 0.3, color=GRID, lw=0.6)
    ax.set_yticks(ys)
    ax.set_yticklabels([lab for lab, _, _ in rows])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.set_title(title)


def labelled(datasets: list, counts: dict) -> list:
    return [(k, f"{lab} ({counts[k]['variants']:,})") for k, lab in datasets]


def figure_forest() -> None:
    pairs, pooled, counts = load2()
    tools = [("SpliceAI", "a  SpliceAI (raw)"), ("SpliceAI masked", "b  SpliceAI (masked)"),
             ("Pangolin", "c  Pangolin"), ("CADD", "d  CADD")]
    fig, axes = plt.subplots(1, 4, figsize=(180 * MM, 72 * MM), sharey=True, sharex=True)
    for ax, (t, title) in zip(axes, tools):
        forest(ax, pairs, pooled, t, labelled(DATASETS, counts), labelled(COMPARISON, counts), title)
    for ax in axes[1:]:
        ax.tick_params(axis="y", labelleft=False)
    axes[0].set_xlim(-0.4, 0.25)
    fig.supxlabel("AUROC difference, AVI − tool (95% interval). Squares sized by random-effects weight; pooled: "
                  "DerSimonian-Laird diamond, Hartung-Knapp line", fontsize=6.5, color=INK2, y=0.02)
    fig.tight_layout(w_pad=0.8, rect=(0, 0.04, 1, 1))
    save(fig, "figure4_part2_forest")


def position_panel(ax, bp: dict, tools: list) -> None:
    offs = {t: (len(tools) - 1) / 2 * 0.22 - i * 0.22 for i, t in enumerate(tools)}
    labels = []
    for r, (key, label) in enumerate(CLASSES):
        yy = len(CLASSES) - 1 - r
        ks = sorted({v["k"] for t, v in bp.get(key, {}).items() if t.split(" - ")[1] in tools})
        labels.append(f"{label} (k = {'/'.join(map(str, ks))})")
        for t in tools:
            p = bp.get(key, {}).get(f"AVI - {t}")
            if not p:
                continue
            hk = p.get("ci_hartung_knapp")
            if hk and np.all(np.isfinite(hk)):
                ax.plot(hk, [yy + offs[t]] * 2, color=COLOR[t], lw=0.5)
            ax.plot(p["ci"], [yy + offs[t]] * 2, color=COLOR[t], lw=1.6, solid_capstyle="round")
            ax.plot(p["estimate"], yy + offs[t], "o", ms=3.8, mfc="white" if t == "SpliceAI masked" else COLOR[t],
                    mec=COLOR[t], mew=1.0)
    ax.axvline(0, color=BASE, lw=0.7, zorder=0)
    clean_rows_axis(ax, labels)


def figure_position() -> None:
    bp = json.load(open(R2 / "by_position_pooled.json", encoding="utf-8"))
    tools = ["SpliceAI", "SpliceAI masked", "Pangolin", "CADD"]
    fig, ax = plt.subplots(figsize=(125 * MM, 60 * MM))
    position_panel(ax, bp, tools)
    ax.set_xlabel("Pooled AUROC difference, AVI − tool (thick: DerSimonian-Laird 95%; thin: Hartung-Knapp)")
    ax.legend(handles=tool_handles(tools), loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2, handlelength=1.4)
    fig.tight_layout()
    save(fig, "figure5_part2_by_position")


def figure_spliceai2() -> None:
    from figures_part1 import region_rows
    p1 = pd.read_csv(R1 / "paired_differences.tsv", sep="\t")
    d1 = pd.read_csv(R1 / "discrimination.tsv", sep="\t")
    pairs, pooled, counts = load2()
    bp = json.load(open(R2 / "by_position_pooled.json", encoding="utf-8"))
    fig = plt.figure(figsize=(180 * MM, 120 * MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 0.8])
    a, b, c = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0])
    diff_rows(a, region_rows(p1, d1, ["SpliceAI2"]), ["SpliceAI2"], emphasise_last=True)
    a.set_xlabel("AUROC difference, AVI − SpliceAI2 (95% interval)")
    a.set_title("a  Part 1, clinical labels")
    datasets = [d for d in DATASETS if d[0] != "opensplice"]
    comparison = [("opensplice", "OpenSplice (evaluated by SpliceAI2 developers)")] + COMPARISON
    forest(b, pairs, pooled, "SpliceAI2", labelled(datasets, counts), labelled(comparison, counts),
           "b  Part 2, assays")
    b.set_xlabel("AUROC difference, AVI − SpliceAI2 (95% interval)")
    position_panel(c, bp, ["SpliceAI2"])
    c.set_xlabel("Pooled AUROC difference, AVI − SpliceAI2 (thick: DL 95%; thin: HK)")
    c.set_title("c  Part 2 by position relative to the exon (OpenSplice not pooled)")
    fig.tight_layout(w_pad=2.5, h_pad=2.0)
    save(fig, "figureS1_spliceai2")


__all__ = ["figure_forest", "figure_position", "figure_spliceai2"]
