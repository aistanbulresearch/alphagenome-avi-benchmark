"""Part 1 figures: discrimination (Figure 2), calibration (Figure 3) and the missense CADD thresholds on
intronic variants (Figure S2)."""
from __future__ import annotations

import json

import pandas as pd

from figures_common import (BASE, COLOR, GRID, INK, INK2, MAIN, MM, R1, TARGETS, ci, clean_rows_axis,
                            diff_rows, np, plt, roc_points, save, tool_handles)

REGIONS = [("donor region", "Donor region"), ("acceptor region", "Acceptor region"),
           ("near intronic", "Near intronic"), ("deep intronic", "Deep intronic"),
           (MAIN, "All intronic ≥3 bp (main comparison)")]
LEVELS = [("supporting pathogenic", "Supporting", "supporting", False),
          ("moderate pathogenic", "Moderate", "moderate", False),
          ("strong pathogenic", "Strong", "strong", False),
          ("very_strong pathogenic", "Very strong", "very_strong", False),
          ("supporting benign", "Supporting (benign)", "supporting", True),
          ("moderate benign", "Moderate (benign)", "moderate", True)]
LO_LIM, HI_LIM = 1e-3, 3000


def region_rows(pairs: pd.DataFrame, disc: pd.DataFrame, tools: list) -> list:
    rows = []
    for key, label in REGIONS:
        vals = {}
        for t in tools:
            r = pairs[(pairs.set == "main") & (pairs.region == key) & (pairs.pair == f"AVI - {t}")]
            if len(r):
                vals[t] = (r["difference"].iloc[0], *ci(r["ci"].iloc[0]))
        npos = int(disc[(disc.set == "main") & (disc.region == key) & (disc.tool == "AVI")]["n_pos"].iloc[0])
        rows.append((f"{label} ({npos} P/LP)", vals))
    return rows


def figure_discrimination() -> None:
    import part1_analysis as pa
    import part1_data as pdata
    v = pdata.load()
    mask = pdata.sets(v)["main"] & pdata.regions(v)[MAIN]
    disc = pd.read_csv(R1 / "discrimination.tsv", sep="\t")
    pairs = pd.read_csv(R1 / "paired_differences.tsv", sep="\t")
    fig, (a, b) = plt.subplots(1, 2, figsize=(180 * MM, 72 * MM), gridspec_kw={"width_ratios": [1, 1.25]})
    for t in ["AVI", "SpliceAI", "Pangolin", "CADD", "phyloP 100-way"]:
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
    tools = ["SpliceAI", "Pangolin", "CADD"]
    diff_rows(b, region_rows(pairs, disc, tools), tools, emphasise_last=True)
    b.set_xlabel("AUROC difference, AVI − tool (95% gene-bootstrap interval)")
    b.set_title("b  AVI − tool, by region")
    b.legend(handles=tool_handles(tools), loc="upper left", bbox_to_anchor=(0, -0.17), ncol=3, columnspacing=1.2,
             handlelength=1.4)
    fig.tight_layout(w_pad=2.5)
    save(fig, "figure2_part1_discrimination")


def calibration_panel(ax, cal: dict, tools: list) -> None:
    """Interval LRs with 95% intervals against each level's target. Open-ended interval: arrow. Interval with
    P/LP but no B/LB: no point estimate; a triangle marks the exact lower bound."""
    n = len(LEVELS)
    offs = {t: (len(tools) - 1) / 2 * 0.22 - i * 0.22 for i, t in enumerate(tools)}
    for r, (key, _label, lvl, benign) in enumerate(LEVELS):
        yy = n - 1 - r
        target = 1 / TARGETS[lvl] if benign else TARGETS[lvl]
        ax.plot([target, target], [yy - 0.38, yy + 0.38], color=INK2, lw=1.2, solid_capstyle="butt", zorder=1)
        for t in tools:
            e = cal[t]["evaluated_on_main"]["interval"].get(key)
            if not e or e["lr"] is None or np.isnan(e["lr"]):
                continue
            y0 = yy + offs[t]
            lo, hi = e["ci"]
            lo_d = max(lo, LO_LIM) if lo > 0 else LO_LIM
            hi_d = min(hi, HI_LIM) if np.isfinite(hi) else HI_LIM
            ax.plot([lo_d, hi_d], [y0] * 2, color=COLOR[t], lw=1.0, solid_capstyle="round", zorder=2)
            if not np.isfinite(hi) or hi > HI_LIM:
                ax.annotate("", xy=(HI_LIM * 1.3, y0), xytext=(HI_LIM, y0), annotation_clip=False,
                            arrowprops={"arrowstyle": "-|>", "color": COLOR[t], "lw": 0.8, "mutation_scale": 6})
            if lo < LO_LIM:
                ax.annotate("", xy=(LO_LIM / 1.3, y0), xytext=(LO_LIM, y0), annotation_clip=False,
                            arrowprops={"arrowstyle": "-|>", "color": COLOR[t], "lw": 0.8, "mutation_scale": 6})
            if np.isinf(e["lr"]):
                ax.plot(lo_d, y0, ">", ms=4.5, color=COLOR[t], zorder=3)
                ax.text(lo_d / 1.25, y0, f"{e['n_pos']} P/LP, no B/LB", ha="right", va="center", fontsize=5.5,
                        color=INK2, zorder=4, bbox={"boxstyle": "square,pad=0.1", "fc": "white", "ec": "none"})
            else:
                ax.plot(e["lr"], y0, "o", ms=3.6, color=COLOR[t], mec="white", mew=0.4, zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(LO_LIM, HI_LIM)
    ax.axvline(1, color=BASE, lw=0.7, zorder=0)
    ax.axhline(1.5, color=GRID, lw=0.6, zorder=0)
    clean_rows_axis(ax, [lab for _, lab, _, _ in LEVELS])
    ax.set_xlabel("Likelihood ratio of the score interval on the main set (log scale)")


def evidence_share(cal: dict, tool: str) -> tuple:
    iv = cal[tool]["evaluated_on_main"]["interval"]
    pos = sum(e["fraction_pos"] for k, e in iv.items() if k.endswith("pathogenic"))
    neg = sum(e["fraction_neg"] for k, e in iv.items() if k.endswith("benign"))
    return pos, neg


def figure_calibration() -> None:
    cal = json.load(open(R1 / "calibration.json", encoding="utf-8"))
    tools = ["AVI", "SpliceAI", "CADD"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(180 * MM, 66 * MM), gridspec_kw={"width_ratios": [1.6, 1]})
    calibration_panel(a, cal, tools)
    handles = [plt.Line2D([], [], color=COLOR[t], marker="o", ms=3.6, lw=1.0, label=t) for t in tools]
    handles.append(plt.Line2D([], [], color=INK2, lw=1.2, marker="|", ms=7, ls="none", label="ClinGen target"))
    if any(np.isinf(cal[t]["evaluated_on_main"]["interval"].get(k, {}).get("lr", 0) or 0)
           for t in tools for k, *_ in LEVELS):
        handles.append(plt.Line2D([], [], color=INK2, marker=">", ms=4.5, ls="none",
                                  label="exact lower bound, no B/LB"))
    a.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -0.25), ncol=3, handlelength=1.4)
    a.set_title("a  Evidence levels, each tested within its own score interval")
    for i, t in enumerate(tools):
        pos, neg = evidence_share(cal, t)
        strong = sum(e["fraction_pos"] for k, e in cal[t]["evaluated_on_main"]["interval"].items()
                     if k.startswith(("strong", "very_strong")) and k.endswith("pathogenic"))
        yy = len(tools) - 1 - i
        for dy, share, marker, fill in ((0.25, pos, "o", True), (0.0, strong, "D", True), (-0.25, neg, "o", False)):
            b.plot([0, share * 100], [yy + dy] * 2, color=COLOR[t], lw=1.0)
            b.plot(share * 100, yy + dy, marker, ms=4 if marker == "o" else 3.4, mfc=COLOR[t] if fill else "white",
                   mec=COLOR[t], mew=1.0)
            text = f"{share:.0%}" if share > 0 or dy >= 0 else "no benign level reached"
            b.text(share * 100 + 3, yy + dy, text, va="center", fontsize=6, color=INK2)
    clean_rows_axis(b, tools)
    b.set_xlim(0, 112)
    b.set_xticks([0, 25, 50, 75, 100])
    b.set_xlabel("Share of variants receiving evidence (%)")
    hb = [plt.Line2D([], [], color=INK2, marker="o", ms=4, ls="none", label="P/LP with any pathogenic evidence"),
          plt.Line2D([], [], color=INK2, marker="D", ms=3.4, ls="none", label="P/LP with strong evidence"),
          plt.Line2D([], [], color=INK2, marker="o", ms=4, mfc="white", ls="none", label="B/LB with benign evidence")]
    b.legend(handles=hb, loc="upper left", bbox_to_anchor=(0, -0.25), ncol=1, handlelength=1.0)
    b.set_title("b  Evidence coverage (prior 0.0441)")
    fig.tight_layout(w_pad=2.5)
    save(fig, "figure3_part1_calibration")


def figure_cadd_missense() -> None:
    pub = pd.read_csv(R1 / "published_thresholds.tsv", sep="\t")
    cadd = pub[pub["tool"].str.startswith("CADD")]
    fig, ax = plt.subplots(figsize=(140 * MM, 58 * MM))
    target = {"PP3 supporting": 2.406, "PP3 moderate": 5.790, "BP4 supporting": 1 / 2.406,
              "BP4 moderate": 1 / 5.790, "BP4 strong": 1 / 33.53}
    labels = []
    rows = list(cadd.itertuples())
    for i, r in enumerate(rows):
        yy = len(rows) - 1 - i
        lo, hi = ci(r.ci)
        ax.plot([target[r.level]] * 2, [yy - 0.35, yy + 0.35], color=INK2, lw=1.2)
        ax.plot([max(lo, 0.01), hi], [yy] * 2, color=COLOR["CADD"], lw=1.0)
        ax.plot(r.lr, yy, "o", ms=4, color=COLOR["CADD"])
        tr = ax.get_yaxis_transform()
        ax.text(1.03, yy, f"{r.fraction_pos:.0%}", transform=tr, va="center", ha="left", fontsize=6, color=INK2)
        ax.text(1.17, yy, f"{r.fraction_neg:.1%}", transform=tr, va="center", ha="left", fontsize=6, color=INK2)
        rng = f"({r.lo:g}, {r.hi:g}]" if r.side == "le" else f"[{r.lo:g}, {r.hi:g})"
        labels.append(f"{r.level}  {rng.replace('-inf', '−∞').replace('inf', '∞')}")
    tr = ax.get_yaxis_transform()
    ax.text(1.03, len(rows) - 0.4, "P/LP", transform=tr, fontsize=6, color=INK, fontweight="bold")
    ax.text(1.17, len(rows) - 0.4, "B/LB", transform=tr, fontsize=6, color=INK, fontweight="bold")
    ax.set_xscale("log")
    ax.set_xlim(0.01, 1000)
    ax.axvline(1, color=BASE, lw=0.7, zorder=0)
    clean_rows_axis(ax, labels)
    ax.set_xlabel("Likelihood ratio on intronic and UTR variants of the main set (log scale)")
    ax.set_title("Missense CADD intervals (Pejaver et al. 2022) on intronic variants", pad=14)
    handles = [plt.Line2D([], [], color=COLOR["CADD"], marker="o", ms=4, lw=1.0, label="CADD, interval LR (95%)"),
               plt.Line2D([], [], color=INK2, lw=1.2, marker="|", ms=7, ls="none", label="ClinGen target")]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -0.3), ncol=2, handlelength=1.4)
    fig.tight_layout()
    save(fig, "figureS2_cadd_missense_thresholds")


__all__ = ["figure_discrimination", "figure_calibration", "figure_cadd_missense", "region_rows", "REGIONS"]
