"""Shared style, colours and drawing helpers of the manuscript figures. Every tool keeps one colour in every
figure; the palette was checked for colour-vision-deficiency separation (all checks pass; three colours are
below 3:1 contrast on white, so legends and labels carry identity as well)."""
from __future__ import annotations

import ast

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from paths import ROOT  # noqa: E402

R1, R2 = ROOT / "results" / "part1", ROOT / "results" / "part2"
OUT = ROOT / "manuscript" / "figures"
MAIN = "intronic >=3 bp (main)"
COLOR = {"AVI": "#2a78d6", "SpliceAI": "#eb6834", "SpliceAI masked": "#eb6834", "Pangolin": "#1baf7a",
         "CADD": "#eda100", "phyloP 100-way": "#e87ba4", "SpliceAI2": "#4a3aa7"}
INK, INK2, MUTED, GRID, BASE, BAND = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#f0efec"
TARGETS = {"supporting": 2.406, "moderate": 5.790, "strong": 33.53, "very_strong": 1124.0}
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
    """A stored interval as floats; 'inf' and 'nan' kept."""
    if isinstance(text, (list, tuple)):
        return tuple(float(x) for x in text)
    s = str(text).replace("nan", "'nan'").replace("inf", "'inf'")
    return tuple(float(x) for x in ast.literal_eval(s))


def save(fig, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("written", OUT / f"{name}.pdf")


def clean_rows_axis(ax, labels: list) -> None:
    """Row labels on the y axis (top to bottom), no y spine, vertical hairline grid."""
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels[::-1])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_axisbelow(True)


def diff_rows(ax, rows: list, tools: list, emphasise_last: bool = False) -> None:
    """Dot-and-interval rows: rows = [(label, {tool: (estimate, lo, hi)})], drawn top to bottom. With
    `emphasise_last` the last row is moved to the top, on a light band, with a bold label."""
    if emphasise_last:
        rows = [rows[-1]] + rows[:-1]
    offsets = {t: (len(tools) - 1) / 2 * 0.22 - i * 0.22 for i, t in enumerate(tools)}
    n = len(rows)
    for r, (_, vals) in enumerate(rows):
        yy = n - 1 - r
        for t in tools:
            if t not in vals:
                continue
            est, lo, hi = vals[t]
            hollow = t == "SpliceAI masked"
            ax.plot([lo, hi], [yy + offsets[t]] * 2, color=COLOR[t], lw=1.0, solid_capstyle="round", zorder=2)
            ax.plot(est, yy + offsets[t], "o", ms=3.6, mfc="white" if hollow else COLOR[t], mec=COLOR[t],
                    mew=1.0, zorder=3)
    ax.axvline(0, color=BASE, lw=0.7, zorder=1)
    clean_rows_axis(ax, [lab for lab, _ in rows])
    if emphasise_last:
        ax.axhspan(n - 1.5, n - 0.5, color=BAND, zorder=0, lw=0)
        ax.get_yticklabels()[-1].set_fontweight("bold")
        ax.get_yticklabels()[-1].set_color(INK)


def tool_handles(tools: list, prefix: str = "AVI − ") -> list:
    return [plt.Line2D([], [], color=COLOR[t], marker="o", ms=3.6, lw=1.0,
                       mfc="white" if t == "SpliceAI masked" else COLOR[t], label=f"{prefix}{t}") for t in tools]


def roc_points(s: np.ndarray, y: np.ndarray, w: np.ndarray) -> tuple:
    """Weighted ROC curve: false and true positive rates at every distinct threshold, descending."""
    order = np.argsort(-s, kind="mergesort")
    s, y, w = s[order], y[order], w[order]
    last = np.r_[np.nonzero(np.diff(s))[0], len(s) - 1]
    tp, fp = np.cumsum(w * (y == 1))[last], np.cumsum(w * (y == 0))[last]
    return np.r_[0, fp / fp[-1]], np.r_[0, tp / tp[-1]]


__all__ = ["R1", "R2", "OUT", "MAIN", "COLOR", "INK", "INK2", "MUTED", "GRID", "BASE", "BAND", "TARGETS", "MM",
           "style", "ci", "save", "clean_rows_axis", "diff_rows", "tool_handles", "roc_points", "plt", "np"]
