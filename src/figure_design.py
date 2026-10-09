"""Figure 1: study design. Counts come from the variant list summary, the calibration sample and the Part 2
pooled-table counts, so the diagram follows the data."""
from __future__ import annotations

import json

import pandas as pd
from matplotlib.patches import FancyBboxPatch

from figures_common import BAND, BASE, INK, INK2, MM, R1, R2, plt, save
from paths import DERIVED

NAMES = {"braun2018_mst1r_ex11": "MST1R exon 11 (minigene)", "compass": "COMPASS (minigene libraries)",
         "cortes_lopez2022_cd19": "CD19 (minigene)", "kitzman_pou1f1_ex2": "POU1F1 exon 2 (minigene)",
         "kitzman_wt1_ex9": "WT1 exon 9 (minigene)", "opensplice": "OpenSplice (MPRA, 599 exons)",
         "mfass": "MFASS",
         "opensplice, developers' benchmark positions": "OpenSplice at AVI benchmark positions"}


def box(ax, x, y, w, h, text, fill="white", bold_first=True, size=6.3):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.012", fc=fill, ec=BASE,
                                lw=0.7, transform=ax.transAxes, clip_on=False))
    lines = text.split("\n")
    ax.text(x + 0.012, y + h - 0.018, lines[0], transform=ax.transAxes, va="top", ha="left", fontsize=size,
            color=INK, fontweight="bold" if bold_first else "normal")
    if len(lines) > 1:
        ax.text(x + 0.012, y + h - 0.052, "\n".join(lines[1:]), transform=ax.transAxes, va="top", ha="left",
                fontsize=size - 0.5, color=INK2, linespacing=1.35)


def arrow(ax, x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), xycoords="axes fraction",
                arrowprops={"arrowstyle": "-|>", "color": INK2, "lw": 0.7, "mutation_scale": 7})


def figure_design() -> None:
    v = pd.read_csv(DERIVED / "variants" / "clinvar_variants.tsv.gz", sep="\t", dtype=str, keep_default_na=False,
                    usecols=["new_label", "exclusion"])
    excl = v["exclusion"].value_counts().to_dict()
    n_all, n_new = len(v), int((v["new_label"] == "yes").sum())
    calib = pd.read_csv(DERIVED / "variants" / "calibration_sample.tsv", sep="\t", dtype=str)["label"].value_counts()
    cov = pd.read_csv(R1 / "coverage.tsv", sep="\t")
    n_main = int(cov[(cov.set == "main") & (cov.tool == "AVI")]["variants"].iloc[0])
    disc = pd.read_csv(R1 / "discrimination.tsv", sep="\t")
    main_row = disc[(disc.set == "main") & (disc.region == "intronic >=3 bp (main)") & (disc.tool == "AVI")].iloc[0]
    counts = json.load(open(R2 / "pooled.json", encoding="utf-8"))["counts"]

    fig, ax = plt.subplots(figsize=(180 * MM, 125 * MM))
    ax.axis("off")
    ax.text(0.0, 0.985, "Part 1  Clinical labels added after the AlphaGenome training snapshot (2025-06-15)",
            transform=ax.transAxes, fontsize=7.2, fontweight="bold", color=INK, va="top")
    box(ax, 0.0, 0.80, 0.40, 0.13, f"ClinVar, October 2026 release\n{n_all:,} variants classified P/LP or B/LB "
                                   "(GRCh38)")
    box(ax, 0.44, 0.69, 0.53, 0.24,
        "Excluded\n"
        f"{excl.get('labelled before', 0):,} already P/LP or B/LB before the snapshot\n"
        "        (ClinVar 2025-06-15 VCF; 2025-07 and 2024-07 releases)\n"
        f"{excl.get(chr(100) + 'evelopers' + chr(39) + ' benchmark or submitter', 0):,} at AVI developers' "
        "benchmark positions or from a development submitter\n"
        f"{excl.get('unresolved historical identity (same position)', 0):,} older record of unresolved identity "
        "at the same position\n"
        f"{sum(n for k, n in excl.items() if k.startswith('identity collision')):,} identity collisions",
        fill=BAND)
    arrow(ax, 0.20, 0.80, 0.20, 0.74)
    arrow(ax, 0.20, 0.775, 0.44, 0.775)
    box(ax, 0.0, 0.62, 0.40, 0.12, f"New-label set\n{n_new:,} variants")
    arrow(ax, 0.20, 0.62, 0.20, 0.56)
    box(ax, 0.0, 0.40, 0.40, 0.16, f"Main set\n{n_main:,} variants: SNV, ≥1 review star, MANE Select name,\n"
                                   "untouched and validation chromosomes,\nrare in gnomAD v4.1 (FAF95 < 0.001)")
    arrow(ax, 0.20, 0.40, 0.20, 0.34)
    box(ax, 0.0, 0.22, 0.40, 0.12, f"Main comparison\nintronic ≥3 bp from an exon: {int(main_row['n_pos']):,} P/LP, "
                                   f"{int(main_row['n_neg']):,} B/LB", fill=BAND)
    box(ax, 0.44, 0.40, 0.53, 0.20,
        "Calibration sample (thresholds derived here)\n"
        f"labelled variants outside the new-label set: all {calib.get('P/LP', 0):,} P/LP and {calib.get('B/LB', 0):,} "
        "B/LB\n(random, stratified by region and chromosome group), non-canonical\nintronic and UTR, ≥1 star, rare. "
        "Thresholds tested on the main set.")
    arrow(ax, 0.44, 0.47, 0.40, 0.47)

    ax.text(0.0, 0.17, "Part 2  Experimental splicing measurements (no predictor selected the variants)",
            transform=ax.transAxes, fontsize=7.2, fontweight="bold", color=INK, va="top")
    main = ["braun2018_mst1r_ex11", "compass", "cortes_lopez2022_cd19", "kitzman_pou1f1_ex2", "kitzman_wt1_ex9",
            "opensplice"]

    def line(k):
        return f"{NAMES[k]}: {counts[k]['variants']:,} ({counts[k]['splice_altering']:,} splice-altering)"
    box(ax, 0.0, -0.13, 0.97, 0.26, "Six main datasets (pooled by random effects); comparison sets reported "
                                   "separately")
    nl = chr(10)
    ax.text(0.012, 0.075, nl.join(line(k) for k in main[:3]), transform=ax.transAxes, va="top", fontsize=5.8,
            color=INK2, linespacing=1.4)
    ax.text(0.50, 0.075, nl.join(line(k) for k in main[3:]), transform=ax.transAxes, va="top", fontsize=5.8,
            color=INK2, linespacing=1.4)
    other = [k for k in ("mfass", "opensplice, developers' benchmark positions") if k in counts]
    ax.text(0.012, -0.04, "Comparison: " + "; ".join(line(k) for k in other) + "." + nl + "Alleles of the Part 1 "
            "new-label set and at AVI developers' benchmark positions were removed from the main datasets.",
            transform=ax.transAxes, va="top", fontsize=5.8, color=INK2, linespacing=1.4)
    ax.text(0.0, -0.155, "Tools: AVI against SpliceAI (raw and masked), Pangolin, CADD v1.7 and conservation "
            "(phyloP 100-way; AVI's inputs phastCons 470-way and phyloP 241-way)." + nl + "Added after the main "
            "results: SpliceAI2.", transform=ax.transAxes, va="top", fontsize=6, color=INK2, linespacing=1.4)
    save(fig, "figure1_study_design")


__all__ = ["figure_design"]
