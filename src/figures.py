"""Manuscript figures (PDF and PNG in manuscript/figures/), from results/part1, results/part2 and, for the
ROC curves, the Part 1 score table.

    Figure 1  study design
    Figure 2  Part 1 discrimination: ROC curves (main region) and AVI − tool AUROC differences by region
    Figure 3  Part 1 calibration: interval likelihood ratios against the ClinGen targets; evidence coverage
    Figure 4  Part 2 forest plots, pooled by random effects (DerSimonian-Laird and Hartung-Knapp)
    Figure 5  Part 2 by position relative to the exon (additional analysis)
    Figure S1 SpliceAI2 (added comparator): Part 1 by region, Part 2 per dataset, Part 2 by position
    Figure S2 published missense CADD thresholds evaluated on intronic variants (additional analysis)

Usage:
    python src/figures.py [--only 1 2 ...]
"""
from __future__ import annotations

import argparse
import sys

from figure_design import figure_design
from figures_common import style
from figures_part1 import figure_cadd_missense, figure_calibration, figure_discrimination
from figures_part2 import figure_forest, figure_position, figure_spliceai2

FIGURES = {"1": figure_design, "2": figure_discrimination, "3": figure_calibration, "4": figure_forest,
           "5": figure_position, "S1": figure_spliceai2, "S2": figure_cadd_missense}


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
