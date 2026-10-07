"""ClinGen calibration of a score (Pejaver et al. 2022, Am J Hum Genet 109:2163).

Local positive likelihood ratio at score s: within the smallest symmetric interval [s - e, s + e] that
holds at least MIN_WINDOW labelled variants,

    lr+(s) = (P/LP in the window / all P/LP) / (B/LB in the window / all B/LB).

This equals the local posterior odds divided by the prior odds once the benign class is weighted to the
prior (Pejaver et al.), so the prior (0.0441) enters only through the evidence levels below. A
one-sided 95% bound comes from bootstrap resampling of the variants (10,000 replicates in the paper).

Pathogenic threshold for a level L: the smallest grid score t such that the lower bound of lr+ is at
least L at every grid score >= t. Benign: the largest t such that the upper bound is at most 1/L at
every grid score <= t. Higher scores are taken as more pathogenic.

Simplification stated in the Methods: the paper's windows also required at least 3% of a set of rare
gnomAD variants; this study has no such reference set, so only the labelled-variant count is used.
"""
from __future__ import annotations

import numpy as np

PRIOR = 0.0441
LEVELS = {"supporting": 2.406, "moderate": 5.790, "strong": 33.53, "very_strong": 1124.0}
MIN_WINDOW = 100
GRID = 400
BOOT = 10_000
SEED = 20261007


def grid_points(scores: np.ndarray, n: int = GRID) -> np.ndarray:
    return np.unique(np.quantile(scores, np.linspace(0.001, 0.999, n)))


def windows(scores: np.ndarray, grid: np.ndarray, k: int = MIN_WINDOW) -> np.ndarray:
    """Half-width e at each grid point: the k-th smallest distance to a labelled score."""
    k = min(k, len(scores))
    return np.asarray([np.partition(np.abs(scores - g), k - 1)[k - 1] for g in grid])


def local_lr(scores: np.ndarray, labels: np.ndarray, grid: np.ndarray, half: np.ndarray) -> np.ndarray:
    pos, neg = np.sort(scores[labels == 1]), np.sort(scores[labels == 0])
    lo, hi = grid - half, grid + half

    def count(x):
        return np.searchsorted(x, hi, side="right") - np.searchsorted(x, lo, side="left")
    p, b = count(pos) / max(len(pos), 1), count(neg) / max(len(neg), 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(b > 0, p / b, np.inf)


def bootstrap_bounds(scores: np.ndarray, labels: np.ndarray, grid: np.ndarray, half: np.ndarray,
                     replicates: int = BOOT, seed: int = SEED) -> tuple[np.ndarray, np.ndarray]:
    """One-sided 95% bounds of lr+ per grid point: (5th percentile, 95th percentile)."""
    rng = np.random.default_rng(seed)
    n = len(scores)
    reps = np.empty((replicates, len(grid)))
    for r in range(replicates):
        idx = rng.integers(0, n, n)
        reps[r] = local_lr(scores[idx], labels[idx], grid, half)
    # no interpolation: an infinite lr+ (no B/LB in the window) must not turn a bound into NaN
    return np.quantile(reps, 0.05, axis=0, method="lower"), np.quantile(reps, 0.95, axis=0, method="higher")


def thresholds(grid: np.ndarray, lower: np.ndarray, upper: np.ndarray, levels: dict | None = None) -> dict:
    """{level: (pathogenic threshold or None, benign threshold or None)} for the given LR targets."""
    out = {}
    for name, lr in (levels or LEVELS).items():
        ok_p = lower >= lr
        path = None
        for i in range(len(grid)):
            if ok_p[i:].all():
                path = float(grid[i])
                break
        ok_b = upper <= 1 / lr
        benign = None
        for i in range(len(grid) - 1, -1, -1):
            if ok_b[: i + 1].all():
                benign = float(grid[i])
                break
        out[name] = (path, benign)
    return out


def calibrate(scores: np.ndarray, labels: np.ndarray, replicates: int = BOOT, seed: int = SEED) -> dict:
    keep = ~np.isnan(scores)
    s, y = scores[keep], labels[keep]
    grid = grid_points(s)
    half = windows(s, grid)
    lower, upper = bootstrap_bounds(s, y, grid, half, replicates, seed)
    return {"n_pos": int((y == 1).sum()), "n_neg": int((y == 0).sum()),
            "thresholds": thresholds(grid, lower, upper),
            "curve": {"grid": grid.tolist(), "lr": local_lr(s, y, grid, half).tolist(),
                      "lower": lower.tolist(), "upper": upper.tolist()}}


def likelihood_ratio(scores: np.ndarray, labels: np.ndarray, t: float, side: str,
                     weights: np.ndarray | None = None) -> float:
    """LR of the evidence 'score >= t' (side 'ge') or 'score <= t' (side 'le') on a test set."""
    w = np.ones(len(scores)) if weights is None else weights.astype(float)
    hit = scores >= t if side == "ge" else scores <= t
    p = (w * hit * (labels == 1)).sum() / max((w * (labels == 1)).sum(), 1e-12)
    b = (w * hit * (labels == 0)).sum() / max((w * (labels == 0)).sum(), 1e-12)
    return float(p / b) if b > 0 else float("inf")


def levels_for_prior(prior: float) -> dict:
    """ClinGen likelihood-ratio targets for another prior, by the rule behind Pejaver et al. (2022):
    LR = c^(1/8), c^(1/4), c^(1/2), c for supporting to very strong, with c the smallest value at which
    the weakest 'likely pathogenic' combination (exponent 6/8) reaches a posterior of 0.90, i.e.
    c = (9 (1 - prior) / prior)^(4/3). This gives c = 350 at prior 0.1 (Tavtigian et al. 2018) and
    1,131 at 0.0441 (the paper reports 1,124)."""
    big_c = (9 * (1 - prior) / prior) ** (4 / 3)
    return {"supporting": big_c ** 0.125, "moderate": big_c ** 0.25, "strong": big_c ** 0.5, "very_strong": big_c}


def thresholds_from_curve(curve: dict, levels: dict) -> dict:
    """Thresholds for other level targets from a stored curve (grid and one-sided bootstrap bounds)."""
    grid, lower, upper = (np.asarray(curve[k], dtype=float) for k in ("grid", "lower", "upper"))
    return thresholds(grid, lower, upper, levels)


def interval_lr(scores: np.ndarray, labels: np.ndarray, lo: float, hi: float, side: str,
                weights: np.ndarray | None = None) -> float:
    """LR of a score interval on a test set: pathogenic side lo <= s < hi, benign side lo < s <= hi."""
    w = np.ones(len(scores)) if weights is None else weights.astype(float)
    hit = (scores >= lo) & (scores < hi) if side == "ge" else (scores > lo) & (scores <= hi)
    p = (w * hit * (labels == 1)).sum() / max((w * (labels == 1)).sum(), 1e-12)
    b = (w * hit * (labels == 0)).sum() / max((w * (labels == 0)).sum(), 1e-12)
    return float(p / b) if b > 0 else (float("inf") if p > 0 else float("nan"))


def level_intervals(thresholds: dict) -> list:
    """[(level, side, lo, hi)]: each evidence level is the score interval between its threshold and the
    next stronger defined threshold (ClinGen: a variant gets the level of the interval it falls in)."""
    order = list(LEVELS)
    out = []
    path = [(lvl, thresholds[lvl][0]) for lvl in order if thresholds[lvl][0] is not None]
    for i, (lvl, lo) in enumerate(path):
        hi = path[i + 1][1] if i + 1 < len(path) else float("inf")
        out.append((lvl, "ge", lo, hi))
    ben = [(lvl, thresholds[lvl][1]) for lvl in order if thresholds[lvl][1] is not None]
    for i, (lvl, hi) in enumerate(ben):
        lo = ben[i + 1][1] if i + 1 < len(ben) else float("-inf")
        out.append((lvl, "le", lo, hi))
    return out


__all__ = ["PRIOR", "LEVELS", "calibrate", "likelihood_ratio", "interval_lr", "level_intervals", "levels_for_prior", "thresholds_from_curve", "thresholds", "local_lr", "windows", "grid_points"]
