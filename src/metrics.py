"""Discrimination metrics with integer weights, and the gene-cluster bootstrap.

A bootstrap replicate resamples genes with replacement; a variant then carries the number of times its
gene was drawn as its weight. AUROC and AUPRC are computed with those weights, so one sort per tool and
set serves every replicate.
"""
from __future__ import annotations

import numpy as np

SEED = 20261007
REPLICATES = 2000


def _order(scores: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Indices that sort scores ascending, and the group id of each tied block in that order."""
    order = np.argsort(scores, kind="mergesort")
    s = scores[order]
    group = np.concatenate([[0], np.cumsum(s[1:] != s[:-1])])
    return order, group


def auroc(scores: np.ndarray, labels: np.ndarray, weights: np.ndarray | None = None,
          prepared: tuple | None = None) -> float:
    """Weighted AUROC (ties count one half); labels 1 = positive. NaN if a class has no weight."""
    w = np.ones(len(scores)) if weights is None else weights.astype(float)
    order, group = prepared or _order(scores)
    lab, ww = labels[order], w[order]
    pos_w = np.bincount(group, weights=ww * (lab == 1))
    neg_w = np.bincount(group, weights=ww * (lab == 0))
    n_pos, n_neg = pos_w.sum(), neg_w.sum()
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    neg_below = np.cumsum(neg_w) - neg_w
    return float((pos_w * (neg_below + 0.5 * neg_w)).sum() / (n_pos * n_neg))


def auprc(scores: np.ndarray, labels: np.ndarray, weights: np.ndarray | None = None,
          prepared: tuple | None = None) -> float:
    """Weighted average precision over the distinct score thresholds (descending)."""
    w = np.ones(len(scores)) if weights is None else weights.astype(float)
    order, group = prepared or _order(scores)
    lab, ww = labels[order], w[order]
    pos_w = np.bincount(group, weights=ww * (lab == 1))[::-1]
    all_w = np.bincount(group, weights=ww)[::-1]
    tp, cum = np.cumsum(pos_w), np.cumsum(all_w)
    if tp[-1] == 0:
        return float("nan")
    precision = np.divide(tp, cum, out=np.zeros_like(tp), where=cum > 0)
    return float((precision * pos_w).sum() / tp[-1])


def cluster_weights(clusters: np.ndarray, replicates: int = REPLICATES, seed: int = SEED):
    """Yield, per replicate, the weight of each variant: how often its cluster was drawn."""
    uniq, inverse = np.unique(clusters, return_inverse=True)
    rng = np.random.default_rng(seed)
    for _ in range(replicates):
        drawn = np.bincount(rng.integers(0, len(uniq), len(uniq)), minlength=len(uniq))
        yield drawn[inverse]


def clusters_of(genes: list) -> np.ndarray:
    """Gene ids as clusters; a variant without a single gene id is its own cluster."""
    out = []
    for i, g in enumerate(genes):
        g = (g or "").strip()
        out.append(g if g and g not in ("-", "-1") and ";" not in g else f"variant:{i}")
    return np.asarray(out)


def interval(values: list, level: float = 0.95) -> tuple[float, float]:
    v = np.asarray([x for x in values if not np.isnan(x)])
    if len(v) == 0:
        return float("nan"), float("nan")
    a = (1 - level) / 2
    return float(np.quantile(v, a)), float(np.quantile(v, 1 - a))


__all__ = ["SEED", "REPLICATES", "auroc", "auprc", "cluster_weights", "clusters_of", "interval"]
