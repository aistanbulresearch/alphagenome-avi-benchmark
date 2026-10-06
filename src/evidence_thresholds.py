"""Evidence-strength thresholds for computational evidence (Bayesian ACMG/AMP framework).

Primary thresholds: Pejaver et al., AJHG 2022, Table 1 (prior probability of pathogenicity 0.0441):
positive likelihood ratios 2.406 (supporting), 5.790 (moderate), 33.53 (strong), 1124 (very strong).
Benign thresholds are the reciprocals.

Alternative threshold set for the sensitivity analysis: Tavtigian et al., Genet Med 2018 (prior 0.10),
very-strong odds of pathogenicity 350, with the supporting, moderate and strong levels obtained by the
exponents 1/8, 1/4 and 1/2. The study reports results under both sets; it does not derive a new set for any
other prior. Posterior probabilities for any prior are computed with Bayes' rule.
"""
from __future__ import annotations

from dataclasses import dataclass

LEVELS = ("supporting", "moderate", "strong", "very_strong")
EXPONENTS = {"supporting": 1 / 8, "moderate": 1 / 4, "strong": 1 / 2, "very_strong": 1.0}

PEJAVER_2022_PRIOR = 0.0441
PEJAVER_2022_LR = {"supporting": 2.406, "moderate": 5.790, "strong": 33.53, "very_strong": 1124.0}
# Posterior values printed in Pejaver 2022 Table 1 (four significant digits) that the code must reproduce.
PEJAVER_2022_POSTERIOR_PRINTED = {"supporting": 0.0999, "very_strong": 0.9811}

TAVTIGIAN_2018_PRIOR = 0.10
TAVTIGIAN_2018_ODDS_VERY_STRONG = 350.0


@dataclass(frozen=True)
class ThresholdSet:
    """Likelihood-ratio thresholds of one framework."""

    name: str
    prior: float
    lr_pathogenic: dict

    def lr_benign(self) -> dict:
        """Benign (LR-) thresholds: reciprocals of the pathogenic thresholds."""
        return {level: 1.0 / value for level, value in self.lr_pathogenic.items()}


def posterior_from_lr(lr: float, prior: float) -> float:
    """Posterior probability of pathogenicity given a likelihood ratio and a prior (Bayes' rule)."""
    if lr <= 0 or not 0 < prior < 1:
        raise ValueError("likelihood ratio must be positive and prior in (0, 1)")
    return lr * prior / ((lr - 1.0) * prior + 1.0)


def scaled_thresholds(odds_very_strong: float) -> dict:
    """Thresholds of the four levels from the very-strong odds (exponents 1/8, 1/4, 1/2, 1)."""
    if odds_very_strong <= 1:
        raise ValueError("very-strong odds must exceed 1")
    return {level: odds_very_strong ** EXPONENTS[level] for level in LEVELS}


PRIMARY = ThresholdSet("Pejaver 2022, Table 1", PEJAVER_2022_PRIOR, dict(PEJAVER_2022_LR))
TAVTIGIAN_2018 = ThresholdSet("Tavtigian 2018", TAVTIGIAN_2018_PRIOR,
                              scaled_thresholds(TAVTIGIAN_2018_ODDS_VERY_STRONG))
SENSITIVITY_PRIORS = (0.10, 0.01)


def posterior_table(threshold_set: ThresholdSet, priors=None) -> dict:
    """Posterior probability at each pathogenic LR threshold for the set's own prior and extra priors."""
    priors = (threshold_set.prior,) + tuple(priors or ())
    return {prior: {level: posterior_from_lr(lr, prior) for level, lr in threshold_set.lr_pathogenic.items()}
            for prior in priors}


def strength_of(lr_lower_bound: float, threshold_set: ThresholdSet = PRIMARY) -> str | None:
    """Highest pathogenic level whose threshold the lower confidence bound of LR+ reaches, or None."""
    reached = None
    for level in LEVELS:
        if lr_lower_bound >= threshold_set.lr_pathogenic[level]:
            reached = level
    return reached


__all__ = [
    "LEVELS", "PRIMARY", "TAVTIGIAN_2018", "SENSITIVITY_PRIORS", "ThresholdSet",
    "posterior_from_lr", "posterior_table", "scaled_thresholds", "strength_of",
    "PEJAVER_2022_POSTERIOR_PRINTED",
]
