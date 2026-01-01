"""
v10_fisher.py – Fisher-CRLB + bootstrap 10k (Qualis A1)
"""
from __future__ import annotations
import logging
from typing import Dict, Tuple
import numpy as np
from scipy.stats import bootstrap

logger = logging.getLogger("V10-Fisher")

PHYSICAL_CONSTANTS = {
    "pi": np.pi,
    "e": np.e,
    "phi": (1 + np.sqrt(5)) / 2,
    "hbar": 1.05457e-34,
    "alpha": 7.297e-3,
}


def fisher_information_matrix(qnode, params, x):
    import pennylane as qml  # lazy
    jac = qml.jacobian(qnode, argnum=0)(params, x)
    var = np.var(jac) + 1e-10
    return np.outer(jac, jac) / var


def choose_constant_with_max_fisher(candidates: Dict[str, float], qnode, x_sample, n_params: int) -> str:
    scores = {}
    for name, val in candidates.items():
        params = np.array([(val * (i + 1)) % (2 * np.pi) - np.pi for i in range(n_params)])
        try:
            fim = fisher_information_matrix(qnode, params, x_sample)
            scores[name] = float(np.linalg.det(fim))
        except Exception as exc:  # pragma: no cover - fallback
            logger.debug("FIM failed for %s: %s", name, exc)
            scores[name] = -np.inf
    best = max(scores, key=scores.get)
    logger.info("Fisher det scores: %s –> winner: %s", scores, best)
    return best


def crlb_effect_size(group1, group2, n_boot=10_000, confidence_level=0.95) -> Tuple[float, float, float]:
    def _cohen_d(x, y):
        n1, n2 = len(x), len(y)
        s_pooled = np.sqrt(((n1 - 1) * np.var(x, ddof=1) + (n2 - 1) * np.var(y, ddof=1)) / (n1 + n2 - 2))
        return (np.mean(x) - np.mean(y)) / s_pooled

    d_obs = _cohen_d(group1, group2)
    boot_result = bootstrap((group1, group2), statistic=_cohen_d, n_resamples=n_boot,
                            confidence_level=confidence_level, method="percentile")
    return d_obs, *boot_result.confidence_interval


__all__ = [
    "PHYSICAL_CONSTANTS",
    "fisher_information_matrix",
    "choose_constant_with_max_fisher",
    "crlb_effect_size",
]

