"""
Fisher-CRLB constant picker + bootstrap 10k (v10.1-A1)
- Cohen's d com bootstrap 10k
- Bound de Cramér-Rao para intervalo inferior
- Seletor de constante física (π, e, φ, ℏ, α)
"""
from __future__ import annotations
import logging
from typing import Dict, Iterable, Tuple
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


def _cohen_d(x: np.ndarray, y: np.ndarray) -> float:
    n1, n2 = len(x), len(y)
    s_pooled = np.sqrt(((n1 - 1) * np.var(x, ddof=1) + (n2 - 1) * np.var(y, ddof=1)) / (n1 + n2 - 2))
    return (np.mean(x) - np.mean(y)) / s_pooled


def crlb_effect_size(group1: np.ndarray, group2: np.ndarray,
                     n_boot: int = 10_000, confidence_level: float = 0.95) -> Dict[str, float]:
    """Cohen d com bootstrap 10k + bound de Fisher."""
    d_obs = _cohen_d(group1, group2)
    boot_result = bootstrap((group1, group2), statistic=_cohen_d,
                            n_resamples=n_boot, confidence_level=confidence_level,
                            method="percentile")
    d_lb, d_ub = boot_result.confidence_interval
    fisher_lb = d_obs - 1 / np.sqrt(len(group1) + len(group2) - 2)
    return {"d_obs": d_obs, "d_lb": d_lb, "d_ub": d_ub, "fisher_lb": fisher_lb}


def fisher_information_score(values: np.ndarray) -> float:
    """Score simples para escolher constante com maior variância informativa."""
    if len(values) < 2:
        return 0.0
    return float(np.var(values) / max(np.mean(values) ** 2 + 1e-9, 1e-9))


def choose_constant_with_max_fisher(constants: Dict[str, float] | None = None,
                                     samples: int = 64,
                                     noise_scale: float = 0.05) -> Tuple[str, float]:
    """Seleciona constante física com maior score Fisher aproximado.

    Gera amostras perturbadas por ruído multiplicativo e escolhe a que maximiza
    fisher_information_score. Retorna par (nome, valor).
    """
    constants = constants or PHYSICAL_CONSTANTS
    scores = {}
    for name, value in constants.items():
        perturbed = value * (1 + noise_scale * np.random.randn(samples))
        scores[name] = fisher_information_score(perturbed)
    best = max(scores, key=scores.get)
    logger.info("Fisher constant selected: %s (score=%.4f)", best, scores[best])
    return best, constants[best]


__all__ = ["crlb_effect_size", "choose_constant_with_max_fisher", "PHYSICAL_CONSTANTS"]
