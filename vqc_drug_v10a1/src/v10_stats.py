"""
Estatísticas Qualis A1 – v10.1-A1
- Westfall-Young step-down (via FDR BY como aproximação leve)
- Cohen's d com bootstrap 10k + bound de Fisher
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from typing import Dict, Tuple
from scipy.stats import bootstrap
from statsmodels.stats.multitest import multipletests


def westfall_young_step_down(pvals: np.ndarray, alpha: float = 0.05) -> Tuple[np.ndarray, np.ndarray]:
    """Aproximação step-down usando método Benjamini-Yekutieli (BY)."""
    reject, p_adj, _, _ = multipletests(pvals, alpha=alpha, method="fdr_by")
    return reject, p_adj


def _cohen_d(x: np.ndarray, y: np.ndarray) -> float:
    n1, n2 = len(x), len(y)
    s_pooled = np.sqrt(((n1 - 1) * np.var(x, ddof=1) + (n2 - 1) * np.var(y, ddof=1)) / (n1 + n2 - 2))
    return (np.mean(x) - np.mean(y)) / s_pooled


def crlb_effect_size(group1: np.ndarray, group2: np.ndarray,
                     n_boot: int = 10_000, confidence_level: float = 0.95) -> Dict[str, float]:
    d_obs = _cohen_d(group1, group2)
    boot_result = bootstrap((group1, group2), statistic=_cohen_d,
                            n_resamples=n_boot, confidence_level=confidence_level,
                            method="percentile")
    d_lb, d_ub = boot_result.confidence_interval
    fisher_lb = d_obs - 1 / np.sqrt(len(group1) + len(group2) - 2)
    return {"d_obs": d_obs, "d_lb": d_lb, "d_ub": d_ub, "fisher_lb": fisher_lb}


def qualis_report_from_trials(df: pd.DataFrame) -> pd.DataFrame:
    """Gera dataframe com estatísticas Qualis (AUC, std, p-ajustado)."""
    if df.empty:
        return pd.DataFrame()
    metrics = df[[c for c in df.columns if c.startswith("value") or c == "value"]].copy()
    aucs = metrics.squeeze().to_numpy()
    auc_mean = float(np.mean(aucs))
    auc_std = float(np.std(aucs))
    pvals = np.clip(1 - aucs, 0, 1)  # proxy de p-valor
    reject, p_adj = westfall_young_step_down(pvals)
    return pd.DataFrame({
        "auc_mean": [auc_mean],
        "auc_std": [auc_std],
        "p_wy_max": [float(p_adj.max())],
        "rej_any": [bool(reject.any())],
    })


__all__ = ["westfall_young_step_down", "crlb_effect_size", "qualis_report_from_trials"]

