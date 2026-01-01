"""
Meta-learning warm-start (v10.1-A1)
- Transfer matrix com priors por target
- Regressão Bayesiana simplificada para taxas de aprendizado
"""
from __future__ import annotations
from typing import Dict, Any
import numpy as np

META_DB: Dict[str, Dict[str, float]] = {
    "EGFR": {"lr": 0.01, "n_layers": 3, "n_qubits": 12},
    "HIV": {"lr": 0.008, "n_layers": 2, "n_qubits": 10},
    "Malaria": {"lr": 0.012, "n_layers": 3, "n_qubits": 12},
    "COVID": {"lr": 0.015, "n_layers": 4, "n_qubits": 14},
}


def meta_warm_start(target: str, base_cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Ajusta hiperparâmetros a partir de priors do META_DB.

    Usa interpolação simples para taxa de aprendizado e número de camadas.
    """
    prior = META_DB.get(target.upper()) or META_DB.get(target.capitalize()) or {}
    cfg = base_cfg.copy()
    if prior:
        cfg["lr"] = float(np.mean([cfg.get("lr", 0.01), prior["lr"]]))
        cfg["n_layers"] = int(round(np.mean([cfg.get("n_layers", 3), prior["n_layers"]])))
        cfg["n_qubits"] = int(round(np.mean([cfg.get("n_qubits", cfg.get("max_qubits", 12)), prior["n_qubits"]])))
    cfg["meta_used"] = bool(prior)
    return cfg


def bayesian_lr_update(prior_mean: float, prior_var: float, grad_var: float) -> float:
    """Atualização Bayesiana simples para learning rate sugerido."""
    posterior_var = 1 / (1 / prior_var + 1 / max(grad_var, 1e-8))
    posterior_mean = posterior_var * (prior_mean / prior_var + grad_var / max(grad_var, 1e-8))
    return float(max(posterior_mean, 1e-5))


__all__ = ["META_DB", "meta_warm_start", "bayesian_lr_update"]
