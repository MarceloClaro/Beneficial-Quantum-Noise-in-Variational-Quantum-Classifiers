"""
Zero-Noise Extrapolation on-the-fly (v10.1-A1)
- Richardson 3 níveis
- Gatilho via curva de aprendizado
"""
from __future__ import annotations
import numpy as np
from typing import Iterable, Tuple


def richardson_extrapolate(values: Iterable[float], scale_factors: Iterable[float]) -> float:
    """Aplica Richardson (3 pontos) para estimar valor em ruído zero.

    Espera-se len(values)==len(scale_factors)==3.
    """
    vals = list(values)
    scales = list(scale_factors)
    if len(vals) != 3 or len(scales) != 3:
        raise ValueError("Richardson requer exatamente 3 pontos")
    # Sistema linear simples
    a1 = (scales[2] * vals[0] - scales[0] * vals[2]) / (scales[2] - scales[0])
    a2 = (scales[2] * vals[1] - scales[1] * vals[2]) / (scales[2] - scales[1])
    return float((a1 + a2) / 2)


def should_trigger_zne(loss_history: Iterable[float], window: int = 5, tol: float = 1e-4) -> bool:
    """Dispara ZNE se curva de perda estabilizar (derivada baixa)."""
    losses = list(loss_history)
    if len(losses) < window + 1:
        return False
    recent = losses[-window:]
    grad = np.abs(np.diff(recent)).mean()
    return bool(grad < tol)


def generate_zne_curve(fn, x, noise_scales=(1.0, 1.5, 2.0)) -> Tuple[np.ndarray, float]:
    """Avalia função fn em 3 níveis de ruído simulados (multiplica x)."""
    vals = []
    for s in noise_scales:
        vals.append(float(fn(x * s)))
    zero_noise = richardson_extrapolate(vals, noise_scales)
    return np.array(vals), zero_noise


__all__ = ["richardson_extrapolate", "should_trigger_zne", "generate_zne_curve"]

