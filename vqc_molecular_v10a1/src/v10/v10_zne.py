"""
v10_zne.py – ZNE 3-level Richardson (Qualis A1)
"""
from __future__ import annotations
import logging
from typing import Tuple
import numpy as np

logger = logging.getLogger("V10-ZNE")


def richardson_extrapolate(values: Tuple[float, float, float], scales: Tuple[float, float, float]) -> float:
    y1, y2, y3 = values
    s1, s2, s3 = scales
    return y1 + (y2 - y1) * (s1 / (s2 - s1)) + (y3 - y2) * (s2 / (s3 - s2))


def zne_on_fly(loss_noisy: float, loss_1_5x: float, loss_2x: float) -> float:
    return richardson_extrapolate((loss_noisy, loss_1_5x, loss_2x), (1.0, 1.5, 2.0))


def should_trigger_zne(loss_history, window: int = 5, tol: float = 1e-4) -> bool:
    losses = list(loss_history)
    if len(losses) < window + 1:
        return False
    recent = losses[-window:]
    grad = np.abs(np.diff(recent)).mean()
    return bool(grad < tol)


__all__ = ["richardson_extrapolate", "zne_on_fly", "should_trigger_zne"]
