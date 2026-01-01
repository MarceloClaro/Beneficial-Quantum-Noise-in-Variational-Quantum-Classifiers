"""
Parameter-Sharing Quantum Kernel (PS-QK) – v10.1-A1
- Reduz parâmetros em ~50% via compartilhamento
- Ablation friendly
"""
from __future__ import annotations
import numpy as np
from typing import Tuple


def psqk_shared_weights(theta: np.ndarray, share_rate: float = 0.5) -> Tuple[np.ndarray, np.ndarray]:
    """Aplica compartilhamento de parâmetros.

    Retorna (theta_compartilhado, mask_aplicada).
    """
    theta = np.asarray(theta, dtype=float)
    n = len(theta)
    k = max(1, int(n * share_rate))
    mask = np.zeros_like(theta, dtype=bool)
    mask[:k] = True
    shared_value = float(theta[:k].mean()) if k > 0 else 0.0
    theta_shared = theta.copy()
    theta_shared[mask] = shared_value
    return theta_shared, mask


__all__ = ["psqk_shared_weights"]
