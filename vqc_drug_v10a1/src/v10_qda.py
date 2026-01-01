"""
Quantum Data Augmentation (QDA) – v10.1-A1
Compatível com classical shadows simples.
"""
from __future__ import annotations
import numpy as np
from typing import Tuple


def qda_augment_shadows(shadows: np.ndarray, n_aug: int = 2, noise: float = 0.01) -> np.ndarray:
    """Gera amostras augmentadas a partir de classical shadows (mock)."""
    shadows = np.asarray(shadows, dtype=float)
    aug = [shadows]
    for _ in range(n_aug):
        perturb = shadows + noise * np.random.randn(*shadows.shape)
        aug.append(perturb)
    return np.vstack(aug)


def shadow_statistics(shadows: np.ndarray) -> Tuple[float, float]:
    shadows = np.asarray(shadows, dtype=float)
    return float(shadows.mean()), float(shadows.std())


__all__ = ["qda_augment_shadows", "shadow_statistics"]
