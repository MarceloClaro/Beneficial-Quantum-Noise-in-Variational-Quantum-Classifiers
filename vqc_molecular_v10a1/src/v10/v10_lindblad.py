"""
v10_lindblad.py – Lindblad-optimal schedule (Qualis A1)
"""
from __future__ import annotations
import logging
from typing import Dict, Tuple
import numpy as np

logger = logging.getLogger("V10-Lindblad")

DEFAULT_T1_T2: Dict[str, Tuple[float, float]] = {
    "ibmq_jakarta": (50e-6, 70e-6),
    "default": (50e-6, 70e-6),
}


class LindbladOptimalSchedule:
    def __init__(self, T1: float = 50e-6, T2: float = 70e-6, gate_time: float = 20e-9):
        self.gamma_T1 = 1 - np.exp(-gate_time / T1)
        self.gamma_T2 = 1 - np.exp(-gate_time / T2)
        self.gate_time = gate_time

    def level(self, epoch: int, n_epochs: int) -> float:
        t = epoch / max(1, n_epochs - 1)
        return float(self.gamma_T1 * np.exp(-t) + self.gamma_T2 * np.exp(-2 * t))

    def lookup_table(self, n_epochs: int) -> np.ndarray:
        return np.array([self.level(e, n_epochs) for e in range(n_epochs)])

    def summary(self, n_epochs: int) -> Dict[str, float]:
        lut = self.lookup_table(n_epochs)
        return {
            "min": float(lut.min()),
            "max": float(lut.max()),
            "mean": float(lut.mean()),
            "gamma_T1": float(self.gamma_T1),
            "gamma_T2": float(self.gamma_T2),
        }


def load_backend_schedule(backend: str = "default", gate_time: float = 20e-9) -> "LindbladOptimalSchedule":
    T1, T2 = DEFAULT_T1_T2.get(backend, DEFAULT_T1_T2["default"])
    return LindbladOptimalSchedule(T1=T1, T2=T2, gate_time=gate_time)


__all__ = ["LindbladOptimalSchedule", "load_backend_schedule", "DEFAULT_T1_T2"]

