"""
v10_psqk.py – Parameter-Shared Quantum Kernel (Qualis A1)
"""
from __future__ import annotations
import numpy as np


def psqk_circuit(params, x, n_qubits: int, n_layers: int):
    import pennylane as qml  # lazy
    for _ in range(n_layers):
        for q in range(n_qubits):
            qml.RY(params[q % len(params)], wires=q)
        for q in range(n_qubits - 1):
            qml.CNOT(wires=[q, q + 1])
    return qml.expval(qml.PauliZ(0))


def psqk_shared_weights(theta: np.ndarray, share_rate: float = 0.5):
    theta = np.asarray(theta, dtype=float)
    n = len(theta)
    k = max(1, int(n * share_rate))
    mask = np.zeros_like(theta, dtype=bool)
    mask[:k] = True
    shared_value = float(theta[:k].mean()) if k > 0 else 0.0
    theta_shared = theta.copy()
    theta_shared[mask] = shared_value
    return theta_shared, mask


__all__ = ["psqk_circuit", "psqk_shared_weights"]
