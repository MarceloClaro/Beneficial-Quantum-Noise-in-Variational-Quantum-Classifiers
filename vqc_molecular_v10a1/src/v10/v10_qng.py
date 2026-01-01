"""
v10_qng.py – Exact QNG with block-diagonal metric (Qualis A1)
"""
from __future__ import annotations
import numpy as np


def qng_metric_block_diag(jacobian: np.ndarray, block_size: int = 2) -> np.ndarray:
    jac = np.atleast_2d(jacobian)
    n_params = jac.shape[1]
    blocks = []
    for start in range(0, n_params, block_size):
        end = min(start + block_size, n_params)
        j_block = jac[:, start:end]
        g_block = j_block.T @ j_block
        blocks.append(g_block)
    return np.block([[blocks[i] if i == j else np.zeros_like(blocks[min(i, j)])
                      for j in range(len(blocks))] for i in range(len(blocks))])


def condition_number(metric: np.ndarray) -> float:
    vals = np.linalg.svd(metric, compute_uv=False)
    if len(vals) == 0 or vals.min() == 0:
        return float("inf")
    return float(vals.max() / vals.min())


__all__ = ["qng_metric_block_diag", "condition_number"]
