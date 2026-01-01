"""
v10_elw.py – Entanglement-Layer-Wise (Qualis A1)
"""
from __future__ import annotations
import logging
from typing import List
import numpy as np

logger = logging.getLogger("V10-ELW")


def entangling_power(circuit, params, x_sample):
    try:
        import pennylane as qml  # type: ignore
        return qml.math.entangling_power(circuit, params, x_sample)
    except Exception as exc:  # pragma: no cover
        logger.debug("entangling_power fallback: %s", exc)
        return 0.0


def choose_elw_template(candidates: List[str], circuit_fn, params_slice, x_sample) -> str:
    entropies = []
    for c in candidates:
        try:
            val = entangling_power(lambda p, x: circuit_fn(p, x, arch=c), params_slice, x_sample)
        except Exception as exc:  # pragma: no cover
            logger.debug("ELW failed for %s: %s", c, exc)
            val = -np.inf
        entropies.append(val)
    best = candidates[int(np.argmax(entropies))] if candidates else ""
    logger.info("ELW entropies: %s –> winner: %s", dict(zip(candidates, entropies)), best)
    return best


__all__ = ["entangling_power", "choose_elw_template"]
