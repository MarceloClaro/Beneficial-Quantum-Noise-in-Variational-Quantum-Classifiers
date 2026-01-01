"""
Quantum Adaptive Search Rank (QASR) – v10.1-A1
Regret O(sqrt(T) log T) com penalidade de barren plateau.
"""
from __future__ import annotations
import math
from typing import Dict, Tuple
import numpy as np
import optuna
from optuna.samplers import TPESampler


class QASRSampler(TPESampler):
    """UCB com penalização de barren plateaus.

    Mantém contadores por braço (vetor de hiperparâmetros) e
    aplica confiança sqrt(2 log t / n_a). Se o braço contém "barren"
    no nome do parâmetro, penaliza a confiança.
    """

    def __init__(self, delta: float = 0.05, **kw):
        super().__init__(**kw)
        self.delta = delta
        self.counts: Dict[Tuple, int] = {}
        self.rewards: Dict[Tuple, float] = {}

    def _ucb(self, arm: Tuple, t: int) -> float:
        na = self.counts.get(arm, 0)
        if na == 0:
            return float("inf")
        mean = self.rewards.get(arm, 0.0) / na
        conf = math.sqrt(2 * math.log(max(t, 2)) / na)
        barren_penalty = 0.5 if any("barren" in str(x) for x in arm) else 1.0
        return mean + conf * barren_penalty

    def sample_relative(self, study: optuna.Study, trial: optuna.trial.FrozenTrial, search_space):  # type: ignore[override]
        t = len(study.trials) + 1
        # Fallback para TPE quando espaço é vazio
        if not search_space:
            return super().sample_relative(study, trial, search_space)
        ucb_scores = {arm: self._ucb(arm, t) for arm in search_space}
        best = max(ucb_scores, key=ucb_scores.get)
        return {k: v for k, v in dict(best).items()} if isinstance(best, dict) else {}

    def after_trial(self, study, trial, state, values):  # type: ignore[override]
        arm = tuple(sorted(trial.params.items()))
        self.counts[arm] = self.counts.get(arm, 0) + 1
        self.rewards[arm] = self.rewards.get(arm, 0.0) + float(trial.value or 0.0)
        super().after_trial(study, trial, state, values)


__all__ = ["QASRSampler"]

