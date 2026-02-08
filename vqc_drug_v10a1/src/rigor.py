"""Ferramentas de validação estatística rigorosa e reproduzível para estudos A1."""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Any

import numpy as np
import pandas as pd
from scipy import stats


@dataclass
class ValidationSummary:
    n_trials: int
    auc_mean: float
    auc_median: float
    auc_std: float
    ci95_low: float
    ci95_high: float
    pvalue_vs_random: float
    effect_size_vs_random: float
    permutation_pvalue: float
    stability_cv: float
    seed: int


def _extract_trial_auc(df: pd.DataFrame) -> np.ndarray:
    if "value" not in df.columns:
        raise ValueError("DataFrame de trials deve conter coluna 'value' (AUC por trial).")
    auc = pd.to_numeric(df["value"], errors="coerce").dropna().to_numpy(dtype=float)
    if auc.size < 3:
        raise ValueError("São necessários ao menos 3 trials válidos para validação estatística.")
    return auc


def _bootstrap_ci(values: np.ndarray, confidence: float = 0.95, n_resamples: int = 5000, seed: int = 42) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, values.size, size=(n_resamples, values.size))
    means = values[idx].mean(axis=1)
    lo = float(np.quantile(means, (1 - confidence) / 2))
    hi = float(np.quantile(means, 1 - (1 - confidence) / 2))
    return lo, hi


def _permutation_pvalue(values: np.ndarray, baseline: float, n_permutations: int = 5000, seed: int = 42) -> float:
    centered = values - baseline
    observed = float(centered.mean())
    rng = np.random.default_rng(seed)
    signs = rng.choice([-1.0, 1.0], size=(n_permutations, values.size))
    permuted = (centered * signs).mean(axis=1)
    pvalue = (np.abs(permuted) >= abs(observed)).mean()
    return float(max(pvalue, 1.0 / n_permutations))


def build_validation_report(
    df: pd.DataFrame,
    *,
    alpha: float = 0.05,
    random_baseline_auc: float = 0.5,
    n_permutations: int = 5000,
    seed: int = 42,
) -> Dict[str, Any]:
    """Gera relatório estatístico completo para publicação reproduzível."""
    auc = _extract_trial_auc(df)

    ci_low, ci_high = _bootstrap_ci(auc, confidence=1 - alpha, seed=seed)
    t_stat, pvalue = stats.ttest_1samp(auc, popmean=random_baseline_auc, alternative="greater")

    pooled_std = auc.std(ddof=1)
    effect_size = (auc.mean() - random_baseline_auc) / pooled_std if pooled_std > 0 else 0.0

    summary = ValidationSummary(
        n_trials=int(auc.size),
        auc_mean=float(auc.mean()),
        auc_median=float(np.median(auc)),
        auc_std=float(auc.std(ddof=1)),
        ci95_low=ci_low,
        ci95_high=ci_high,
        pvalue_vs_random=float(pvalue),
        effect_size_vs_random=float(effect_size),
        permutation_pvalue=_permutation_pvalue(auc, random_baseline_auc, n_permutations, seed),
        stability_cv=float(auc.std(ddof=1) / (auc.mean() + 1e-12)),
        seed=seed,
    )

    return {
        "summary": asdict(summary),
        "hypotheses": {
            "h0": f"AUC média <= {random_baseline_auc:.2f}",
            "h1": f"AUC média > {random_baseline_auc:.2f}",
            "alpha": alpha,
            "ttest_statistic": float(t_stat),
        },
        "reproducibility": {
            "n_permutations": n_permutations,
            "seed": seed,
            "deterministic": True,
        },
        "quality_gates": {
            "gate_pvalue": bool(summary.pvalue_vs_random < alpha),
            "gate_permutation": bool(summary.permutation_pvalue < alpha),
            "gate_ci_above_random": bool(summary.ci95_low > random_baseline_auc),
            "gate_stability_cv_lt_10pct": bool(summary.stability_cv < 0.10),
        },
    }


def write_validation_artifacts(report: Dict[str, Any], out_dir: str | Path = ".") -> None:
    """Persiste artefatos JSON + Markdown para rastreabilidade em submissão A1."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    json_path = out_dir / "02_validation_report.json"
    md_path = out_dir / "02_validation_report.md"

    with json_path.open("w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    s = report["summary"]
    qg = report["quality_gates"]
    lines = [
        "# Relatório de Validação Estatística (A1)",
        "",
        "## Resumo",
        f"- Trials válidos: **{s['n_trials']}**",
        f"- AUC média: **{s['auc_mean']:.4f}** (IC95% {s['ci95_low']:.4f}, {s['ci95_high']:.4f})",
        f"- p-valor (teste unilateral vs aleatório): **{s['pvalue_vs_random']:.3e}**",
        f"- p-valor (permutação com sinal): **{s['permutation_pvalue']:.3e}**",
        f"- Tamanho de efeito (Cohen d): **{s['effect_size_vs_random']:.4f}**",
        f"- Estabilidade (CV): **{100*s['stability_cv']:.2f}%**",
        "",
        "## Gates de Qualidade",
        f"- p-valor < alpha: **{qg['gate_pvalue']}**",
        f"- permutação < alpha: **{qg['gate_permutation']}**",
        f"- IC95% totalmente acima do baseline: **{qg['gate_ci_above_random']}**",
        f"- CV < 10%: **{qg['gate_stability_cv_lt_10pct']}**",
    ]
    md_path.write_text("\n".join(lines), encoding="utf-8")


__all__ = ["build_validation_report", "write_validation_artifacts"]
