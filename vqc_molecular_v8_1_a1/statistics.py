# ============================================================================
# statistics.py
# Testes estatísticos múltiplos, effect sizes, e intervalos de confiança
# ============================================================================
import numpy as np
import pandas as pd
from scipy import stats
from typing import Dict, Tuple, List
import json

def cohen_d_with_bootstrap_ci(x: np.ndarray,
                               y: np.ndarray,
                               ci: float = 0.95,
                               n_bootstrap: int = 10000) -> Dict:
    """
    Calcula Cohen d com intervalo de confiança via bootstrap.

    Args:
        x: Array VQC AUCs
        y: Array baseline AUCs
        ci: Confiança (default 0.95)
        n_bootstrap: Iterações bootstrap (default 10000)

    Returns:
        Dict com d, CI, e interpretação
    """
    x = np.asarray(x)
    y = np.asarray(y)

    # Cohen d observado
    pooled_std = np.sqrt(((len(x)-1)*np.var(x, ddof=1) + (len(y)-1)*np.var(y, ddof=1)) /
                         (len(x) + len(y) - 2))
    d_obs = (np.mean(x) - np.mean(y)) / pooled_std if pooled_std > 0 else 0

    # Bootstrap
    d_boot = []
    np.random.seed(42)
    for _ in range(n_bootstrap):
        x_boot = np.random.choice(x, size=len(x), replace=True)
        y_boot = np.random.choice(y, size=len(y), replace=True)
        pooled_std_boot = np.sqrt(((len(x_boot)-1)*np.var(x_boot, ddof=1) +
                                    (len(y_boot)-1)*np.var(y_boot, ddof=1)) /
                                   (len(x_boot) + len(y_boot) - 2))
        if pooled_std_boot > 0:
            d_boot.append((np.mean(x_boot) - np.mean(y_boot)) / pooled_std_boot)
        else:
            d_boot.append(0)

    # Intervalo de confiança (percentil)
    alpha_ci = (1 - ci) / 2 * 100
    ci_low = np.percentile(d_boot, alpha_ci)
    ci_high = np.percentile(d_boot, 100 - alpha_ci)

    # Interpretação
    if abs(d_obs) < 0.2:
        interpretation = "negligible"
    elif abs(d_obs) < 0.5:
        interpretation = "small"
    elif abs(d_obs) < 0.8:
        interpretation = "medium"
    else:
        interpretation = "large"

    return {
        "cohen_d": float(d_obs),
        "ci_low": float(ci_low),
        "ci_high": float(ci_high),
        "ci_level": ci,
        "n_bootstrap": n_bootstrap,
        "interpretation": interpretation,
        "significant_at_0_05": ci_low > 0  # CI não inclui zero
    }


def hedges_g(x: np.ndarray, y: np.ndarray) -> float:
    """Calcula Hedges g (versão não-viesada de Cohen d)."""
    n1, n2 = len(x), len(y)
    pooled_std = np.sqrt(((n1-1)*np.var(x, ddof=1) + (n2-1)*np.var(y, ddof=1)) /
                         (n1 + n2 - 2))
    d = (np.mean(x) - np.mean(y)) / pooled_std if pooled_std > 0 else 0

    # Correção de viés
    J = 1 - (3 / (4 * (n1 + n2 - 2) - 1))
    return d * J


def glass_delta(x: np.ndarray, y: np.ndarray) -> float:
    """Calcula Glass delta (usa SD do grupo controle)."""
    return (np.mean(x) - np.mean(y)) / np.std(y, ddof=1)


def bonferroni_holm_correction(pvalues: List[float],
                                alpha: float = 0.05) -> Dict:
    """
    Aplica correção de Bonferroni-Holm para testes múltiplos.

    Args:
        pvalues: Lista de p-values
        alpha: Nível de significância global

    Returns:
        DataFrame com p-values originais, corrigidos, e decisão
    """
    n_tests = len(pvalues)
    sorted_idx = np.argsort(pvalues)
    sorted_pvals = np.asarray(pvalues)[sorted_idx]

    # Bonferroni-Holm: cada teste i é testado em α/(n-i+1)
    critical_values = alpha / np.arange(n_tests, 0, -1)

    # Descobrir ponto de parada (primeira que não rejeita)
    reject = np.ones(n_tests, dtype=bool)
    for i, (p, cv) in enumerate(zip(sorted_pvals, critical_values)):
        if p > cv:
            reject[i:] = False
            break

    # Reordenar para índice original
    reject_original = np.zeros(n_tests, dtype=bool)
    reject_original[sorted_idx] = reject

    adjusted_pvals = np.minimum(1.0, sorted_pvals * np.arange(n_tests, 0, -1))
    adjusted_original = np.zeros(n_tests)
    adjusted_original[sorted_idx] = adjusted_pvals

    return {
        "n_tests": n_tests,
        "original_pvalues": pvalues,
        "adjusted_pvalues": list(adjusted_original),
        "critical_values": list(critical_values),
        "reject_H0": list(reject_original),
        "n_significant": int(np.sum(reject_original)),
        "fdr_control": "strong (Bonferroni-Holm)"
    }


def fdr_benjamini_hochberg(pvalues: List[float],
                           fdr_level: float = 0.05) -> Dict:
    """
    Controle FDR (False Discovery Rate) via Benjamini-Hochberg.

    Args:
        pvalues: Lista de p-values
        fdr_level: FDR target (default 0.05)

    Returns:
        Dict com resultados
    """
    n_tests = len(pvalues)
    sorted_idx = np.argsort(pvalues)
    sorted_pvals = np.asarray(pvalues)[sorted_idx]

    # BH: p(i) ≤ (i/m) * FDR
    rank = np.arange(1, n_tests + 1)
    threshold = rank / n_tests * fdr_level

    # Encontrar maior i tal que p(i) ≤ threshold(i)
    reject_mask = sorted_pvals <= threshold
    if not np.any(reject_mask):
        largest_i = 0
    else:
        largest_i = np.max(np.where(reject_mask)[0]) + 1

    reject = np.zeros(n_tests, dtype=bool)
    reject[:largest_i] = True

    # Reordenar
    reject_original = np.zeros(n_tests, dtype=bool)
    reject_original[sorted_idx] = reject

    return {
        "n_tests": n_tests,
        "fdr_level": fdr_level,
        "n_rejections": int(np.sum(reject_original)),
        "reject_H0": list(reject_original),
        "method": "Benjamini-Hochberg",
        "control_type": "FDR (less conservative than FWER)"
    }


def ttest_with_correction(vqc_aucs: List[float],
                          baseline_aucs: List[float],
                          method: str = "bonferroni_holm",
                          alpha: float = 0.05) -> Dict:
    """
    Teste t independente com correção para múltiplos testes.

    Args:
        vqc_aucs: Lista de AUCs (VQC)
        baseline_aucs: Lista de AUCs (baseline)
        method: 'bonferroni_holm' ou 'benjamini_hochberg'
        alpha: Nível de significância

    Returns:
        Dict com resultados completos
    """
    x = np.asarray(vqc_aucs)
    y = np.asarray(baseline_aucs)

    # Teste t bilateral
    t_stat, p_bilateral = stats.ttest_ind(x, y)

    # Teste t unilateral (superior: VQC > baseline)
    p_unilateral = stats.ttest_ind(x, y).pvalue / 2 if t_stat > 0 else 1 - stats.ttest_ind(x, y).pvalue / 2

    # Effect sizes
    d_info = cohen_d_with_bootstrap_ci(x, y, ci=0.95, n_bootstrap=10000)
    g = hedges_g(x, y)
    delta = glass_delta(x, y)

    # Correção múltipla
    pvals = [p_bilateral]  # aqui teríamos mais p-values em análises múltiplas
    if method == "bonferroni_holm":
        correction = bonferroni_holm_correction(pvals, alpha=alpha)
    else:
        correction = fdr_benjamini_hochberg(pvals, fdr_level=alpha)

    return {
        "comparison": "VQC vs Baseline",
        "n_vqc": len(x),
        "n_baseline": len(y),
        "mean_vqc": float(np.mean(x)),
        "mean_baseline": float(np.mean(y)),
        "mean_difference": float(np.mean(x) - np.mean(y)),
        "std_vqc": float(np.std(x, ddof=1)),
        "std_baseline": float(np.std(y, ddof=1)),
        "ttest": {
            "t_statistic": float(t_stat),
            "p_bilateral": float(p_bilateral),
            "p_unilateral": float(p_unilateral),
        },
        "effect_sizes": {
            "cohen_d": d_info["cohen_d"],
            "cohen_d_ci": [d_info["ci_low"], d_info["ci_high"]],
            "interpretation": d_info["interpretation"],
            "hedges_g": float(g),
            "glass_delta": float(delta)
        },
        "multiple_comparison": correction,
        "conclusion": "SUPERIOR" if correction["reject_H0"][0] and t_stat > 0 else "NON-SIGNIFICANT"
    }


def comprehensive_statistics_report(df: pd.DataFrame,
                                    vqc_col: str = "vqc_auc",
                                    baseline_col: str = "baseline_auc",
                                    group_col: str = None) -> str:
    """
    Gera relatório completo de estatísticas.

    Args:
        df: DataFrame com colunas AUC
        vqc_col: Nome coluna VQC
        baseline_col: Nome coluna baseline
        group_col: Coluna para subgrupos (ex: "target" ou "fold")

    Returns:
        String formatada para relatório
    """
    report = "\n" + "="*80 + "\n"
    report += "RELATÓRIO ESTATÍSTICO COMPLETO\n"
    report += "="*80 + "\n\n"

    if group_col is None:
        # Análise global
        results = ttest_with_correction(
            df[vqc_col].tolist(),
            df[baseline_col].tolist()
        )
        report += f"ANÁLISE GLOBAL:\n"
        report += f"  N total: {results['n_vqc']}\n"
        report += f"  VQC: {results['mean_vqc']:.4f} ± {results['std_vqc']:.4f}\n"
        report += f"  Baseline: {results['mean_baseline']:.4f} ± {results['std_baseline']:.4f}\n"
        report += f"  Diferença: {results['mean_difference']:.4f}\n"
        report += f"  Cohen d: {results['effect_sizes']['cohen_d']:.3f} [{results['effect_sizes']['cohen_d_ci'][0]:.3f}, {results['effect_sizes']['cohen_d_ci'][1]:.3f}]\n"
        report += f"  p-value: {results['ttest']['p_bilateral']:.6f}\n"
        report += f"  Conclusão: {results['conclusion']}\n"
    else:
        # Análise por grupo
        for group in df[group_col].unique():
            subset = df[df[group_col] == group]
            results = ttest_with_correction(
                subset[vqc_col].tolist(),
                subset[baseline_col].tolist()
            )
            report += f"\n{group}:\n"
            report += f"  N: {len(subset)}\n"
            report += f"  VQC: {results['mean_vqc']:.4f} ± {results['std_vqc']:.4f}\n"
            report += f"  Cohen d: {results['effect_sizes']['cohen_d']:.3f}\n"
            report += f"  p-value: {results['ttest']['p_bilateral']:.6f}\n"

    report += "\n" + "="*80 + "\n"
    return report


if __name__ == "__main__":
    # Exemplo
    vqc = [0.92, 0.91, 0.93, 0.90, 0.94]
    baseline = [0.88, 0.87, 0.89, 0.86, 0.90]

    print(comprehensive_statistics_report(
        pd.DataFrame({"vqc_auc": vqc, "baseline_auc": baseline})
    ))

    results = ttest_with_correction(vqc, baseline)
    print(json.dumps(results, indent=2))
