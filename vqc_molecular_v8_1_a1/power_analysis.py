# ============================================================================
# power_analysis.py
# Análise de poder pré-experimento para determinação de tamanho amostral
# ============================================================================
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Dict, List
from scipy import stats

def required_sample_size(effect_size: float = 0.35,
                         alpha: float = 0.05,
                         power: float = 0.8,
                         ratio: float = 1.0) -> int:
    """
    Calcula n por grupo usando t-test independente.

    Args:
        effect_size: Cohen d (default 0.35 = small-to-medium effect)
        alpha: Nível de significância bilateral (default 0.05)
        power: Potência desejada (1 - β, default 0.8)
        ratio: Razão n2/n1 (default 1 = grupos iguais)

    Returns:
        Tamanho amostral por grupo (ceil)
    """
    # Fórmula de Neyman
    # n = [(1 + 1/ratio) * (z_alpha + z_beta)^2] / effect_size^2
    z_alpha = stats.norm.ppf(1 - alpha/2)  # bilateral
    z_beta = stats.norm.ppf(power)
    numerator = (1 + 1/ratio) * (z_alpha + z_beta) ** 2
    denominator = effect_size ** 2
    n = numerator / denominator
    return int(np.ceil(n))


def power_curve(effect_sizes: List[float],
                alpha: float = 0.05,
                n_per_group: List[int] = None) -> Dict:
    """
    Calcula curva de poder (power vs effect size ou power vs n).

    Args:
        effect_sizes: Lista de effect sizes (Cohen d)
        alpha: Nível de significância
        n_per_group: Lista de tamanhos amostrais (se None, calcula para d=0.35)

    Returns:
        Dict com resultados
    """
    if n_per_group is None:
        n_per_group = [required_sample_size(d, alpha=alpha) for d in effect_sizes]

    results = {
        "effect_sizes": effect_sizes,
        "n_per_group": n_per_group,
        "powers": [],
        "alpha": alpha
    }

    z_alpha = stats.norm.ppf(1 - alpha/2)

    for d, n in zip(effect_sizes, n_per_group):
        # Potência = Φ(|d|*sqrt(n/2) - z_α)
        noncentrality = np.abs(d) * np.sqrt(n / 2)
        power = stats.norm.cdf(noncentrality - z_alpha)
        results["powers"].append(power)

    return results


def plot_power_curve(effect_sizes: List[float] = None,
                     alpha: float = 0.05,
                     n_values: List[int] = None,
                     title: str = "Power Analysis",
                     figsize: Tuple = (8, 6),
                     dpi: int = 600,
                     output_file: str = None) -> str:
    """
    Plota curva de poder em alta resolução (600 dpi).

    Args:
        effect_sizes: Lista de Cohen d (default: 0.2 a 0.8)
        alpha: Nível de significância
        n_values: Lista de tamanhos amostrais (para ficar de outro lado)
        title: Título da figura
        figsize: Tamanho (inches)
        dpi: Resolução
        output_file: Arquivo de saída (default: power_analysis.png)

    Returns:
        Caminho do arquivo salvo
    """
    if output_file is None:
        output_file = "power_analysis.png"

    if effect_sizes is None:
        effect_sizes = np.linspace(0.1, 1.0, 20)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # Curva 1: Power vs Effect Size (n fixo)
    if n_values is None:
        n_values = [50, 100, 200, 500]

    for n in n_values:
        results = power_curve([d for d in effect_sizes], alpha=alpha, n_per_group=[n]*len(effect_sizes))
        ax.plot(results["effect_sizes"], results["powers"],
                label=f"n={n}", lw=2, marker='o', markersize=4, alpha=0.8)

    # Linha de referência (power = 0.8)
    ax.axhline(0.8, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='Power = 0.8')
    ax.axhline(0.9, color='orange', linestyle='--', linewidth=1, alpha=0.5, label='Power = 0.9')

    # Formatação
    ax.set_xlabel("Effect Size (Cohen d)", fontsize=12, fontweight='bold')
    ax.set_ylabel("Power (1 - β)", fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_ylim([0, 1.05])
    ax.set_xlim([0, max(effect_sizes)])
    ax.grid(True, alpha=0.3)
    ax.legend(loc='lower right', framealpha=0.9)

    # Adicionar texto com configuração
    textstr = f"α = {alpha}\nHypothesis: two-tailed"
    ax.text(0.98, 0.05, textstr, transform=ax.transAxes,
            fontsize=10, verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close()

    print(f"✅ Curva de poder salva: {output_file}")
    return output_file


def sample_size_table(effect_sizes: List[float] = None,
                      alpha: float = 0.05,
                      power: float = 0.8) -> str:
    """
    Gera tabela de tamanho amostral recomendado.

    Returns:
        String formatada como tabela ASCII
    """
    if effect_sizes is None:
        effect_sizes = [0.20, 0.35, 0.50, 0.65, 0.80]  # Cohen: small, small-med, med, med-large, large

    # Header
    table = "\n" + "="*80 + "\n"
    table += f"SAMPLE SIZE TABLE (α={alpha}, Power={power})\n"
    table += "="*80 + "\n"
    table += f"{'Cohen d':<12} {'Description':<20} {'n per group':<15} {'Total N':<15}\n"
    table += "-"*80 + "\n"

    descriptions = {
        0.20: "Small effect",
        0.35: "Small-to-medium",
        0.50: "Medium effect",
        0.65: "Medium-to-large",
        0.80: "Large effect",
    }

    for d in effect_sizes:
        n = required_sample_size(d, alpha=alpha, power=power)
        desc = descriptions.get(d, "Custom")
        table += f"{d:<12.2f} {desc:<20} {n:<15} {n*2:<15}\n"

    table += "="*80 + "\n\n"
    return table


def power_summary_json(effect_size: float = 0.35,
                       alpha: float = 0.05,
                       power: float = 0.8) -> Dict:
    """Retorna resumo de análise de poder em JSON."""
    n = required_sample_size(effect_size, alpha=alpha, power=power)
    return {
        "analysis_type": "independent_samples_ttest",
        "hypothesis": "two_tailed",
        "effect_size": effect_size,
        "effect_size_interpretation": "small_to_medium",
        "alpha": alpha,
        "desired_power": power,
        "n_per_group": n,
        "total_n": n * 2,
        "recommendation": f"Planeje coletar pelo menos {n} amostras por grupo (total {n*2})"
    }


if __name__ == "__main__":
    # Exemplo de uso
    print("\n" + "="*80)
    print("POWER ANALYSIS PRÉ-EXPERIMENTO")
    print("="*80)

    # Tabela de tamanho amostral
    print(sample_size_table(effect_sizes=[0.20, 0.35, 0.50, 0.80], alpha=0.05, power=0.8))

    # Calcular n para efeito observado
    effect = 0.35
    n = required_sample_size(effect, alpha=0.05, power=0.8)
    print(f"Para Cohen d = {effect}:")
    print(f"  → n por grupo = {n}")
    print(f"  → Total N = {n*2}")

    # Plotar curva de poder
    plot_power_curve(
        effect_sizes=np.linspace(0.1, 1.0, 20),
        alpha=0.05,
        n_values=[50, 100, 200, 500],
        output_file="power_analysis_example.png"
    )

    # JSON
    import json
    summary = power_summary_json(0.35, 0.05, 0.8)
    print("\n" + json.dumps(summary, indent=2))
