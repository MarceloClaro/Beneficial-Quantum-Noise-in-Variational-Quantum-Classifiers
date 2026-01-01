# ============================================================================
# figures.py
# Geração de figuras em alta resolução (600 dpi) para publicação
# ============================================================================
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Optional
from sklearn.metrics import roc_curve, auc

# Estilos de publicação
plt.style.use('default')
sns.set_palette("husl")

# Formato de fontes para publicação
plt.rcParams.update({
    'font.size': 10,
    'font.family': 'serif',
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.dpi': 100,  # renderização
    'savefig.dpi': 600,  # salva em 600 dpi
    'lines.linewidth': 1.5,
    'lines.markersize': 6
})


def create_output_dir(output_dir: str = "04_figuras_publicacao") -> str:
    """Cria diretório para figuras."""
    os.makedirs(output_dir, exist_ok=True)
    return output_dir


def fig_power_curve(effect_sizes: List[float] = None,
                    n_per_group: List[int] = None,
                    alpha: float = 0.05,
                    output_dir: str = "04_figuras_publicacao",
                    filename: str = "fig1_power_curve.png") -> str:
    """
    Figura 1: Curva de poder (Power vs Effect Size).

    Args:
        effect_sizes: Lista de Cohen d (default 0.1 a 1.0)
        n_per_group: Tamanhos amostrais para diferentes curvas
        alpha: Nível de significância
        output_dir: Diretório de saída
        filename: Nome do arquivo

    Returns:
        Caminho do arquivo
    """
    if effect_sizes is None:
        effect_sizes = np.linspace(0.1, 1.0, 30)

    if n_per_group is None:
        n_per_group = [50, 100, 200, 500]

    os.makedirs(output_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(4, 3.5))

    # Calcular poder para cada n
    from scipy import stats
    for n in n_per_group:
        powers = []
        for d in effect_sizes:
            z_alpha = stats.norm.ppf(1 - alpha/2)
            noncentrality = np.abs(d) * np.sqrt(n / 2)
            power = stats.norm.cdf(noncentrality - z_alpha)
            powers.append(power)

        ax.plot(effect_sizes, powers, marker='o', markersize=3,
                label=f'$n$ = {n}', linewidth=1.5, alpha=0.85)

    # Linhas de referência
    ax.axhline(0.80, color='red', linestyle='--', linewidth=1, alpha=0.5, label='Power = 0.80')
    ax.axhline(0.90, color='orange', linestyle='--', linewidth=1, alpha=0.4, label='Power = 0.90')

    # Formatação
    ax.set_xlabel('Effect Size (Cohen $d$)', fontweight='bold')
    ax.set_ylabel('Statistical Power', fontweight='bold')
    ax.set_title('Power Analysis Curves', fontweight='bold', fontsize=12)
    ax.set_ylim([0, 1.05])
    ax.set_xlim([0, 1.0])
    ax.grid(True, alpha=0.25, linestyle=':', linewidth=0.5)
    ax.legend(loc='lower right', framealpha=0.95, fontsize=9)

    # Anotação
    ax.text(0.05, 0.05, f'$\\alpha$ = {alpha} (two-tailed)',
            transform=ax.transAxes, fontsize=9,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    plt.tight_layout()
    filepath = os.path.join(output_dir, filename)
    plt.savefig(filepath, dpi=600, bbox_inches='tight', format='png')
    plt.close()

    print(f"✅ {filename} saved ({os.path.getsize(filepath)/1024:.1f} KB)")
    return filepath


def fig_roc_comparison(vqc_aucs: List[float],
                       baseline_aucs: List[float],
                       y_true: np.ndarray,
                       y_pred_vqc: np.ndarray,
                       y_pred_baseline: np.ndarray,
                       output_dir: str = "04_figuras_publicacao",
                       filename: str = "fig2_roc_comparison.png") -> str:
    """
    Figura 2: Curvas ROC (VQC vs Baseline).

    Args:
        vqc_aucs: Lista de AUCs (VQC) para múltiplos folds
        baseline_aucs: Lista de AUCs (Baseline)
        y_true: Labels verdadeiros
        y_pred_vqc: Predições VQC (probabilidades)
        y_pred_baseline: Predições Baseline
        output_dir: Diretório de saída
        filename: Nome do arquivo

    Returns:
        Caminho do arquivo
    """
    os.makedirs(output_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(5, 4.5))

    # Curvas ROC
    fpr_vqc, tpr_vqc, _ = roc_curve(y_true, y_pred_vqc)
    fpr_baseline, tpr_baseline, _ = roc_curve(y_true, y_pred_baseline)

    auc_vqc = auc(fpr_vqc, tpr_vqc)
    auc_baseline = auc(fpr_baseline, tpr_baseline)

    ax.plot(fpr_vqc, tpr_vqc, lw=2, label=f'VQC (AUC = {auc_vqc:.3f})',
            color='#1f77b4', alpha=0.85)
    ax.plot(fpr_baseline, tpr_baseline, lw=2, label=f'Baseline (AUC = {auc_baseline:.3f})',
            color='#ff7f0e', alpha=0.85)

    # Diagonal (random classifier)
    ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.3, label='Random')

    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])
    ax.set_xlabel('False Positive Rate', fontweight='bold')
    ax.set_ylabel('True Positive Rate', fontweight='bold')
    ax.set_title('ROC Curves: VQC vs Baseline', fontweight='bold', fontsize=12)
    ax.legend(loc='lower right', framealpha=0.95, fontsize=10)
    ax.grid(True, alpha=0.25, linestyle=':', linewidth=0.5)

    plt.tight_layout()
    filepath = os.path.join(output_dir, filename)
    plt.savefig(filepath, dpi=600, bbox_inches='tight', format='png')
    plt.close()

    print(f"✅ {filename} saved ({os.path.getsize(filepath)/1024:.1f} KB)")
    return filepath


def fig_forest_plot(effect_sizes: Dict[str, Dict],
                    targets: List[str],
                    output_dir: str = "04_figuras_publicacao",
                    filename: str = "fig3_forest_effect_sizes.png") -> str:
    """
    Figura 3: Forest plot com effect sizes (Cohen d) e IC 95%.

    Args:
        effect_sizes: Dict com estrutura {target: {cohen_d, ci_low, ci_high}}
        targets: Lista de targets (ex: ['EGFR', 'HIV', 'Malaria', 'COVID'])
        output_dir: Diretório de saída
        filename: Nome do arquivo

    Returns:
        Caminho do arquivo
    """
    os.makedirs(output_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(5, 4))

    # Preparar dados
    n_studies = len(targets)
    y_pos = np.arange(n_studies)

    cohen_ds = []
    ci_lows = []
    ci_highs = []

    for target in targets:
        if target in effect_sizes:
            d = effect_sizes[target].get('cohen_d', 0)
            ci_low = effect_sizes[target].get('ci_low', d)
            ci_high = effect_sizes[target].get('ci_high', d)
        else:
            d, ci_low, ci_high = 0, 0, 0

        cohen_ds.append(d)
        ci_lows.append(ci_low)
        ci_highs.append(ci_high)

    # Erros (diferença do ponto central)
    errors = [
        np.array(cohen_ds) - np.array(ci_lows),
        np.array(ci_highs) - np.array(cohen_ds)
    ]

    # Plotar pontos e ICs
    colors = plt.cm.Set2(np.linspace(0, 1, n_studies))
    ax.scatter(cohen_ds, y_pos, s=100, color=colors, zorder=3, alpha=0.8)

    for i, (d, low, high) in enumerate(zip(cohen_ds, ci_lows, ci_highs)):
        ax.plot([low, high], [i, i], color=colors[i], lw=2, alpha=0.6)
        # Marcas de IC
        ax.plot([low, low], [i-0.1, i+0.1], color=colors[i], lw=1.5)
        ax.plot([high, high], [i-0.1, i+0.1], color=colors[i], lw=1.5)

    # Linha de nulidade (d=0)
    ax.axvline(0, color='red', linestyle='--', lw=1, alpha=0.4, label='No effect')

    # Formatação
    ax.set_yticks(y_pos)
    ax.set_yticklabels(targets)
    ax.set_xlabel("Effect Size (Cohen's $d$)", fontweight='bold')
    ax.set_title('Effect Sizes with 95% Confidence Intervals', fontweight='bold', fontsize=12)
    ax.grid(True, alpha=0.2, axis='x', linestyle=':')
    ax.legend(loc='best', framealpha=0.95)

    plt.tight_layout()
    filepath = os.path.join(output_dir, filename)
    plt.savefig(filepath, dpi=600, bbox_inches='tight', format='png')
    plt.close()

    print(f"✅ {filename} saved ({os.path.getsize(filepath)/1024:.1f} KB)")
    return filepath


def fig_optuna_history(study,
                       output_dir: str = "04_figuras_publicacao",
                       filename: str = "fig4_optuna_trials.png") -> str:
    """
    Figura 4: Histórico de otimização Optuna.

    Args:
        study: Optuna Study object
        output_dir: Diretório de saída
        filename: Nome do arquivo

    Returns:
        Caminho do arquivo
    """
    os.makedirs(output_dir, exist_ok=True)

    trials_df = study.trials_dataframe()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    # Gráfico 1: Value vs Trial
    ax1.scatter(trials_df['number'], trials_df['value'], alpha=0.6, s=30, color='#1f77b4')
    ax1.plot(trials_df['number'], trials_df['value'].cummax(), 'r-', lw=2, alpha=0.7, label='Best so far')
    ax1.set_xlabel('Trial Number', fontweight='bold')
    ax1.set_ylabel('Objective Value (ROC-AUC)', fontweight='bold')
    ax1.set_title('Optimization History', fontweight='bold')
    ax1.grid(True, alpha=0.25, linestyle=':')
    ax1.legend(loc='lower right')

    # Gráfico 2: Distribuição de values
    ax2.hist(trials_df['value'], bins=20, alpha=0.7, color='#1f77b4', edgecolor='black')
    ax2.axvline(trials_df['value'].mean(), color='red', linestyle='--', lw=2, label=f'Mean = {trials_df["value"].mean():.3f}')
    ax2.set_xlabel('Objective Value (ROC-AUC)', fontweight='bold')
    ax2.set_ylabel('Frequency', fontweight='bold')
    ax2.set_title('Distribution of Trial Values', fontweight='bold')
    ax2.grid(True, alpha=0.25, axis='y', linestyle=':')
    ax2.legend(loc='upper left')

    plt.tight_layout()
    filepath = os.path.join(output_dir, filename)
    plt.savefig(filepath, dpi=600, bbox_inches='tight', format='png')
    plt.close()

    print(f"✅ {filename} saved ({os.path.getsize(filepath)/1024:.1f} KB)")
    return filepath


def generate_all_figures(study, effect_sizes_dict: Dict, targets: List[str],
                         y_test, y_pred_vqc, y_pred_baseline,
                         output_dir: str = "04_figuras_publicacao"):
    """Gera todas as figuras de uma vez."""
    print(f"\n📊 GERANDO FIGURAS PUBLICAÇÃO ({output_dir}/)...")

    files_generated = []

    # Fig 1: Power curve
    fig1 = fig_power_curve(output_dir=output_dir)
    files_generated.append(fig1)

    # Fig 2: ROC comparison
    try:
        fig2 = fig_roc_comparison([], [], y_test, y_pred_vqc, y_pred_baseline, output_dir=output_dir)
        files_generated.append(fig2)
    except:
        print(f"⚠️  Pulando Fig 2 (dados insuficientes)")

    # Fig 3: Forest plot
    try:
        fig3 = fig_forest_plot(effect_sizes_dict, targets, output_dir=output_dir)
        files_generated.append(fig3)
    except:
        print(f"⚠️  Pulando Fig 3 (dados insuficientes)")

    # Fig 4: Optuna history
    try:
        fig4 = fig_optuna_history(study, output_dir=output_dir)
        files_generated.append(fig4)
    except:
        print(f"⚠️  Pulando Fig 4 (study não disponível)")

    print(f"\n✅ {len(files_generated)} figuras geradas com sucesso!")
    return files_generated


if __name__ == "__main__":
    # Teste
    print("Gerando figuras de exemplo...")
    fig_power_curve()
    print("✅ Exemplo concluído")

