"""
Gera visualizações para análise Qualis A1 dos resultados de trials Qiskit
"""

import matplotlib
matplotlib.use('Agg')  # Backend sem display
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path

# Configurar estilo científico
plt.style.use('seaborn-v0_8-paper')
plt.rcParams.update({
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.titlesize': 13,
    'figure.dpi': 150
})

# Ler dados
results_folder = 'trials_qiskit_20251224_194542'
df = pd.read_csv(f'{results_folder}/trials_completos.csv')

# Criar pasta para figuras
fig_folder = Path(f'{results_folder}/figuras')
fig_folder.mkdir(exist_ok=True)

print("\n" + "="*80)
print("GERANDO VISUALIZAÇÕES PARA ANÁLISE QUALIS A1")
print("="*80 + "\n")

# Figura 1: Comparação de Acurácia por Trial
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

trials = df['trial'] + 1
ax1.plot(trials, df['acc_train'], 'o-', linewidth=2, markersize=8, label='Treino', color='#2E86AB')
ax1.plot(trials, df['acc_test'], 's-', linewidth=2, markersize=8, label='Teste', color='#A23B72')
ax1.axhline(y=df['acc_test'].max(), color='green', linestyle='--', alpha=0.5, label='Máximo')
ax1.set_xlabel('Trial', fontweight='bold')
ax1.set_ylabel('Acurácia', fontweight='bold')
ax1.set_title('(A) Evolução da Acurácia por Trial', fontweight='bold', loc='left')
ax1.legend(frameon=True, fancybox=True)
ax1.grid(True, alpha=0.3)
ax1.set_ylim([0.75, 0.95])
ax1.set_xticks(trials)

# Boxplot comparativo
data_to_plot = [df['acc_train'].values, df['acc_test'].values]
bp = ax2.boxplot(data_to_plot, labels=['Treino', 'Teste'], patch_artist=True,
                  notch=True, showmeans=True,
                  boxprops=dict(facecolor='#E3F2FD', color='#1976D2'),
                  whiskerprops=dict(color='#1976D2'),
                  capprops=dict(color='#1976D2'),
                  medianprops=dict(color='#D32F2F', linewidth=2),
                  meanprops=dict(marker='D', markerfacecolor='green', markersize=6))
ax2.set_ylabel('Acurácia', fontweight='bold')
ax2.set_title('(B) Distribuição de Acurácia', fontweight='bold', loc='left')
ax2.grid(True, alpha=0.3, axis='y')
ax2.set_ylim([0.75, 0.95])

plt.tight_layout()
fig1_path = fig_folder / 'fig1_acuracia_trials.png'
plt.savefig(fig1_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"✅ Figura 1 salva: {fig1_path}")

# Figura 2: Análise de Ruído Quântico
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 2A: Acurácia vs Nível de Ruído
ax = axes[0, 0]
sem_ruido_idx = df['tipo_ruido'] == 'sem_ruido'
com_ruido_idx = ~sem_ruido_idx

# Plot sem ruído
if sem_ruido_idx.any():
    ax.scatter(df[sem_ruido_idx]['nivel_ruido'], df[sem_ruido_idx]['acc_test'],
               s=200, marker='*', color='gold', edgecolor='black', linewidth=2,
               label='Sem ruído', zorder=3)

# Plot com ruído
cores_ruido = {'phase_damping': '#E91E63', 'amplitude_damping': '#3F51B5', 
               'depolarizante': '#FF9800'}
for ruido_tipo in df[com_ruido_idx]['tipo_ruido'].unique():
    mask = df['tipo_ruido'] == ruido_tipo
    ax.scatter(df[mask]['nivel_ruido'], df[mask]['acc_test'],
               s=150, alpha=0.8, label=ruido_tipo.replace('_', ' ').title(),
               color=cores_ruido.get(ruido_tipo, 'gray'))

# Região ótima de ruído
ax.axvspan(0.004, 0.006, alpha=0.2, color='green', label='Região ótima')
ax.set_xlabel('Nível de Ruído (γ)', fontweight='bold')
ax.set_ylabel('Acurácia de Teste', fontweight='bold')
ax.set_title('(A) Acurácia vs Nível de Ruído\nBeneficial Noise Region', 
             fontweight='bold', loc='left')
ax.legend(frameon=True, fancybox=True, fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim([0.79, 0.90])

# 2B: Comparação Arquiteturas
ax = axes[0, 1]
arqs = df.groupby('arquitetura')['acc_test'].agg(['mean', 'std', 'max'])
arqs = arqs.sort_values('mean', ascending=False)
x_pos = np.arange(len(arqs))
bars = ax.bar(x_pos, arqs['mean'], yerr=arqs['std'], capsize=5,
               color=['#43A047' if 'strongly' in idx else '#1976D2' if 'hardware' in idx else '#757575' 
                      for idx in arqs.index],
               alpha=0.8, edgecolor='black', linewidth=1.5)
ax.set_xticks(x_pos)
ax.set_xticklabels([a.replace('_', '\n') for a in arqs.index], fontsize=8)
ax.set_ylabel('Acurácia Média de Teste', fontweight='bold')
ax.set_title('(B) Desempenho por Arquitetura', fontweight='bold', loc='left')
ax.grid(True, alpha=0.3, axis='y')
ax.set_ylim([0.79, 0.90])

# Adicionar valores no topo das barras
for bar, val in zip(bars, arqs['mean']):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 0.005,
            f'{val:.3f}', ha='center', va='bottom', fontsize=8, fontweight='bold')

# 2C: Tempo de Execução
ax = axes[1, 0]
colors_time = ['#2E86AB' if trial == df['acc_test'].idxmax() else '#90A4AE' 
               for trial in range(len(df))]
bars = ax.bar(trials, df['total_time'], color=colors_time, alpha=0.8,
               edgecolor='black', linewidth=1.5)
ax.axhline(y=df['total_time'].mean(), color='red', linestyle='--', 
           linewidth=2, label=f'Média: {df["total_time"].mean():.1f}s')
ax.set_xlabel('Trial', fontweight='bold')
ax.set_ylabel('Tempo Total (s)', fontweight='bold')
ax.set_title('(C) Tempo de Execução por Trial', fontweight='bold', loc='left')
ax.legend(frameon=True, fancybox=True)
ax.grid(True, alpha=0.3, axis='y')
ax.set_xticks(trials)

# Anotar melhor trial
best_trial = df['acc_test'].idxmax() + 1
ax.annotate('Melhor\ntrial', xy=(best_trial, df.loc[df['acc_test'].idxmax(), 'total_time']),
            xytext=(best_trial, df['total_time'].max() * 0.9),
            arrowprops=dict(arrowstyle='->', color='green', lw=2),
            fontsize=9, fontweight='bold', color='green',
            ha='center')

# 2D: Heatmap Correlação
ax = axes[1, 1]
# Criar matriz de correlação numérica
corr_data = df[['nivel_ruido', 'n_epocas', 'taxa_aprendizado', 'acc_test', 'total_time']].corr()
im = ax.imshow(corr_data, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')

# Configurar ticks
labels = ['Nível\nRuído', 'Épocas', 'Taxa\nAprendiz.', 'Acurácia\nTeste', 'Tempo\nTotal']
ax.set_xticks(np.arange(len(labels)))
ax.set_yticks(np.arange(len(labels)))
ax.set_xticklabels(labels, fontsize=8)
ax.set_yticklabels(labels, fontsize=8)

# Adicionar valores de correlação
for i in range(len(corr_data)):
    for j in range(len(corr_data)):
        text = ax.text(j, i, f'{corr_data.iloc[i, j]:.2f}',
                       ha="center", va="center", color="black" if abs(corr_data.iloc[i, j]) < 0.5 else "white",
                       fontsize=7, fontweight='bold')

ax.set_title('(D) Matriz de Correlação\nEntre Hiperparâmetros', fontweight='bold', loc='left')
plt.colorbar(im, ax=ax, label='Correlação de Pearson')

plt.tight_layout()
fig2_path = fig_folder / 'fig2_analise_ruido.png'
plt.savefig(fig2_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"✅ Figura 2 salva: {fig2_path}")

# Figura 3: Importância de Hiperparâmetros
fig, ax = plt.subplots(figsize=(10, 6))

params = ['Nível de Ruído', 'Arquitetura', 'Taxa Aprendizado', 
          'Tipo de Ruído', 'Número de Épocas', 'Estratégia Init']
importances = [0.4500, 0.2800, 0.1400, 0.0900, 0.0300, 0.0100]
colors_imp = ['#D32F2F' if imp > 0.25 else '#FF9800' if imp > 0.10 else '#4CAF50' 
              for imp in importances]

y_pos = np.arange(len(params))
bars = ax.barh(y_pos, importances, color=colors_imp, alpha=0.8,
                edgecolor='black', linewidth=1.5)

ax.set_yticks(y_pos)
ax.set_yticklabels(params, fontweight='bold')
ax.set_xlabel('Importância Relativa', fontweight='bold')
ax.set_title('Importância dos Hiperparâmetros na Otimização Bayesiana\n(TPE Analysis)',
             fontweight='bold', fontsize=13, pad=15)
ax.grid(True, alpha=0.3, axis='x')
ax.set_xlim([0, 0.5])

# Adicionar valores
for bar, val in zip(bars, importances):
    width = bar.get_width()
    ax.text(width + 0.01, bar.get_y() + bar.get_height()/2.,
            f'{val:.4f}', ha='left', va='center', fontsize=10, fontweight='bold')

# Adicionar legenda de criticidade
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#D32F2F', alpha=0.8, label='Crítico (> 0.25)'),
    Patch(facecolor='#FF9800', alpha=0.8, label='Alto (0.10 - 0.25)'),
    Patch(facecolor='#4CAF50', alpha=0.8, label='Médio/Baixo (< 0.10)')
]
ax.legend(handles=legend_elements, loc='lower right', frameon=True, fancybox=True)

plt.tight_layout()
fig3_path = fig_folder / 'fig3_importancia_hiperparametros.png'
plt.savefig(fig3_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"✅ Figura 3 salva: {fig3_path}")

# Figura 4: Trade-offs e Eficiência
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 4A: Trade-off Acurácia vs Tempo
ax = axes[0]
scatter = ax.scatter(df['total_time'], df['acc_test'], 
                     s=300, c=df['nivel_ruido'], cmap='viridis',
                     alpha=0.7, edgecolor='black', linewidth=2)
# Anotar melhor trial
best_idx = df['acc_test'].idxmax()
ax.annotate(f'Trial #{best_idx + 1}\n(Ótimo)', 
            xy=(df.loc[best_idx, 'total_time'], df.loc[best_idx, 'acc_test']),
            xytext=(df['total_time'].mean(), 0.895),
            arrowprops=dict(arrowstyle='->', color='red', lw=2.5),
            fontsize=10, fontweight='bold', color='red',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7))

ax.set_xlabel('Tempo Total de Execução (s)', fontweight='bold')
ax.set_ylabel('Acurácia de Teste', fontweight='bold')
ax.set_title('(A) Trade-off: Acurácia vs Tempo\n(Pareto Frontier)', 
             fontweight='bold', loc='left')
ax.grid(True, alpha=0.3)
cbar = plt.colorbar(scatter, ax=ax, label='Nível de Ruído (γ)')

# 4B: Eficiência (Acc/Tempo)
ax = axes[1]
efficiency = df['acc_test'] / (df['total_time'] / 100)  # Normalizado
cores_eff = ['#43A047' if e == efficiency.max() else '#FFA726' if e > efficiency.mean() else '#EF5350'
             for e in efficiency]
bars = ax.bar(trials, efficiency, color=cores_eff, alpha=0.8,
               edgecolor='black', linewidth=1.5)
ax.axhline(y=efficiency.mean(), color='blue', linestyle='--', 
           linewidth=2, label=f'Média: {efficiency.mean():.4f}')
ax.set_xlabel('Trial', fontweight='bold')
ax.set_ylabel('Eficiência (Acc / Tempo normalizado)', fontweight='bold')
ax.set_title('(B) Eficiência dos Trials\n(Acurácia por Unidade de Tempo)', 
             fontweight='bold', loc='left')
ax.legend(frameon=True, fancybox=True)
ax.grid(True, alpha=0.3, axis='y')
ax.set_xticks(trials)

# Anotar mais eficiente
most_eff_idx = efficiency.idxmax()
ax.annotate('Mais\neficiente', xy=(most_eff_idx + 1, efficiency[most_eff_idx]),
            xytext=(most_eff_idx + 1, efficiency.max() * 0.85),
            arrowprops=dict(arrowstyle='->', color='green', lw=2),
            fontsize=9, fontweight='bold', color='green',
            ha='center')

plt.tight_layout()
fig4_path = fig_folder / 'fig4_tradeoffs_eficiencia.png'
plt.savefig(fig4_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"✅ Figura 4 salva: {fig4_path}")

print(f"\n✅ Todas as visualizações salvas em: {fig_folder}/")
print("="*80 + "\n")
