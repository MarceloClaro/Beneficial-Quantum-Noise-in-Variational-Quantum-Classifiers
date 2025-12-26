#!/usr/bin/env python3
"""
Script para gerar visualizações dos resultados QAOA
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# Diretório de resultados
results_dir = sys.argv[1] if len(sys.argv) > 1 else "qaoa_framework/qaoa_framework/results/qaoa_run_20251226_231733"

# Carregar resultados
df = pd.read_csv(f"{results_dir}/resultados.csv")
resumo = pd.read_csv(f"{results_dir}/resumo.csv")

print("="*80)
print("RESULTADOS DO EXPERIMENTO QAOA")
print("="*80)
print(f"\nDiretório: {results_dir}\n")

print("RESUMO ESTATÍSTICO:")
print(resumo.to_string(index=False))
print()

print("\nRESULTADOS INDIVIDUAIS:")
print(df.to_string(index=False))
print()

# Análise de ruído benéfico
print("\n" + "="*80)
print("ANÁLISE DE RUÍDO BENÉFICO")
print("="*80)

sem_ruido = df[df['tipo_ruido'] == 'sem_ruido']['approx_ratio']
com_ruido = df[df['tipo_ruido'] == 'depolarizing']['approx_ratio']

print(f"\nApproximation Ratio SEM ruído:")
print(f"  Média: {sem_ruido.mean():.4f} ± {sem_ruido.std():.4f}")
print(f"  Min: {sem_ruido.min():.4f}, Max: {sem_ruido.max():.4f}")

print(f"\nApproximation Ratio COM ruído (depolarizing, p=0.005):")
print(f"  Média: {com_ruido.mean():.4f} ± {com_ruido.std():.4f}")
print(f"  Min: {com_ruido.min():.4f}, Max: {com_ruido.max():.4f}")

diferenca = ((com_ruido.mean() - sem_ruido.mean()) / sem_ruido.mean()) * 100
if diferenca > 0:
    print(f"\n✅ RUÍDO BENÉFICO DETECTADO!")
    print(f"   Melhoria de {diferenca:+.2f}% com ruído")
else:
    print(f"\n❌ Ruído prejudicial detectado")
    print(f"   Degradação de {diferenca:.2f}%")

# Criar visualizações
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Análise de Ruído Benéfico em QAOA', fontsize=16, fontweight='bold')

# 1. Comparação de Approximation Ratio
ax1 = axes[0, 0]
df_plot = df.groupby('tipo_ruido')['approx_ratio'].agg(['mean', 'std']).reset_index()
x = range(len(df_plot))
ax1.bar(x, df_plot['mean'], yerr=df_plot['std'], capsize=5, 
        color=['#3498db', '#e74c3c'], alpha=0.7)
ax1.set_xticks(x)
ax1.set_xticklabels(df_plot['tipo_ruido'], rotation=0)
ax1.set_ylabel('Approximation Ratio')
ax1.set_title('Comparação: Sem Ruído vs. Com Ruído')
ax1.grid(axis='y', alpha=0.3)
for i, (mean, std) in enumerate(zip(df_plot['mean'], df_plot['std'])):
    ax1.text(i, mean + std + 0.01, f'{mean:.3f}', ha='center', fontweight='bold')

# 2. Energia Final
ax2 = axes[0, 1]
df_plot_e = df.groupby('tipo_ruido')['energia_final'].agg(['mean', 'std']).reset_index()
ax2.bar(x, df_plot_e['mean'], yerr=df_plot_e['std'], capsize=5,
        color=['#3498db', '#e74c3c'], alpha=0.7)
ax2.set_xticks(x)
ax2.set_xticklabels(df_plot_e['tipo_ruido'], rotation=0)
ax2.set_ylabel('Energia Final')
ax2.set_title('Energia do Hamiltoniano')
ax2.grid(axis='y', alpha=0.3)

# 3. Distribuição de Approximation Ratio
ax3 = axes[1, 0]
sem_ruido_vals = df[df['tipo_ruido'] == 'sem_ruido']['approx_ratio']
com_ruido_vals = df[df['tipo_ruido'] == 'depolarizing']['approx_ratio']
ax3.boxplot([sem_ruido_vals, com_ruido_vals], labels=['Sem Ruído', 'Com Ruído'])
ax3.set_ylabel('Approximation Ratio')
ax3.set_title('Distribuição dos Resultados')
ax3.grid(axis='y', alpha=0.3)

# 4. Tempo de Execução
ax4 = axes[1, 1]
df_plot_t = df.groupby('tipo_ruido')['tempo_execucao'].agg(['mean', 'std']).reset_index()
ax4.bar(x, df_plot_t['mean'], yerr=df_plot_t['std'], capsize=5,
        color=['#2ecc71', '#95a5a6'], alpha=0.7)
ax4.set_xticks(x)
ax4.set_xticklabels(df_plot_t['tipo_ruido'], rotation=0)
ax4.set_ylabel('Tempo (s)')
ax4.set_title('Tempo de Execução')
ax4.grid(axis='y', alpha=0.3)

plt.tight_layout()
output_file = f"{results_dir}/visualizacao.png"
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"\n📊 Visualização salva em: {output_file}")
print("="*80)
