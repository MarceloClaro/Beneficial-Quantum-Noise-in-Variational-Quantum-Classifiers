#!/usr/bin/env python3
"""
Tarefa 11: Script de Demonstração QAOA
Demonstra uso completo do framework QAOA com análise de ruído benéfico.

Uso:
    python demo_qaoa_rapido.py
    
Autor: Framework QAOA
Data: 2025-12-26
"""

import sys
from pathlib import Path

# Adicionar diretório de scripts ao path
sys.path.insert(0, str(Path(__file__).parent / 'scripts'))

import numpy as np
import pandas as pd

print("="*80)
print(" "*20 + "DEMONSTRAÇÃO QAOA 100 QUBITS")
print(" "*15 + "Análise de Ruído Quântico Benéfico")
print("="*80 + "\n")

# Verificar dependências
print("Verificando dependências...")
try:
    from scripts.problem_generator import gerar_problema_benchmark
    from scripts.circuit_builder import QAOACircuitBuilder
    from scripts.noise_models import criar_noise_model_qiskit, MODELOS_RUIDO
    print("✓ Módulos do framework carregados")
except ImportError as e:
    print(f"❌ Erro ao importar módulos: {e}")
    sys.exit(1)

try:
    from qiskit import transpile
    from qiskit_aer import AerSimulator
    print("✓ Qiskit disponível")
except ImportError:
    print("❌ Qiskit não disponível. Instale com: pip install qiskit qiskit-aer")
    sys.exit(1)

print("\n" + "-"*80)
print("DEMO 1: PROBLEMA MAXCUT PEQUENO (10 QUBITS)")
print("-"*80 + "\n")

# 1. Criar problema
print("[1/4] Gerando problema MaxCut...")
problem = gerar_problema_benchmark(
    graph_type='erdos_renyi',
    n_nodes=10,
    edge_probability=0.5,
    seed=42
)

print(f"  • Nós: {problem.n_nodes}")
print(f"  • Arestas: {len(problem.edges)}")
print(f"  • Valor ótimo: {problem.optimal_value:.2f}")

# 2. Construir circuito
print("\n[2/4] Construindo circuito QAOA...")
builder = QAOACircuitBuilder(
    n_qubits=problem.n_nodes,
    p_layers=2,
    hamiltonian=problem.to_hamiltonian()
)

qc = builder.build()
stats = builder.get_circuit_stats()

print(f"  • P-layers: 2")
print(f"  • Depth: {stats['depth']}")
print(f"  • Gates: {stats['total_gates']} (CX: {stats['cx_gates']})")
print(f"  • Parâmetros: {stats['parameters']}")

# 3. Executar sem ruído
print("\n[3/4] Executando QAOA SEM RUÍDO...")

from scipy.optimize import minimize

def evaluate_simple(params, circuit, problem, simulator, shots=512):
    """Avaliação simplificada."""
    bound_circuit = circuit.bind_parameters(
        {builder.params[i]: params[i] for i in range(len(params))}
    )
    
    transpiled = transpile(bound_circuit, simulator)
    job = simulator.run(transpiled, shots=shots)
    result = job.result()
    counts = result.get_counts()
    
    # Calcular energia
    energia = 0.0
    total = sum(counts.values())
    
    for bitstring, count in counts.items():
        cut = [int(b) for b in bitstring[::-1]]
        energia_config = 0.0
        
        for (i, j), coef in problem.to_hamiltonian().items():
            zi_zj = 1 if cut[i] == cut[j] else -1
            energia_config += coef * zi_zj
        
        energia += (count / total) * energia_config
    
    return energia

# Simulador sem ruído
simulator_clean = AerSimulator(method='automatic')
params_init = builder.initialize_parameters('heuristic', seed=42)

result_clean = minimize(
    lambda p: evaluate_simple(p, qc, problem, simulator_clean),
    params_init,
    method='COBYLA',
    options={'maxiter': 50}
)

energia_clean = result_clean.fun
approx_ratio_clean = abs(energia_clean / problem.optimal_value)

print(f"  • Energia final: {energia_clean:.4f}")
print(f"  • Approximation ratio: {approx_ratio_clean:.4f}")
print(f"  • Iterações: {result_clean.nit}")

# 4. Executar com ruído
print("\n[4/4] Executando QAOA COM RUÍDO (depolarizing, 0.001)...")

noise_model = criar_noise_model_qiskit('depolarizing', 0.001)
simulator_noisy = AerSimulator(noise_model=noise_model, method='automatic')

result_noisy = minimize(
    lambda p: evaluate_simple(p, qc, problem, simulator_noisy),
    params_init,
    method='COBYLA',
    options={'maxiter': 50}
)

energia_noisy = result_noisy.fun
approx_ratio_noisy = abs(energia_noisy / problem.optimal_value)

print(f"  • Energia final: {energia_noisy:.4f}")
print(f"  • Approximation ratio: {approx_ratio_noisy:.4f}")
print(f"  • Iterações: {result_noisy.nit}")

# Análise
print("\n" + "-"*80)
print("ANÁLISE COMPARATIVA")
print("-"*80 + "\n")

print(f"Sem ruído:   approx_ratio = {approx_ratio_clean:.4f}")
print(f"Com ruído:   approx_ratio = {approx_ratio_noisy:.4f}")

diferenca = approx_ratio_noisy - approx_ratio_clean
diferenca_pct = (diferenca / approx_ratio_clean) * 100

print(f"Diferença:   {diferenca:+.4f} ({diferenca_pct:+.2f}%)")

if diferenca > 0.01:
    print("\n✅ RUÍDO BENÉFICO DETECTADO!")
    print(f"   O ruído melhorou o approximation ratio em {abs(diferenca_pct):.2f}%")
elif diferenca < -0.01:
    print("\n⚠️  RUÍDO PREJUDICIAL neste caso")
    print(f"   O ruído piorou o approximation ratio em {abs(diferenca_pct):.2f}%")
else:
    print("\n➖ IMPACTO NEUTRO")
    print("   O ruído teve impacto mínimo (< 1%)")

print("\n" + "-"*80)
print("DEMO 2: TESTE DE ESCALABILIDADE (100 QUBITS - ESTRUTURA)")
print("-"*80 + "\n")

print("[1/2] Criando problema com 100 qubits...")
problem_large = gerar_problema_benchmark(
    graph_type='erdos_renyi',
    n_nodes=100,
    edge_probability=0.1,
    seed=42
)

print(f"  • Nós: {problem_large.n_nodes}")
print(f"  • Arestas: {len(problem_large.edges)}")
print(f"  • Valor aproximado: {problem_large.optimal_value:.2f}")

print("\n[2/2] Construindo circuito QAOA para 100 qubits...")
builder_large = QAOACircuitBuilder(
    n_qubits=100,
    p_layers=3,
    hamiltonian=problem_large.to_hamiltonian()
)

qc_large = builder_large.build(measurements=False)  # Sem medições para análise
stats_large = builder_large.get_circuit_stats()

print(f"  • P-layers: 3")
print(f"  • Depth: {stats_large['depth']}")
print(f"  • Gates: {stats_large['total_gates']} (CX: {stats_large['cx_gates']})")
print(f"  • Parâmetros: {stats_large['parameters']}")

print("\n⚠️  Nota: Execução completa com 100 qubits requer:")
print("   - MPS backend (method='matrix_product_state')")
print("   - Memória: ~16GB RAM recomendado")
print("   - Tempo: ~1-2 horas para otimização completa")

print("\n" + "="*80)
print("DEMONSTRAÇÃO CONCLUÍDA")
print("="*80 + "\n")

print("📚 Próximos Passos:")
print("  1. Execute main.py para experimento completo com YAML config")
print("  2. Use hyperparameter_tuning.py para otimização Bayesiana")
print("  3. Use visualization.py para gerar gráficos dos resultados")
print()
print("📖 Documentação: qaoa_framework/README.md")
print("⚙️  Configuração: qaoa_framework/configs/experiment_qaoa.yaml")
print()
