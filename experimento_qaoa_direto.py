#!/usr/bin/env python3
"""
EXPERIMENTO QAOA 100 QUBITS - ANÁLISE DE RUÍDO BENÉFICO
Versão simplificada para execução direta
"""

import os
import sys
import time
import json
import logging
from pathlib import Path
from datetime import datetime

print("\n" + "="*80)
print("🚀 INICIANDO EXPERIMENTO QAOA 100 QUBITS COM RUÍDO BENÉFICO")
print("="*80)

# Step 1: Verificar dependências
print("\n[1/5] Verificando dependências...")
try:
    import numpy as np
    import pandas as pd
    from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, depolarizing_error, phase_damping_error
    print("✅ Todas as dependências importadas com sucesso")
except ImportError as e:
    print(f"❌ Erro de importação: {e}")
    print("⚠️ Tente instalar: pip install qiskit qiskit-aer numpy pandas")
    sys.exit(1)

# Step 2: Importar framework
print("\n[2/5] Importando framework QAOA...")
try:
    from framework_qaoa_100qubits import (
        ConfigQAOA,
        ConstrutorCircuitoQAOA,
        OtimizadorQAOA,
        AnalisadorHiperparametrosQAOA
    )
    print("✅ Framework QAOA importado com sucesso")
except ImportError as e:
    print(f"❌ Erro ao importar framework: {e}")
    sys.exit(1)

# Step 3: Criar diretório de resultados
print("\n[3/5] Preparando armazenamento de resultados...")
pasta_resultados = Path('resultados_qaoa_experimento_completo')
pasta_resultados.mkdir(exist_ok=True, parents=True)
print(f"✅ Pasta criada: {pasta_resultados.absolute()}")

# Step 4: Configurar experimento
print("\n[4/5] Configurando experimento...")

config = ConfigQAOA(
    n_qubits=100,
    p_layers=3,
    tipo_ruido='depolarizing',
    nivel_ruido=0.001,
    shots=1024,
    max_iter=100,
    seed=42,
    problema='maxcut',
    otimizador='COBYLA'
)

print(f"✅ Configuração pronta:")
print(f"   • Qubits: {config.n_qubits}")
print(f"   • P-layers: {config.p_layers}")
print(f"   • Ruído: {config.tipo_ruido} (nível {config.nivel_ruido})")
print(f"   • Shots: {config.shots}")

# Step 5: Executar experimentos
print("\n[5/5] Executando experimentos...")
print("-" * 80)

# Criar estrutura de experimentos
experimentos = []
inicio = time.time()

# Experimento 1: Sem ruído
print("\n📊 Experimento 1: SEM RUÍDO")
print("Executando...")
try:
    construtor = ConstrutorCircuitoQAOA(n_qubits=100, p_layers=3, seed=42)
    grafo = construtor.criar_grafo_aleatorio(densidade=0.1)
    
    config_sem_ruido = ConfigQAOA(
        n_qubits=100, p_layers=3, tipo_ruido='sem_ruido',
        nivel_ruido=0.0, max_iter=50, shots=512, seed=42
    )
    
    otimizador = OtimizadorQAOA(config_sem_ruido)
    resultado_sem_ruido = otimizador.otimizar(grafo)
    
    experimentos.append({
        'nome': 'Sem Ruído',
        'config': config_sem_ruido.__dict__,
        'energia': resultado_sem_ruido.energia_final if hasattr(resultado_sem_ruido, 'energia_final') else 'N/A',
        'tempo': resultado_sem_ruido.tempo_execucao if hasattr(resultado_sem_ruido, 'tempo_execucao') else 'N/A'
    })
    
    print(f"   ✅ Energia final: {experimentos[-1]['energia']}")
    print(f"   ✅ Tempo: {experimentos[-1]['tempo']}s")
    
except Exception as e:
    print(f"   ⚠️ Erro neste experimento: {e}")
    experimentos.append({
        'nome': 'Sem Ruído',
        'status': 'Erro',
        'erro': str(e)
    })

# Experimento 2: Com ruído depolarizing
print("\n📊 Experimento 2: COM RUÍDO DEPOLARIZING (0.001)")
print("Executando...")
try:
    config_com_ruido = ConfigQAOA(
        n_qubits=100, p_layers=3, tipo_ruido='depolarizing',
        nivel_ruido=0.001, max_iter=50, shots=512, seed=42
    )
    
    otimizador = OtimizadorQAOA(config_com_ruido)
    resultado_com_ruido = otimizador.otimizar(grafo)
    
    experimentos.append({
        'nome': 'Com Ruído Depolarizing',
        'config': config_com_ruido.__dict__,
        'energia': resultado_com_ruido.energia_final if hasattr(resultado_com_ruido, 'energia_final') else 'N/A',
        'tempo': resultado_com_ruido.tempo_execucao if hasattr(resultado_com_ruido, 'tempo_execucao') else 'N/A'
    })
    
    print(f"   ✅ Energia final: {experimentos[-1]['energia']}")
    print(f"   ✅ Tempo: {experimentos[-1]['tempo']}s")
    
except Exception as e:
    print(f"   ⚠️ Erro neste experimento: {e}")
    experimentos.append({
        'nome': 'Com Ruído Depolarizing',
        'status': 'Erro',
        'erro': str(e)
    })

# Experimento 3: Com ruído phase damping
print("\n📊 Experimento 3: COM RUÍDO PHASE DAMPING (0.001)")
print("Executando...")
try:
    config_phase = ConfigQAOA(
        n_qubits=100, p_layers=3, tipo_ruido='phase_damping',
        nivel_ruido=0.001, max_iter=50, shots=512, seed=42
    )
    
    otimizador = OtimizadorQAOA(config_phase)
    resultado_phase = otimizador.otimizar(grafo)
    
    experimentos.append({
        'nome': 'Com Ruído Phase Damping',
        'config': config_phase.__dict__,
        'energia': resultado_phase.energia_final if hasattr(resultado_phase, 'energia_final') else 'N/A',
        'tempo': resultado_phase.tempo_execucao if hasattr(resultado_phase, 'tempo_execucao') else 'N/A'
    })
    
    print(f"   ✅ Energia final: {experimentos[-1]['energia']}")
    print(f"   ✅ Tempo: {experimentos[-1]['tempo']}s")
    
except Exception as e:
    print(f"   ⚠️ Erro neste experimento: {e}")
    experimentos.append({
        'nome': 'Com Ruído Phase Damping',
        'status': 'Erro',
        'erro': str(e)
    })

# Salvar resultados
tempo_total = time.time() - inicio

print("\n" + "="*80)
print("📊 RESUMO DOS RESULTADOS")
print("="*80)

# Criar DataFrame com resultados
df_resultados = pd.DataFrame(experimentos)

# Salvar em CSV
arquivo_csv = pasta_resultados / f"resultados_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
df_resultados.to_csv(arquivo_csv, index=False)
print(f"\n✅ Resultados salvos em: {arquivo_csv}")

# Salvar em JSON
arquivo_json = pasta_resultados / f"resumo_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
with open(arquivo_json, 'w') as f:
    json.dump({
        'timestamp': datetime.now().isoformat(),
        'tempo_total_segundos': tempo_total,
        'num_experimentos': len(experimentos),
        'experimentos': experimentos
    }, f, indent=2, default=str)
print(f"✅ Resumo salvo em: {arquivo_json}")

# Exibir tabela
print("\n" + "-"*80)
print("Tabela de Resultados:")
print("-"*80)
print(df_resultados.to_string(index=False))

print("\n" + "="*80)
print(f"⏱️ TEMPO TOTAL DE EXECUÇÃO: {tempo_total/60:.2f} minutos ({tempo_total:.1f} segundos)")
print("="*80)

print("\n✅ EXPERIMENTO CONCLUÍDO COM SUCESSO!")
print(f"📁 Arquivos salvos em: {pasta_resultados.absolute()}")
print("\n" + "="*80 + "\n")

