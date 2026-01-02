#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXPERIMENTO SIMPLES - Framework V8 - Modo Rapido
"""

import os
import sys
import time
import json
from pathlib import Path
from datetime import datetime

WORKSPACE = Path(__file__).parent
sys.path.insert(0, str(WORKSPACE))
RESULTS_DIR = WORKSPACE / "experimento_completo_resultados"
RESULTS_DIR.mkdir(exist_ok=True)

print("=" * 80)
print("FRAMEWORK V8 - TESTE SIMPLIFICADO E RAPIDO")
print("=" * 80)
print(f"Data: {datetime.now().strftime('%d de %B de %Y as %H:%M:%S')}\n")

# ============================================================================
# TESTE 1: VERIFICAR IMPORTACOES
# ============================================================================
print("[1/3] Verificando importacoes...")

try:
    import pennylane as qml
    from sklearn.datasets import load_breast_cancer
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import train_test_split
    from framework_investigativo_completo import ClassificadorVQC
    print("     OK - Todas as bibliotecas importadas\n")
    imports_ok = True
except Exception as e:
    print(f"     ERROR - {str(e)[:60]}\n")
    imports_ok = False

# ============================================================================
# TESTE 2: VERIFICAR DATASET
# ============================================================================
print("[2/3] Verificando dados...")

try:
    cancer = load_breast_cancer()
    X = StandardScaler().fit_transform(cancer.data)
    y = cancer.target
    
    X_train, X_test, y_train, y_test = train_test_split(
        X[:100], y[:100], test_size=0.3, random_state=42
    )
    
    print(f"     OK - Dataset carregado ({X_train.shape[0]} treino, {X_test.shape[0]} teste)\n")
    dataset_ok = True
except Exception as e:
    print(f"     ERROR - {str(e)[:60]}\n")
    dataset_ok = False

# ============================================================================
# TESTE 3: VERIFICAR CLASSIFIADOR
# ============================================================================
print("[3/3] Verificando ClassificadorVQC...\n")

results = []

if imports_ok and dataset_ok:
    print("     Testando com parametros MINIMOS para speed...\n")
    
    try:
        vqc = ClassificadorVQC(
            n_qubits=2,
            n_camadas=1,
            tipo_ruido='sem_ruido',
            n_epocas=1,
            early_stopping=False,
            batch_size=10,
            learning_rate=0.01
        )
        
        print("     [*] Iniciando fit (pode levar alguns segundos)...")
        inicio = time.time()
        
        vqc.fit(X_train[:10], y_train[:10])
        
        duracao = time.time() - inicio
        
        print(f"     [*] Fit completado em {duracao:.2f}s")
        
        acuracia = vqc.score(X_test[:5], y_test[:5])
        print(f"     [*] Acuracia: {acuracia:.4f}\n")
        
        results.append({
            'teste': 'ClassificadorVQC Basico',
            'status': 'PASS',
            'acuracia': float(acuracia),
            'tempo_segundos': float(duracao)
        })
        
    except KeyboardInterrupt:
        print("\n     TIMEOUT - Teste interrompido (executando muito lentamente)\n")
        results.append({
            'teste': 'ClassificadorVQC Basico',
            'status': 'TIMEOUT'
        })
    except Exception as e:
        print(f"     ERROR - {str(e)[:60]}\n")
        results.append({
            'teste': 'ClassificadorVQC Basico',
            'status': 'ERROR',
            'error': str(e)[:60]
        })

# ============================================================================
# SALVAR RESULTADOS
# ============================================================================
print("=" * 80)
print("RESULTADOS")
print("=" * 80 + "\n")

results_json = RESULTS_DIR / "resultados_simples.json"
with open(results_json, 'w') as f:
    json.dump(results, f, indent=2)

print(f"[OK] Resultados salvos: {results_json.name}\n")

# ============================================================================
# RELATORIO FINAL
# ============================================================================
print("=" * 80)
print("RELATORIO")
print("=" * 80 + "\n")

report = "# Framework V8 - Teste Simplificado\n\n"
report += f"Data: {datetime.now().strftime('%d de %B de %Y as %H:%M:%S')}\n\n"

report += "## Verificacoes\n\n"
report += f"- Importacoes: {'OK' if imports_ok else 'ERROR'}\n"
report += f"- Dataset: {'OK' if dataset_ok else 'ERROR'}\n"
report += f"- ClassificadorVQC: {results[0]['status'] if results else 'N/A'}\n\n"

report += "## Conclusao\n\n"
if results and results[0]['status'] == 'PASS':
    report += f"Framework V8 esta operacional!\n"
    report += f"Acuracia: {results[0]['acuracia']:.4f}\n"
    report += f"Tempo treino: {results[0]['tempo_segundos']:.2f}s\n"
elif results and results[0]['status'] == 'TIMEOUT':
    report += "Framework V8 precisa otimizacao (muito lento).\n"
    report += "Recomendacao: usar menos qubits/layers para producao.\n"
else:
    report += "Framework V8 tem problemas.\n"

report_file = RESULTS_DIR / "RELATORIO_SIMPLES.md"
with open(report_file, 'w') as f:
    f.write(report)

print(f"[OK] Relatorio criado: {report_file.name}\n")

# ============================================================================
# RESUMO
# ============================================================================
arquivos = list(RESULTS_DIR.glob("*"))
print("=" * 80)
print("ARQUIVOS GERADOS")
print("=" * 80)
print(f"\nTotal: {len(arquivos)} arquivos\n")

for arq in sorted(arquivos):
    print(f"  - {arq.name}")

print("\n" + "=" * 80)
print("TESTE FINALIZADO")
print("=" * 80 + "\n")
