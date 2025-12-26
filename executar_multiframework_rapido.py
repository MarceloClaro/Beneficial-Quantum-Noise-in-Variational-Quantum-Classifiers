"""
Script rápido para executar o multiframework (PennyLane, Qiskit, Cirq) e gerar resultados comparativos.
Versão otimizada com menos épocas para execução mais rápida.
"""

import os
import sys
import json
import time
import logging
from datetime import datetime
from pathlib import Path
import pandas as pd
import numpy as np

# Configurar logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Carregar configuração
with open('qai_config.json', 'r') as f:
    config = json.load(f)

logger.info("="*80)
logger.info("EXECUÇÃO MULTIFRAMEWORK RÁPIDA - QUALIS A1 ENHANCED")
logger.info(f"Versão: {config['version']}")
logger.info(f"Frameworks habilitados: {', '.join(config['enabled_frameworks'])}")
logger.info("="*80)

# Criar pasta de resultados
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
pasta_resultados = f"resultados_multiframework_{timestamp}"
os.makedirs(pasta_resultados, exist_ok=True)
logger.info(f"\n📁 Pasta de resultados: {pasta_resultados}/")

# Estrutura para armazenar todos os resultados
resultados_multiframework = {
    'timestamp': timestamp,
    'config': config,
    'frameworks': {},
    'comparacao': {}
}

# ============================================================================
# 1. FRAMEWORK PENNYLANE (Baseline)
# ============================================================================
if 'pennylane' in config['enabled_frameworks']:
    logger.info("\n" + "="*80)
    logger.info("1. EXECUTANDO FRAMEWORK PENNYLANE (BASELINE)")
    logger.info("="*80)
    
    try:
        from framework_investigativo_completo import (
            ClassificadorVQC,
            carregar_datasets
        )
        
        datasets = carregar_datasets()
        dataset_moons = datasets['moons']
        
        # Configuração focada para teste muito rápido
        config_pnl = {
            'arquitetura': 'strongly_entangling',
            'tipo_ruido': 'phase_damping',
            'nivel_ruido': 0.005,
            'estrategia_init': 'quantico'
        }
        
        logger.info(f"\n  Configuração: {config_pnl}")
        
        inicio = time.time()
        
        vqc = ClassificadorVQC(
            n_qubits=4,
            n_camadas=2,
            arquitetura=config_pnl['arquitetura'],
            estrategia_init=config_pnl['estrategia_init'],
            seed=config['default_seed']
        )
        
        vqc.tipo_ruido = config_pnl['tipo_ruido']
        vqc.nivel_ruido = config_pnl['nivel_ruido']
        
        vqc.fit(
            dataset_moons['X_train'][:30],  # Subset reduzido
            dataset_moons['y_train'][:30],
            n_epocas=5  # Muito reduzido
        )
        
        acc_test = vqc.score(dataset_moons['X_test'][:15], dataset_moons['y_test'][:15])
        tempo = time.time() - inicio
        
        resultado = {
            'framework': 'pennylane',
            'arquitetura': config_pnl['arquitetura'],
            'tipo_ruido': config_pnl['tipo_ruido'],
            'nivel_ruido': config_pnl['nivel_ruido'],
            'estrategia_init': config_pnl['estrategia_init'],
            'acuracia_teste': float(acc_test),
            'tempo_segundos': float(tempo)
        }
        
        resultados_multiframework['frameworks']['pennylane'] = {
            'status': 'success',
            'resultados': [resultado],
            'melhor_acuracia': float(acc_test)
        }
        
        logger.info(f"  ✅ PennyLane: Acurácia={acc_test:.4f}, Tempo={tempo:.2f}s")
        
    except Exception as e:
        logger.error(f"  ❌ Erro no PennyLane: {e}")
        resultados_multiframework['frameworks']['pennylane'] = {
            'status': 'error',
            'erro': str(e)
        }

# ============================================================================
# 2. FRAMEWORK QISKIT
# ============================================================================
if 'qiskit' in config['enabled_frameworks']:
    logger.info("\n" + "="*80)
    logger.info("2. EXECUTANDO FRAMEWORK QISKIT")
    logger.info("="*80)
    
    try:
        from framework_qiskit import (
            ClassificadorVQCQiskit,
            carregar_datasets as carregar_datasets_qiskit
        )
        
        datasets_qiskit = carregar_datasets_qiskit()
        dataset_moons_q = datasets_qiskit['moons']
        
        config_qiskit = {
            'arquitetura': 'strongly_entangling',
            'tipo_ruido': 'phase_damping',
            'nivel_ruido': 0.005
        }
        
        logger.info(f"\n  Configuração: {config_qiskit}")
        
        inicio = time.time()
        
        vqc = ClassificadorVQCQiskit(
            n_qubits=4,
            n_camadas=2,
            arquitetura=config_qiskit['arquitetura'],
            tipo_ruido=config_qiskit['tipo_ruido'],
            nivel_ruido=config_qiskit['nivel_ruido'],
            n_epocas=5,
            seed=config['default_seed']
        )
        
        vqc.fit(dataset_moons_q['X_train'][:30], dataset_moons_q['y_train'][:30])
        
        acc_test = vqc.score(dataset_moons_q['X_test'][:15], dataset_moons_q['y_test'][:15])
        tempo = time.time() - inicio
        
        resultado = {
            'framework': 'qiskit',
            'arquitetura': config_qiskit['arquitetura'],
            'tipo_ruido': config_qiskit['tipo_ruido'],
            'nivel_ruido': config_qiskit['nivel_ruido'],
            'acuracia_teste': float(acc_test),
            'tempo_segundos': float(tempo)
        }
        
        resultados_multiframework['frameworks']['qiskit'] = {
            'status': 'success',
            'resultados': [resultado],
            'melhor_acuracia': float(acc_test)
        }
        
        logger.info(f"  ✅ Qiskit: Acurácia={acc_test:.4f}, Tempo={tempo:.2f}s")
        
    except Exception as e:
        logger.error(f"  ❌ Erro no Qiskit: {e}")
        resultados_multiframework['frameworks']['qiskit'] = {
            'status': 'error',
            'erro': str(e)
        }

# ============================================================================
# 3. FRAMEWORK CIRQ
# ============================================================================
if 'cirq' in config['enabled_frameworks']:
    logger.info("\n" + "="*80)
    logger.info("3. EXECUTANDO FRAMEWORK CIRQ")
    logger.info("="*80)
    
    try:
        from framework_cirq import (
            ClassificadorVQCCirq,
            gerar_dataset_sintetico
        )
        from sklearn.model_selection import train_test_split
        
        # Gerar dataset sintético reduzido
        X, y = gerar_dataset_sintetico(n_samples=50, n_features=2)
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.3, random_state=config['default_seed']
        )
        
        config_cirq = {
            'ansatz': 'strongly_entangling',
            'modelo_ruido': 'phase_damping',
            'nivel_ruido': 0.005
        }
        
        logger.info(f"\n  Configuração: {config_cirq}")
        
        inicio = time.time()
        
        vqc = ClassificadorVQCCirq(
            n_qubits=4,
            n_camadas=2,
            ansatz=config_cirq['ansatz'],
            modelo_ruido=config_cirq['modelo_ruido'],
            nivel_ruido=config_cirq['nivel_ruido'],
            n_epocas=5,
            n_shots=256,  # Reduzido
            seed=config['default_seed']
        )
        
        vqc.fit(X_train[:30], y_train[:30])
        
        y_pred_test = vqc.predict(X_test[:15])
        
        acc_test = np.mean(y_pred_test == y_test[:15])
        tempo = time.time() - inicio
        
        resultado = {
            'framework': 'cirq',
            'arquitetura': config_cirq['ansatz'],
            'tipo_ruido': config_cirq['modelo_ruido'],
            'nivel_ruido': config_cirq['nivel_ruido'],
            'acuracia_teste': float(acc_test),
            'tempo_segundos': float(tempo)
        }
        
        resultados_multiframework['frameworks']['cirq'] = {
            'status': 'success',
            'resultados': [resultado],
            'melhor_acuracia': float(acc_test)
        }
        
        logger.info(f"  ✅ Cirq: Acurácia={acc_test:.4f}, Tempo={tempo:.2f}s")
        
    except Exception as e:
        logger.error(f"  ❌ Erro no Cirq: {e}")
        resultados_multiframework['frameworks']['cirq'] = {
            'status': 'error',
            'erro': str(e)
        }

# ============================================================================
# 4. CONSOLIDAÇÃO E ANÁLISE COMPARATIVA
# ============================================================================
logger.info("\n" + "="*80)
logger.info("4. CONSOLIDANDO RESULTADOS")
logger.info("="*80)

# Consolidar todos os resultados
todos_resultados = []
for fw_nome, fw_dados in resultados_multiframework['frameworks'].items():
    if fw_dados['status'] == 'success':
        todos_resultados.extend(fw_dados['resultados'])

if todos_resultados:
    df_resultados = pd.DataFrame(todos_resultados)
    
    # Salvar CSV
    csv_path = os.path.join(pasta_resultados, 'resultados_multiframework.csv')
    df_resultados.to_csv(csv_path, index=False)
    logger.info(f"\n  ✓ CSV salvo: {csv_path}")
    
    # Estatísticas comparativas
    logger.info("\n  📊 ESTATÍSTICAS COMPARATIVAS:")
    
    for framework in df_resultados['framework'].unique():
        df_fw = df_resultados[df_resultados['framework'] == framework]
        logger.info(f"\n    {framework.upper()}:")
        logger.info(f"      • Acurácia: {df_fw['acuracia_teste'].iloc[0]:.4f}")
        logger.info(f"      • Tempo: {df_fw['tempo_segundos'].iloc[0]:.2f}s")
    
    # Comparação entre frameworks
    logger.info("\n  🏆 RANKING DE FRAMEWORKS:")
    ranking = df_resultados.sort_values('acuracia_teste', ascending=False)
    for i, row in enumerate(ranking.itertuples(), 1):
        logger.info(f"    {i}. {row.framework.upper()}: {row.acuracia_teste:.4f}")
    
    resultados_multiframework['comparacao'] = {
        'melhor_framework': ranking.iloc[0]['framework'],
        'melhor_acuracia_geral': float(ranking.iloc[0]['acuracia_teste']),
        'todos_frameworks': ranking[['framework', 'acuracia_teste']].to_dict('records')
    }

# ============================================================================
# 5. SALVAR RESULTADOS
# ============================================================================
logger.info("\n" + "="*80)
logger.info("5. SALVANDO RESULTADOS")
logger.info("="*80)

# Salvar JSON completo
json_path = os.path.join(pasta_resultados, 'resultados_completos.json')
with open(json_path, 'w', encoding='utf-8') as f:
    json.dump(resultados_multiframework, f, indent=2, ensure_ascii=False)
logger.info(f"\n  ✓ JSON salvo: {json_path}")

# Gerar manifesto de execução (Task 5 do MegaPrompt)
manifesto = {
    'data_execucao': timestamp,
    'versao_framework': config['version'],
    'seed_utilizada': config['default_seed'],
    'frameworks_executados': list(resultados_multiframework['frameworks'].keys()),
    'total_experimentos': len(todos_resultados) if todos_resultados else 0,
    'pasta_resultados': pasta_resultados,
    'python_version': sys.version,
    'modo_execucao': 'rapido',
    'bibliotecas': {
        'numpy': np.__version__,
        'pandas': pd.__version__
    }
}

manifesto_path = os.path.join(pasta_resultados, 'execution_manifest.json')
with open(manifesto_path, 'w', encoding='utf-8') as f:
    json.dump(manifesto, f, indent=2, ensure_ascii=False)
logger.info(f"  ✓ Manifesto de execução salvo: {manifesto_path}")

# ============================================================================
# 6. RESUMO FINAL
# ============================================================================
logger.info("\n" + "="*80)
logger.info("✅ EXECUÇÃO MULTIFRAMEWORK COMPLETA")
logger.info("="*80)

logger.info(f"\n📁 Resultados salvos em: {pasta_resultados}/")
logger.info(f"\n📊 Resumo:")
logger.info(f"  • Frameworks executados: {len([f for f in resultados_multiframework['frameworks'].values() if f['status'] == 'success'])}/{len(config['enabled_frameworks'])}")
logger.info(f"  • Total de experimentos: {len(todos_resultados) if todos_resultados else 0}")

if 'comparacao' in resultados_multiframework and resultados_multiframework['comparacao']:
    logger.info(f"  • Melhor framework: {resultados_multiframework['comparacao']['melhor_framework'].upper()}")
    logger.info(f"  • Melhor acurácia: {resultados_multiframework['comparacao']['melhor_acuracia_geral']:.4f}")

logger.info("\n" + "="*80)
logger.info("Os resultados estão prontos para serem usados na atualização da documentação.")
logger.info("="*80)

# Retornar código de sucesso
sys.exit(0)
