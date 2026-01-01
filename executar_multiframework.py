"""
Script para executar o multiframework (PennyLane, Qiskit, Cirq) e gerar resultados comparativos.
Aplica as melhorias do MegaPrompt Especializado para QUALIS A1.

Este script executa experimentos focados em todos os frameworks e atualiza os resultados.
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
logger.info("EXECUÇÃO MULTIFRAMEWORK - QUALIS A1 ENHANCED")
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
        
        # Configuração focada para teste rápido
        configs_pennylane = [
            {
                'arquitetura': 'strongly_entangling',
                'tipo_ruido': 'phase_damping',
                'nivel_ruido': 0.005,
                'estrategia_init': 'quantico'
            },
            {
                'arquitetura': 'hardware_efficient',
                'tipo_ruido': 'depolarizante',
                'nivel_ruido': 0.01,
                'estrategia_init': 'matematico'
            },
            {
                'arquitetura': 'strongly_entangling',
                'tipo_ruido': 'sem_ruido',
                'nivel_ruido': 0.0,
                'estrategia_init': 'quantico'
            }
        ]
        
        resultados_pennylane = []
        
        for i, cfg in enumerate(configs_pennylane, 1):
            logger.info(f"\n  [{i}/{len(configs_pennylane)}] Config: {cfg}")
            
            inicio = time.time()
            
            vqc = ClassificadorVQC(
                n_qubits=4,
                n_camadas=2,
                arquitetura=cfg['arquitetura'],
                estrategia_init=cfg['estrategia_init'],
                tipo_ruido=cfg['tipo_ruido'],
                nivel_ruido=cfg['nivel_ruido'],
                n_epocas=10,  # Reduzido para execução rápida
                seed=config['default_seed']
            )
            
            vqc.fit(
                dataset_moons['X_train'],
                dataset_moons['y_train']
            )
            
            acc_train = vqc.score(dataset_moons['X_train'], dataset_moons['y_train'])
            acc_test = vqc.score(dataset_moons['X_test'], dataset_moons['y_test'])
            tempo = time.time() - inicio
            
            resultado = {
                'framework': 'pennylane',
                'arquitetura': cfg['arquitetura'],
                'tipo_ruido': cfg['tipo_ruido'],
                'nivel_ruido': cfg['nivel_ruido'],
                'estrategia_init': cfg['estrategia_init'],
                'acuracia_treino': float(acc_train),
                'acuracia_teste': float(acc_test),
                'tempo_segundos': float(tempo)
            }
            
            resultados_pennylane.append(resultado)
            logger.info(f"    ✓ Acurácia Teste: {acc_test:.4f} | Tempo: {tempo:.2f}s")
        
        resultados_multiframework['frameworks']['pennylane'] = {
            'status': 'success',
            'resultados': resultados_pennylane,
            'melhor_acuracia': max([r['acuracia_teste'] for r in resultados_pennylane])
        }
        
        logger.info(f"\n  ✅ PennyLane completo: {len(resultados_pennylane)} experimentos")
        
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
        
        configs_qiskit = [
            {
                'arquitetura': 'strongly_entangling',
                'tipo_ruido': 'phase_damping',
                'nivel_ruido': 0.005
            },
            {
                'arquitetura': 'hardware_efficient',
                'tipo_ruido': 'depolarizante',
                'nivel_ruido': 0.01
            },
            {
                'arquitetura': 'basico',
                'tipo_ruido': 'sem_ruido',
                'nivel_ruido': 0.0
            }
        ]
        
        resultados_qiskit = []
        
        for i, cfg in enumerate(configs_qiskit, 1):
            logger.info(f"\n  [{i}/{len(configs_qiskit)}] Config: {cfg}")
            
            inicio = time.time()
            
            vqc = ClassificadorVQCQiskit(
                n_qubits=4,
                n_camadas=2,
                arquitetura=cfg['arquitetura'],
                tipo_ruido=cfg['tipo_ruido'],
                nivel_ruido=cfg['nivel_ruido'],
                n_epocas=10,
                seed=config['default_seed']
            )
            
            vqc.fit(dataset_moons_q['X_train'], dataset_moons_q['y_train'])
            
            acc_train = vqc.score(dataset_moons_q['X_train'], dataset_moons_q['y_train'])
            acc_test = vqc.score(dataset_moons_q['X_test'], dataset_moons_q['y_test'])
            tempo = time.time() - inicio
            
            resultado = {
                'framework': 'qiskit',
                'arquitetura': cfg['arquitetura'],
                'tipo_ruido': cfg['tipo_ruido'],
                'nivel_ruido': cfg['nivel_ruido'],
                'acuracia_treino': float(acc_train),
                'acuracia_teste': float(acc_test),
                'tempo_segundos': float(tempo)
            }
            
            resultados_qiskit.append(resultado)
            logger.info(f"    ✓ Acurácia Teste: {acc_test:.4f} | Tempo: {tempo:.2f}s")
        
        resultados_multiframework['frameworks']['qiskit'] = {
            'status': 'success',
            'resultados': resultados_qiskit,
            'melhor_acuracia': max([r['acuracia_teste'] for r in resultados_qiskit])
        }
        
        logger.info(f"\n  ✅ Qiskit completo: {len(resultados_qiskit)} experimentos")
        
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
        
        # Gerar dataset sintético
        X, y = gerar_dataset_sintetico(n_samples=100, n_features=2)
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.3, random_state=config['default_seed']
        )
        
        configs_cirq = [
            {
                'ansatz': 'strongly_entangling',
                'modelo_ruido': 'phase_damping',
                'nivel_ruido': 0.005
            },
            {
                'ansatz': 'hardware_efficient',
                'modelo_ruido': 'depolarizante',
                'nivel_ruido': 0.01
            },
            {
                'ansatz': 'basic_entangler',
                'modelo_ruido': None,
                'nivel_ruido': 0.0
            }
        ]
        
        resultados_cirq = []
        
        for i, cfg in enumerate(configs_cirq, 1):
            logger.info(f"\n  [{i}/{len(configs_cirq)}] Config: {cfg}")
            
            inicio = time.time()
            
            vqc = ClassificadorVQCCirq(
                n_qubits=4,
                n_camadas=2,
                ansatz=cfg['ansatz'],
                modelo_ruido=cfg['modelo_ruido'],
                nivel_ruido=cfg['nivel_ruido'],
                n_epocas=10,
                n_shots=512,
                seed=config['default_seed']
            )
            
            vqc.fit(X_train, y_train)
            
            y_pred_train = vqc.predict(X_train)
            y_pred_test = vqc.predict(X_test)
            
            acc_train = np.mean(y_pred_train == y_train)
            acc_test = np.mean(y_pred_test == y_test)
            tempo = time.time() - inicio
            
            resultado = {
                'framework': 'cirq',
                'arquitetura': cfg['ansatz'],
                'tipo_ruido': cfg['modelo_ruido'] or 'sem_ruido',
                'nivel_ruido': cfg['nivel_ruido'],
                'acuracia_treino': float(acc_train),
                'acuracia_teste': float(acc_test),
                'tempo_segundos': float(tempo)
            }
            
            resultados_cirq.append(resultado)
            logger.info(f"    ✓ Acurácia Teste: {acc_test:.4f} | Tempo: {tempo:.2f}s")
        
        resultados_multiframework['frameworks']['cirq'] = {
            'status': 'success',
            'resultados': resultados_cirq,
            'melhor_acuracia': max([r['acuracia_teste'] for r in resultados_cirq])
        }
        
        logger.info(f"\n  ✅ Cirq completo: {len(resultados_cirq)} experimentos")
        
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

# Consolidar todos os resultados em um DataFrame
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
        logger.info(f"      • Acurácia média: {df_fw['acuracia_teste'].mean():.4f}")
        logger.info(f"      • Acurácia máxima: {df_fw['acuracia_teste'].max():.4f}")
        logger.info(f"      • Tempo médio: {df_fw['tempo_segundos'].mean():.2f}s")
        logger.info(f"      • Experimentos: {len(df_fw)}")
    
    # Análise de ruído benéfico
    logger.info("\n  🔬 ANÁLISE DE RUÍDO BENÉFICO:")
    
    for framework in df_resultados['framework'].unique():
        df_fw = df_resultados[df_resultados['framework'] == framework]
        
        sem_ruido = df_fw[df_fw['tipo_ruido'] == 'sem_ruido']
        com_ruido = df_fw[df_fw['tipo_ruido'] != 'sem_ruido']
        
        if len(sem_ruido) > 0 and len(com_ruido) > 0:
            acc_sem = sem_ruido['acuracia_teste'].mean()
            acc_com = com_ruido['acuracia_teste'].mean()
            delta = acc_com - acc_sem
            
            logger.info(f"\n    {framework.upper()}:")
            logger.info(f"      • Sem ruído: {acc_sem:.4f}")
            logger.info(f"      • Com ruído: {acc_com:.4f}")
            logger.info(f"      • Delta: {delta:+.4f} ({'+' if delta > 0 else ''}{delta*100:.2f}%)")
    
    # Comparação entre frameworks
    logger.info("\n  🏆 RANKING DE FRAMEWORKS:")
    ranking = df_resultados.groupby('framework')['acuracia_teste'].max().sort_values(ascending=False)
    for i, (fw, acc) in enumerate(ranking.items(), 1):
        logger.info(f"    {i}. {fw.upper()}: {acc:.4f}")
    
    resultados_multiframework['comparacao'] = {
        'melhor_framework': ranking.index[0],
        'melhor_acuracia_geral': float(ranking.iloc[0]),
        'ranking': ranking.to_dict()
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
logger.info(f"  • Frameworks executados: {len([f for f in resultados_multiframework['frameworks'].values() if f['status'] == 'success'])}")
logger.info(f"  • Total de experimentos: {len(todos_resultados) if todos_resultados else 0}")

if 'comparacao' in resultados_multiframework and resultados_multiframework['comparacao']:
    logger.info(f"  • Melhor framework: {resultados_multiframework['comparacao']['melhor_framework'].upper()}")
    logger.info(f"  • Melhor acurácia: {resultados_multiframework['comparacao']['melhor_acuracia_geral']:.4f}")

logger.info("\n" + "="*80)
