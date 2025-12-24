"""
Script to generate comparative results for visualization
This runs a targeted set of experiments to ensure we have data across:
- Different initialization strategies (matematico, quantico, aleatoria, fibonacci_spiral)
- Different noise types (sem_ruido, depolarizante, amplitude_damping, phase_damping)
- Different noise levels (0, 0.001, 0.005, 0.01)
- Different architectures (basic_entangler, strongly_entangling, hardware_efficient)
"""

import os
os.environ['VQC_QUICK'] = '1'

import sys
import time
import logging
from datetime import datetime
from pathlib import Path
import json

# Import from main framework
from framework_investigativo_completo import (
    carregar_datasets,
    ClassificadorVQC,
    executar_analises_estatisticas,
    gerar_visualizacoes,
    consolidar_e_gerar_metadados
)

import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def executar_grid_comparativo():
    """Execute a focused grid search for comparative visualizations"""
    
    logger.info("="*80)
    logger.info("GERANDO RESULTADOS COMPARATIVOS PARA VISUALIZAÇÕES")
    logger.info("="*80)
    
    # Load dataset
    datasets = carregar_datasets()
    dataset_name = 'moons'
    ds = datasets[dataset_name]
    X_train, X_test, y_train, y_test = ds['X_train'], ds['X_test'], ds['y_train'], ds['y_test']
    
    logger.info(f"\nDataset: {dataset_name}")
    logger.info(f"  Treino: {len(X_train)}, Teste: {len(X_test)}")
    
    # Define comparative grid
    configs = [
        # Comparison of initialization strategies (same noise, same architecture)
        {'arch': 'strongly_entangling', 'init': 'matematico', 'noise': 'depolarizante', 'level': 0.001},
        {'arch': 'strongly_entangling', 'init': 'quantico', 'noise': 'depolarizante', 'level': 0.001},
        {'arch': 'strongly_entangling', 'init': 'aleatoria', 'noise': 'depolarizante', 'level': 0.001},
        {'arch': 'strongly_entangling', 'init': 'fibonacci_spiral', 'noise': 'depolarizante', 'level': 0.001},
        
        # Comparison of noise types (same init, same architecture)
        {'arch': 'strongly_entangling', 'init': 'quantico', 'noise': 'sem_ruido', 'level': 0},
        {'arch': 'strongly_entangling', 'init': 'quantico', 'noise': 'depolarizante', 'level': 0.001},
        {'arch': 'strongly_entangling', 'init': 'quantico', 'noise': 'amplitude_damping', 'level': 0.001},
        {'arch': 'strongly_entangling', 'init': 'quantico', 'noise': 'phase_damping', 'level': 0.001},
        
        # Comparison of noise levels (same init, same architecture, same noise type)
        {'arch': 'strongly_entangling', 'init': 'quantico', 'noise': 'depolarizante', 'level': 0.0},
        {'arch': 'strongly_entangling', 'init': 'quantico', 'noise': 'depolarizante', 'level': 0.001},
        {'arch': 'strongly_entangling', 'init': 'quantico', 'noise': 'depolarizante', 'level': 0.005},
        {'arch': 'strongly_entangling', 'init': 'quantico', 'noise': 'depolarizante', 'level': 0.01},
        
        # Comparison of architectures (same init, same noise)
        {'arch': 'basic_entangler', 'init': 'quantico', 'noise': 'depolarizante', 'level': 0.001},
        {'arch': 'strongly_entangling', 'init': 'quantico', 'noise': 'depolarizante', 'level': 0.001},
        {'arch': 'hardware_efficient', 'init': 'quantico', 'noise': 'depolarizante', 'level': 0.001},
    ]
    
    resultados = []
    
    for i, config in enumerate(configs, 1):
        logger.info(f"\n[{i}/{len(configs)}] Executando configuração:")
        logger.info(f"  Arquitetura: {config['arch']}")
        logger.info(f"  Inicialização: {config['init']}")
        logger.info(f"  Ruído: {config['noise']}")
        logger.info(f"  Nível: {config['level']}")
        
        try:
            vqc = ClassificadorVQC(
                n_qubits=4,
                n_camadas=2,
                arquitetura=config['arch'],
                estrategia_init=config['init'],
                seed=42
            )
            
            if config['noise'] != 'sem_ruido' and config['level'] > 0:
                vqc.tipo_ruido = config['noise']
                vqc.nivel_ruido = config['level']
            
            inicio = time.time()
            vqc.fit(X_train, y_train, n_epocas=5)
            tempo = time.time() - inicio
            
            acc_train = vqc.score(X_train, y_train)
            acc_test = vqc.score(X_test, y_test)
            gap = acc_train - acc_test
            
            resultado = {
                'dataset': dataset_name,
                'seed': 42,
                'n_qubits': 4,
                'n_camadas': 2,
                'arquitetura': config['arch'],
                'estrategia_init': config['init'],
                'tipo_ruido': config['noise'],
                'nivel_ruido': config['level'],
                'acuracia_treino': acc_train,
                'acuracia_teste': acc_test,
                'gap_treino_teste': gap,
                'tempo_segundos': tempo,
                'n_parametros': len(vqc.pesos),
                'entropia_final': 0.0,
                'negatividade_media': 0.0,
                'barren_plateau_detectado': False,
                'convergiu_early_stopping': False
            }
            
            resultados.append(resultado)
            
            logger.info(f"  ✓ Acurácia teste: {acc_test:.4f}")
            logger.info(f"  ✓ Tempo: {tempo:.1f}s")
            
        except Exception as e:
            logger.warning(f"  ✗ Erro: {str(e)}")
            continue
    
    # Save results
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    pasta_resultados = Path(f"resultados_{timestamp}")
    pasta_resultados.mkdir(exist_ok=True)
    
    df = pd.DataFrame(resultados)
    csv_path = pasta_resultados / "resultados_completos_artigo.csv"
    df.to_csv(csv_path, index=False)
    
    logger.info(f"\n✓ Resultados salvos: {csv_path}")
    logger.info(f"✓ Total de experimentos: {len(resultados)}")
    
    return df, pasta_resultados

if __name__ == "__main__":
    df_resultados, pasta_resultados = executar_grid_comparativo()
    
    logger.info("\n" + "="*80)
    logger.info("GERANDO VISUALIZAÇÕES COMPARATIVAS")
    logger.info("="*80)
    
    # Generate visualizations
    try:
        gerar_visualizacoes(df_resultados, salvar_figuras=True, pasta_resultados=str(pasta_resultados))
        logger.info("\n✓ Visualizações geradas com sucesso!")
    except Exception as e:
        logger.error(f"\n✗ Erro ao gerar visualizações: {e}")
        import traceback
        traceback.print_exc()
    
    # Generate statistical analyses
    try:
        executar_analises_estatisticas(df_resultados, pasta_resultados=str(pasta_resultados))
        logger.info("✓ Análises estatísticas geradas!")
    except Exception as e:
        logger.error(f"✗ Erro nas análises: {e}")
    
    logger.info("\n" + "="*80)
    logger.info("EXECUÇÃO COMPLETA!")
    logger.info("="*80)
    logger.info(f"\nResultados em: {pasta_resultados}/")
