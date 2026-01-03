#!/usr/bin/env python3
"""
QUALIS A1 Enhanced Execution
Executa framework_quantum_advanced_v8.py com análise otimizada para publicação A1
"""

import os
import sys
import json
import time
import logging
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-8s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def import_framework():
    """Importa o framework do v8"""
    try:
        from framework_quantum_advanced_v8 import (
            AdvancedVQC, 
            CircuitArchitecture,
            NoiseModel,
            DatasetLoader,
            DeepChemDatasetLoader,
            ErrorMitigationType,
            AdvancedConfig
        )
        return {
            'AdvancedVQC': AdvancedVQC,
            'CircuitArchitecture': CircuitArchitecture,
            'NoiseModel': NoiseModel,
            'DatasetLoader': DatasetLoader,
            'DeepChemDatasetLoader': DeepChemDatasetLoader,
            'ErrorMitigationType': ErrorMitigationType,
            'AdvancedConfig': AdvancedConfig
        }
    except ImportError as e:
        logger.error(f"✗ Erro ao importar framework: {e}")
        return None

def load_datasets(framework_module):
    """Carrega todos os datasets disponíveis"""
    datasets = {}
    loader = framework_module['DatasetLoader']()
    deepchem_loader = framework_module['DeepChemDatasetLoader']()
    
    # Datasets sklearn
    sklearn_datasets = [
        'IRIS', 'WINE', 'BREAST_CANCER', 'DIGITS', 'DIABETES'
    ]
    
    # Datasets DeepChem
    deepchem_datasets = ['BACE', 'HIV']
    
    logger.info("=" * 80)
    logger.info("CARREGANDO DATASETS")
    logger.info("=" * 80)
    
    for dataset_name in sklearn_datasets:
        try:
            X, y = loader.load(dataset_name)
            datasets[dataset_name] = {'X': X, 'y': y, 'source': 'sklearn'}
            logger.info(f"✓ {dataset_name:20} | Amostras: {len(X):6} | Features: {X.shape[1]:2}")
        except Exception as e:
            logger.warning(f"✗ {dataset_name:20} | Erro: {str(e)[:50]}")
    
    for dataset_name in deepchem_datasets:
        try:
            X, y = deepchem_loader.load(dataset_name)
            datasets[dataset_name] = {'X': X, 'y': y, 'source': 'deepchem'}
            logger.info(f"✓ {dataset_name:20} | Amostras: {len(X):6} | Features: {X.shape[1]:2}")
        except Exception as e:
            logger.warning(f"✗ {dataset_name:20} | Erro: {str(e)[:50]}")
    
    logger.info(f"\n✓ {len(datasets)} datasets carregados com sucesso")
    return datasets

def execute_qualis_a1_optimized(framework_module, datasets):
    """Executa versão otimizada para QUALIS A1"""
    
    AdvancedVQC = framework_module['AdvancedVQC']
    CircuitArchitecture = framework_module['CircuitArchitecture']
    NoiseModel = framework_module['NoiseModel']
    AdvancedConfig = framework_module['AdvancedConfig']
    
    results = []
    
    # Configuração otimizada para publicação
    config = AdvancedConfig(
        n_qubits=6,
        layers=3,
        learning_rate=0.01,
        iterations=100,
        seed=42
    )
    
    logger.info("=" * 80)
    logger.info("EXECUÇÃO OTIMIZADA PARA QUALIS A1")
    logger.info("=" * 80)
    logger.info(f"Qubits: {config.n_qubits}")
    logger.info(f"Camadas: {config.layers}")
    logger.info(f"Taxa de aprendizado: {config.learning_rate}")
    logger.info(f"Iterações: {config.iterations}")
    logger.info(f"Seed: {config.seed}")
    logger.info("")
    
    circuits_to_test = [
        CircuitArchitecture.STRONGLY_ENTANGLING,
        CircuitArchitecture.REAL_AMPLITUDES,
        CircuitArchitecture.EFFICIENT_SU2,
        CircuitArchitecture.HARDWARE_EFFICIENT
    ]
    
    noise_models_to_test = [
        NoiseModel.DEPOLARIZING,
        NoiseModel.AMPLITUDE_DAMPING,
        NoiseModel.PHASE_DAMPING,
        NoiseModel.BIT_FLIP
    ]
    
    datasets_to_test = ['IRIS', 'WINE', 'BREAST_CANCER', 'DIGITS']
    
    total_experiments = len(circuits_to_test) * len(noise_models_to_test) * len(datasets_to_test)
    current = 0
    
    logger.info(f"Total de experimentos: {total_experiments}")
    logger.info("")
    
    for dataset_name in datasets_to_test:
        if dataset_name not in datasets:
            logger.warning(f"Dataset {dataset_name} não disponível")
            continue
        
        X, y = datasets[dataset_name]['X'], datasets[dataset_name]['y']
        
        for circuit in circuits_to_test:
            for noise in noise_models_to_test:
                current += 1
                
                try:
                    # Log do experimento
                    logger.info(f"[{current:3}/{total_experiments}] Executando: {dataset_name:15} | "
                              f"{str(circuit.value):25} | {str(noise.value):20}")
                    
                    # Criar e treinar modelo
                    vqc = AdvancedVQC(config, circuit, noise)
                    start_time = time.time()
                    
                    # Train/test split
                    from sklearn.model_selection import train_test_split
                    X_train, X_test, y_train, y_test = train_test_split(
                        X, y, test_size=0.2, random_state=42
                    )
                    
                    # Treinar
                    vqc.fit(X_train, y_train)
                    
                    # Avaliar
                    train_score = vqc.score(X_train, y_train)
                    test_score = vqc.score(X_test, y_test)
                    elapsed = time.time() - start_time
                    
                    results.append({
                        'dataset': dataset_name,
                        'circuit': str(circuit.value),
                        'noise_model': str(noise.value),
                        'train_accuracy': train_score,
                        'test_accuracy': test_score,
                        'time_seconds': elapsed,
                        'n_samples': len(X),
                        'n_features': X.shape[1]
                    })
                    
                    logger.info(f"           → Train: {train_score:.4f} | Test: {test_score:.4f} | "
                              f"Tempo: {elapsed:.2f}s")
                    
                except Exception as e:
                    logger.error(f"           → Erro: {str(e)[:60]}")
    
    return pd.DataFrame(results)

def generate_publication_tables(results_df):
    """Gera tabelas para publicação"""
    
    logger.info("\n" + "=" * 80)
    logger.info("TABELAS DE RESULTADOS PARA PUBLICAÇÃO")
    logger.info("=" * 80)
    
    # Tabela 1: Melhores resultados por dataset
    logger.info("\n### Tabela 1: Top 5 Melhores Configurações")
    logger.info("-" * 80)
    
    top_5 = results_df.nlargest(5, 'test_accuracy')[
        ['dataset', 'circuit', 'noise_model', 'test_accuracy', 'train_accuracy', 'time_seconds']
    ]
    
    for idx, row in top_5.iterrows():
        logger.info(f"{row['dataset']:15} | {row['circuit']:25} | {row['noise_model']:20} | "
                   f"Test: {row['test_accuracy']:.4f} | Train: {row['train_accuracy']:.4f} | "
                   f"Tempo: {row['time_seconds']:.2f}s")
    
    # Tabela 2: Estatísticas por dataset
    logger.info("\n### Tabela 2: Estatísticas por Dataset")
    logger.info("-" * 80)
    
    dataset_stats = results_df.groupby('dataset').agg({
        'test_accuracy': ['mean', 'std', 'min', 'max'],
        'time_seconds': 'mean'
    })
    
    for dataset in dataset_stats.index:
        stats = dataset_stats.loc[dataset]
        logger.info(f"{dataset:15} | Média: {stats[('test_accuracy', 'mean')]:.4f} ± "
                   f"{stats[('test_accuracy', 'std')]:.4f} | "
                   f"Range: [{stats[('test_accuracy', 'min')]:.4f}, "
                   f"{stats[('test_accuracy', 'max')]:.4f}]")
    
    # Tabela 3: Performance por tipo de circuito
    logger.info("\n### Tabela 3: Performance por Tipo de Circuito")
    logger.info("-" * 80)
    
    circuit_stats = results_df.groupby('circuit')['test_accuracy'].agg(['mean', 'std', 'count'])
    
    for circuit in circuit_stats.index:
        stats = circuit_stats.loc[circuit]
        logger.info(f"{circuit:25} | Média: {stats['mean']:.4f} ± {stats['std']:.4f} | "
                   f"Experimentos: {int(stats['count'])}")
    
    # Tabela 4: Impacto do modelo de ruído
    logger.info("\n### Tabela 4: Impacto do Modelo de Ruído")
    logger.info("-" * 80)
    
    noise_stats = results_df.groupby('noise_model')['test_accuracy'].agg(['mean', 'std', 'count'])
    
    for noise in noise_stats.index:
        stats = noise_stats.loc[noise]
        logger.info(f"{noise:20} | Média: {stats['mean']:.4f} ± {stats['std']:.4f} | "
                   f"Experimentos: {int(stats['count'])}")
    
    return {
        'top_5': top_5,
        'dataset_stats': dataset_stats,
        'circuit_stats': circuit_stats,
        'noise_stats': noise_stats
    }

def save_results(results_df, tables_dict):
    """Salva resultados em arquivo"""
    
    output_dir = Path('resultados_qualis_a1')
    output_dir.mkdir(exist_ok=True)
    
    # Salvar CSV
    csv_path = output_dir / 'qualis_a1_results.csv'
    results_df.to_csv(csv_path, index=False)
    logger.info(f"\n✓ Resultados salvos em: {csv_path}")
    
    # Salvar markdown
    md_path = output_dir / 'QUALIS_A1_RESULTS.md'
    
    with open(md_path, 'w', encoding='utf-8') as f:
        f.write("# Framework Quantum Advanced V8 - Resultados QUALIS A1\n\n")
        f.write(f"Data: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## Resumo Executivo\n\n")
        f.write(f"- Total de experimentos: {len(results_df)}\n")
        f.write(f"- Acurácia média: {results_df['test_accuracy'].mean():.4f}\n")
        f.write(f"- Melhor resultado: {results_df['test_accuracy'].max():.4f}\n")
        f.write(f"- Tempo total: {results_df['time_seconds'].sum():.2f}s\n\n")
        
        f.write("## Top 5 Configurações\n\n")
        f.write(tables_dict['top_5'].to_markdown(index=False))
        
        f.write("\n\n## Estatísticas por Dataset\n\n")
        f.write(str(tables_dict['dataset_stats']))
        
        f.write("\n\n## Performance por Circuito\n\n")
        f.write(str(tables_dict['circuit_stats']))
        
        f.write("\n\n## Impacto do Ruído\n\n")
        f.write(str(tables_dict['noise_stats']))
        
        f.write("\n\n## Tabela Completa de Resultados\n\n")
        f.write(results_df.to_markdown(index=False))
    
    logger.info(f"✓ Markdown salvo em: {md_path}")
    
    return output_dir

def main():
    """Função principal"""
    
    logger.info("")
    logger.info("=" * 80)
    logger.info("FRAMEWORK QUANTUM ADVANCED V8 - EXECUÇÃO OTIMIZADA QUALIS A1")
    logger.info("=" * 80)
    logger.info("")
    
    # Importar framework
    framework_module = import_framework()
    if not framework_module:
        logger.error("Falha ao importar framework")
        sys.exit(1)
    
    # Carregar datasets
    datasets = load_datasets(framework_module)
    if not datasets:
        logger.error("Nenhum dataset carregado")
        sys.exit(1)
    
    # Executar experimentos
    results_df = execute_qualis_a1_optimized(framework_module, datasets)
    
    if results_df.empty:
        logger.error("Nenhum resultado foi obtido")
        sys.exit(1)
    
    # Gerar tabelas
    tables_dict = generate_publication_tables(results_df)
    
    # Salvar resultados
    output_dir = save_results(results_df, tables_dict)
    
    logger.info("")
    logger.info("=" * 80)
    logger.info("✓ EXECUÇÃO QUALIS A1 CONCLUÍDA COM SUCESSO!")
    logger.info("=" * 80)
    logger.info(f"Resultados salvos em: {output_dir}")
    logger.info("")

if __name__ == '__main__':
    main()
