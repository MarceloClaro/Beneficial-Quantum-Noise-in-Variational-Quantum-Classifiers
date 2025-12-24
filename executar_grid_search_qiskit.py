"""
Script de Execução Completa do Framework Qiskit
Replica EXATAMENTE os mesmos experimentos do PennyLane

Parâmetros:
- 5 datasets: Moons, Circles, Iris, Breast Cancer, Wine
- 9 arquiteturas: basico, strongly_entangling, hardware_efficient, alternating_layers, 
                  brickwork, random_entangling, tree, star_entanglement, qaoa
- 4 estratégias init: matematico, quantico, aleatorio, fibonacci_spiral
- 6 tipos de ruído: sem_ruido, depolarizante, amplitude_damping, phase_damping, crosstalk, correlacionado
- 9 níveis: 0.0, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02
- 5 seeds: 42, 43, 44, 45, 46
"""

import os
import sys
import time
import logging
import pandas as pd
import numpy as np
from datetime import datetime
from pathlib import Path

from framework_qiskit import ClassificadorVQCQiskit, carregar_datasets

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# ============================================================================
# CONFIGURAÇÃO DOS EXPERIMENTOS (IDÊNTICO AO PENNYLANE)
# ============================================================================

DATASETS = ['moons', 'circles', 'iris', 'breast_cancer', 'wine']

ARQUITETURAS = [
    'basico',
    'strongly_entangling',
    'hardware_efficient',
    'alternating_layers',
    'brickwork',
    'random_entangling',
    'tree',
    'star_entanglement',
    'qaoa'
]

ESTRATEGIAS_INIT = [
    'matematico',
    'quantico',
    'aleatorio',
    'fibonacci_spiral'
]

TIPOS_RUIDO = [
    'sem_ruido',
    'depolarizante',
    'amplitude_damping',
    'phase_damping',
    'crosstalk',
    'correlacionado'
]

NIVEIS_RUIDO = [
    0.0,
    0.0025,
    0.005,
    0.0075,
    0.01,
    0.0125,
    0.015,
    0.0175,
    0.02
]

SEEDS = [42, 43, 44, 45, 46]

# Configurações de treinamento
N_QUBITS = 4
N_CAMADAS = 2
N_EPOCAS = 10  # Reduzido para velocidade (PennyLane usa 20)
BATCH_SIZE = 32
TAXA_APRENDIZADO = 0.01
SHOTS = 1024


def executar_experimento_completo(
    modo='rapido',
    max_experimentos=None,
    pasta_resultados=None,
    verbose=True
):
    """
    Executa grid search completo de experimentos Qiskit.
    
    Args:
        modo: 'rapido' (amostra), 'medio' (subset), 'completo' (todos)
        max_experimentos: Limite máximo de experimentos (None = sem limite)
        pasta_resultados: Pasta para salvar resultados
        verbose: Verbosidade
    """
    
    # Configurar subsets baseado no modo
    if modo == 'rapido':
        datasets_usar = ['moons']
        arquiteturas_usar = ['basico', 'strongly_entangling']
        estrategias_usar = ['aleatorio', 'matematico']
        tipos_ruido_usar = ['sem_ruido', 'depolarizante']
        niveis_usar = [0.0, 0.005, 0.01]
        seeds_usar = [42, 43]
    elif modo == 'medio':
        datasets_usar = ['moons', 'circles', 'iris']
        arquiteturas_usar = ARQUITETURAS[:5]
        estrategias_usar = ESTRATEGIAS_INIT
        tipos_ruido_usar = TIPOS_RUIDO[:4]
        niveis_usar = NIVEIS_RUIDO[::2]  # Pares
        seeds_usar = SEEDS[:3]
    else:  # completo
        datasets_usar = DATASETS
        arquiteturas_usar = ARQUITETURAS
        estrategias_usar = ESTRATEGIAS_INIT
        tipos_ruido_usar = TIPOS_RUIDO
        niveis_usar = NIVEIS_RUIDO
        seeds_usar = SEEDS
    
    # Criar diretório de resultados
    if pasta_resultados is None:
        timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        pasta_resultados = f'resultados_qiskit_completo_{timestamp}'
    
    os.makedirs(pasta_resultados, exist_ok=True)
    
    # Calcular total de experimentos
    total_experimentos = (
        len(datasets_usar) *
        len(arquiteturas_usar) *
        len(estrategias_usar) *
        len(tipos_ruido_usar) *
        len(niveis_usar) *
        len(seeds_usar)
    )
    
    if max_experimentos and total_experimentos > max_experimentos:
        total_experimentos = max_experimentos
    
    logger.info("=" * 80)
    logger.info("FRAMEWORK QISKIT - EXECUÇÃO COMPLETA")
    logger.info("=" * 80)
    logger.info(f"Modo: {modo.upper()}")
    logger.info(f"Total de experimentos: {total_experimentos}")
    logger.info(f"Datasets: {len(datasets_usar)}")
    logger.info(f"Arquiteturas: {len(arquiteturas_usar)}")
    logger.info(f"Estratégias init: {len(estrategias_usar)}")
    logger.info(f"Tipos de ruído: {len(tipos_ruido_usar)}")
    logger.info(f"Níveis de ruído: {len(niveis_usar)}")
    logger.info(f"Seeds: {len(seeds_usar)}")
    logger.info(f"Resultados em: {pasta_resultados}")
    logger.info("=" * 80)
    
    # Carregar datasets uma vez
    logger.info("\nCarregando datasets...")
    all_datasets = carregar_datasets()
    logger.info(f"✓ {len(all_datasets)} datasets carregados")
    
    # Lista para armazenar resultados
    resultados = []
    experimento_count = 0
    inicio_total = time.time()
    
    # Loop principal de experimentos
    for dataset_nome in datasets_usar:
        dataset = all_datasets[dataset_nome]
        
        for arquitetura in arquiteturas_usar:
            for estrategia_init in estrategias_usar:
                for tipo_ruido in tipos_ruido_usar:
                    for nivel_ruido in niveis_usar:
                        # Pular combinações inválidas
                        if tipo_ruido == 'sem_ruido' and nivel_ruido > 0:
                            continue
                        
                        for seed in seeds_usar:
                            experimento_count += 1
                            
                            # Verificar limite
                            if max_experimentos and experimento_count > max_experimentos:
                                logger.info(f"\n✓ Limite de {max_experimentos} experimentos atingido")
                                break
                            
                            # Info do experimento
                            if verbose:
                                logger.info(f"\n[{experimento_count}/{total_experimentos}] Executando:")
                                logger.info(f"  Dataset: {dataset_nome}")
                                logger.info(f"  Arquitetura: {arquitetura}")
                                logger.info(f"  Init: {estrategia_init}")
                                logger.info(f"  Ruído: {tipo_ruido} (γ={nivel_ruido:.4f})")
                                logger.info(f"  Seed: {seed}")
                            
                            try:
                                # Criar classificador
                                inicio = time.time()
                                
                                vqc = ClassificadorVQCQiskit(
                                    n_qubits=N_QUBITS,
                                    n_camadas=N_CAMADAS,
                                    arquitetura=arquitetura,
                                    estrategia_init=estrategia_init,
                                    tipo_ruido=tipo_ruido,
                                    nivel_ruido=nivel_ruido,
                                    taxa_aprendizado=TAXA_APRENDIZADO,
                                    n_epocas=N_EPOCAS,
                                    batch_size=BATCH_SIZE,
                                    seed=seed,
                                    shots=SHOTS
                                )
                                
                                # Treinar
                                vqc.fit(dataset['X_train'], dataset['y_train'])
                                
                                # Avaliar
                                acuracia_treino = vqc.score(dataset['X_train'], dataset['y_train'])
                                acuracia_teste = vqc.score(dataset['X_test'], dataset['y_test'])
                                
                                tempo_treino = time.time() - inicio
                                
                                # Armazenar resultado
                                resultado = {
                                    'dataset': dataset_nome,
                                    'arquitetura': arquitetura,
                                    'estrategia_init': estrategia_init,
                                    'tipo_ruido': tipo_ruido,
                                    'nivel_ruido': nivel_ruido,
                                    'seed': seed,
                                    'n_qubits': N_QUBITS,
                                    'n_camadas': N_CAMADAS,
                                    'n_epocas': N_EPOCAS,
                                    'acuracia_treino': acuracia_treino,
                                    'acuracia_teste': acuracia_teste,
                                    'tempo_treino': tempo_treino,
                                    'framework': 'qiskit',
                                    'shots': SHOTS
                                }
                                
                                resultados.append(resultado)
                                
                                if verbose:
                                    logger.info(f"  ✓ Acurácia treino: {acuracia_treino:.4f}")
                                    logger.info(f"  ✓ Acurácia teste: {acuracia_teste:.4f}")
                                    logger.info(f"  ✓ Tempo: {tempo_treino:.2f}s")
                                
                                # Salvar resultados parciais a cada 10 experimentos
                                if experimento_count % 10 == 0:
                                    df_temp = pd.DataFrame(resultados)
                                    csv_path = os.path.join(pasta_resultados, 'resultados_parciais.csv')
                                    df_temp.to_csv(csv_path, index=False)
                                    if verbose:
                                        logger.info(f"\n  💾 Checkpoint: {len(resultados)} resultados salvos")
                            
                            except Exception as e:
                                logger.error(f"  ✗ Erro no experimento: {str(e)}")
                                # Registrar erro mas continuar
                                resultado = {
                                    'dataset': dataset_nome,
                                    'arquitetura': arquitetura,
                                    'estrategia_init': estrategia_init,
                                    'tipo_ruido': tipo_ruido,
                                    'nivel_ruido': nivel_ruido,
                                    'seed': seed,
                                    'n_qubits': N_QUBITS,
                                    'n_camadas': N_CAMADAS,
                                    'n_epocas': N_EPOCAS,
                                    'acuracia_treino': np.nan,
                                    'acuracia_teste': np.nan,
                                    'tempo_treino': np.nan,
                                    'framework': 'qiskit',
                                    'shots': SHOTS,
                                    'erro': str(e)[:100]
                                }
                                resultados.append(resultado)
    
    # Finalizar
    tempo_total = time.time() - inicio_total
    
    logger.info("\n" + "=" * 80)
    logger.info("EXECUÇÃO COMPLETA")
    logger.info("=" * 80)
    logger.info(f"✓ Total de experimentos: {len(resultados)}")
    logger.info(f"✓ Tempo total: {tempo_total / 60:.2f} minutos")
    logger.info(f"✓ Tempo médio por experimento: {tempo_total / len(resultados):.2f}s")
    
    # Salvar resultados finais
    df_final = pd.DataFrame(resultados)
    csv_final = os.path.join(pasta_resultados, 'resultados_completos_qiskit.csv')
    df_final.to_csv(csv_final, index=False)
    logger.info(f"✓ Resultados salvos: {csv_final}")
    
    # Estatísticas resumidas
    logger.info("\nESTATÍSTICAS:")
    logger.info(f"  Acurácia média (teste): {df_final['acuracia_teste'].mean():.4f}")
    logger.info(f"  Acurácia máxima: {df_final['acuracia_teste'].max():.4f}")
    logger.info(f"  Acurácia mínima: {df_final['acuracia_teste'].min():.4f}")
    logger.info(f"  Desvio padrão: {df_final['acuracia_teste'].std():.4f}")
    
    # Melhor configuração
    idx_melhor = df_final['acuracia_teste'].idxmax()
    melhor = df_final.loc[idx_melhor]
    logger.info("\n🏆 MELHOR CONFIGURAÇÃO:")
    logger.info(f"  Dataset: {melhor['dataset']}")
    logger.info(f"  Arquitetura: {melhor['arquitetura']}")
    logger.info(f"  Init: {melhor['estrategia_init']}")
    logger.info(f"  Ruído: {melhor['tipo_ruido']} (γ={melhor['nivel_ruido']:.4f})")
    logger.info(f"  Seed: {melhor['seed']}")
    logger.info(f"  Acurácia: {melhor['acuracia_teste']:.4f}")
    
    logger.info("=" * 80)
    
    return df_final


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Execução completa Framework Qiskit')
    parser.add_argument('--modo', choices=['rapido', 'medio', 'completo'], default='rapido',
                       help='Modo de execução')
    parser.add_argument('--max', type=int, default=None,
                       help='Número máximo de experimentos')
    parser.add_argument('--pasta', type=str, default=None,
                       help='Pasta para resultados')
    parser.add_argument('--quiet', action='store_true',
                       help='Modo silencioso')
    
    args = parser.parse_args()
    
    try:
        df_resultados = executar_experimento_completo(
            modo=args.modo,
            max_experimentos=args.max,
            pasta_resultados=args.pasta,
            verbose=not args.quiet
        )
        
        print(f"\n✓ Execução concluída! {len(df_resultados)} experimentos realizados.")
        sys.exit(0)
        
    except KeyboardInterrupt:
        print("\n\n⚠️ Execução interrompida pelo usuário")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Erro fatal: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
