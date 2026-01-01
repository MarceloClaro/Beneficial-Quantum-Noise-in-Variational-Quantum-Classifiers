"""
Script otimizado para execução de 2 horas com geração de visualizações.

Estratégia:
- Executar subset otimizado de experimentos (~100-150 experimentos)
- Gerar visualizações de alta qualidade para cada configuração importante
- Focar em configurações que demonstram ruído benéfico
- Salvar todas as imagens em pasta organizada
"""

import os
import sys
import time
import logging
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path

from framework_qiskit import (
    ClassificadorVQCQiskit,
    carregar_datasets,
    visualizar_bloch_sphere,
    visualizar_state_city_3d,
    visualizar_qsphere,
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# ============================================================================
# CONFIGURAÇÃO OTIMIZADA PARA 2 HORAS
# ============================================================================

# Seleção estratégica para demonstrar ruído benéfico
DATASETS_SELECIONADOS = ['moons', 'circles', 'iris']  # 3 datasets

ARQUITETURAS_SELECIONADAS = [
    'basico',
    'strongly_entangling',
    'hardware_efficient',
]  # 3 arquiteturas

RUIDOS_SELECIONADOS = [
    'sem_ruido',
    'phase_damping',
    'amplitude_damping',
    'depolarizante',
]  # 4 tipos

NIVEIS_SELECIONADOS = [
    0.0,      # Baseline
    0.005,    # Ruído benéfico típico
    0.01,     # Limite benéfico
]  # 3 níveis

SEEDS_SELECIONADAS = [42, 43]  # 2 seeds

# Configurações
N_QUBITS = 4
N_CAMADAS = 2
N_EPOCAS = 5  # Reduzido para velocidade
SHOTS = 512   # Reduzido para velocidade

# Total: 3 × 3 × 4 × 3 × 2 = 216 experimentos (estimado ~30s cada = 1.8h)
# Mais tempo para visualizações = ~2h total


def gerar_visualizacoes_experimento(vqc, dataset, dataset_nome, config_nome, pasta_viz):
    """
    Gera todas as visualizações para um experimento.
    
    Args:
        vqc: Classificador treinado
        dataset: Dataset usado
        dataset_nome: Nome do dataset
        config_nome: Nome da configuração
        pasta_viz: Pasta para salvar visualizações
    """
    try:
        # Pegar exemplo do teste
        x_exemplo = dataset['X_test'][0]
        
        # 1. Esfera de Bloch
        bloch_path = os.path.join(pasta_viz, f'{config_nome}_bloch.png')
        visualizar_bloch_sphere(vqc, x_exemplo, bloch_path)
        logger.info(f"    ✓ Bloch sphere: {bloch_path}")
        
        # 2. State City 3D
        city_path = os.path.join(pasta_viz, f'{config_nome}_city3d.png')
        visualizar_state_city_3d(vqc, x_exemplo, city_path)
        logger.info(f"    ✓ State City 3D: {city_path}")
        
        # 3. Q-Sphere
        qsphere_path = os.path.join(pasta_viz, f'{config_nome}_qsphere.png')
        visualizar_qsphere(vqc, x_exemplo, qsphere_path)
        logger.info(f"    ✓ Q-Sphere: {qsphere_path}")
        
        # 4. Diagrama de circuito
        circuit_path = os.path.join(pasta_viz, f'{config_nome}_circuit.png')
        vqc.get_circuit_diagram(circuit_path)
        logger.info(f"    ✓ Circuito: {circuit_path}")
        
        return 4  # 4 visualizações geradas
    
    except Exception as e:
        logger.error(f"    ✗ Erro ao gerar visualizações: {e}")
        return 0


def executar_experimentos_com_visualizacoes():
    """
    Executa experimentos otimizados com geração de visualizações.
    """
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    pasta_resultados = f'resultados_qiskit_2h_{timestamp}'
    pasta_viz = os.path.join(pasta_resultados, 'visualizacoes')
    
    os.makedirs(pasta_resultados, exist_ok=True)
    os.makedirs(pasta_viz, exist_ok=True)
    
    logger.info("=" * 80)
    logger.info("FRAMEWORK QISKIT - EXECUÇÃO OTIMIZADA 2H COM VISUALIZAÇÕES")
    logger.info("=" * 80)
    logger.info(f"Datasets: {len(DATASETS_SELECIONADOS)}")
    logger.info(f"Arquiteturas: {len(ARQUITETURAS_SELECIONADAS)}")
    logger.info(f"Tipos de ruído: {len(RUIDOS_SELECIONADOS)}")
    logger.info(f"Níveis de ruído: {len(NIVEIS_SELECIONADOS)}")
    logger.info(f"Seeds: {len(SEEDS_SELECIONADAS)}")
    
    total_combinacoes = (len(DATASETS_SELECIONADOS) * len(ARQUITETURAS_SELECIONADAS) *
                        len(RUIDOS_SELECIONADOS) * len(NIVEIS_SELECIONADOS) * 
                        len(SEEDS_SELECIONADAS))
    logger.info(f"Total de experimentos: {total_combinacoes}")
    logger.info(f"Pasta de resultados: {pasta_resultados}")
    logger.info("=" * 80)
    
    # Carregar datasets
    logger.info("\nCarregando datasets...")
    all_datasets = carregar_datasets()
    logger.info("✓ Datasets carregados")
    
    # Listas para resultados
    resultados = []
    total_visualizacoes = 0
    experimento_count = 0
    inicio_total = time.time()
    
    # Loop principal
    for dataset_nome in DATASETS_SELECIONADOS:
        dataset = all_datasets[dataset_nome]
        
        for arquitetura in ARQUITETURAS_SELECIONADAS:
            for tipo_ruido in RUIDOS_SELECIONADOS:
                for nivel_ruido in NIVEIS_SELECIONADOS:
                    
                    # Pular combinação inválida
                    if tipo_ruido == 'sem_ruido' and nivel_ruido > 0:
                        continue
                    
                    for seed in SEEDS_SELECIONADAS:
                        experimento_count += 1
                        
                        # Nome da configuração
                        config_nome = f'{dataset_nome}_{arquitetura}_{tipo_ruido}_g{nivel_ruido:.4f}_s{seed}'
                        
                        logger.info(f"\n[{experimento_count}/{total_combinacoes}] {config_nome}")
                        
                        try:
                            inicio = time.time()
                            
                            # Criar e treinar VQC
                            vqc = ClassificadorVQCQiskit(
                                n_qubits=N_QUBITS,
                                n_camadas=N_CAMADAS,
                                arquitetura=arquitetura,
                                estrategia_init='matematico',
                                tipo_ruido=tipo_ruido,
                                nivel_ruido=nivel_ruido,
                                taxa_aprendizado=0.01,
                                n_epocas=N_EPOCAS,
                                batch_size=32,
                                seed=seed,
                                shots=SHOTS
                            )
                            
                            vqc.fit(dataset['X_train'], dataset['y_train'])
                            
                            # Avaliar
                            acc_treino = vqc.score(dataset['X_train'], dataset['y_train'])
                            acc_teste = vqc.score(dataset['X_test'], dataset['y_test'])
                            
                            tempo_treino = time.time() - inicio
                            
                            logger.info(f"  ✓ Acurácia treino: {acc_treino:.4f}")
                            logger.info(f"  ✓ Acurácia teste: {acc_teste:.4f}")
                            logger.info(f"  ✓ Tempo: {tempo_treino:.1f}s")
                            
                            # Gerar visualizações (a cada 2 experimentos para economizar tempo)
                            n_viz = 0
                            if experimento_count % 2 == 0:
                                logger.info(f"  📊 Gerando visualizações...")
                                n_viz = gerar_visualizacoes_experimento(
                                    vqc, dataset, dataset_nome, config_nome, pasta_viz
                                )
                                total_visualizacoes += n_viz
                            
                            # Salvar resultado
                            resultado = {
                                'experimento': experimento_count,
                                'config_nome': config_nome,
                                'dataset': dataset_nome,
                                'arquitetura': arquitetura,
                                'tipo_ruido': tipo_ruido,
                                'nivel_ruido': nivel_ruido,
                                'seed': seed,
                                'n_qubits': N_QUBITS,
                                'n_camadas': N_CAMADAS,
                                'n_epocas': N_EPOCAS,
                                'acuracia_treino': acc_treino,
                                'acuracia_teste': acc_teste,
                                'tempo_treino': tempo_treino,
                                'visualizacoes_geradas': n_viz,
                                'framework': 'qiskit',
                                'shots': SHOTS
                            }
                            
                            resultados.append(resultado)
                            
                            # Checkpoint a cada 10 experimentos
                            if experimento_count % 10 == 0:
                                df_temp = pd.DataFrame(resultados)
                                csv_temp = os.path.join(pasta_resultados, 'resultados_parciais.csv')
                                df_temp.to_csv(csv_temp, index=False)
                                
                                tempo_decorrido = time.time() - inicio_total
                                tempo_restante = (tempo_decorrido / experimento_count) * (total_combinacoes - experimento_count)
                                
                                logger.info(f"\n  💾 Checkpoint: {experimento_count} experimentos")
                                logger.info(f"  ⏱️ Tempo decorrido: {tempo_decorrido/60:.1f} min")
                                logger.info(f"  ⏱️ Tempo estimado restante: {tempo_restante/60:.1f} min")
                                logger.info(f"  📊 Visualizações geradas: {total_visualizacoes}")
                        
                        except Exception as e:
                            logger.error(f"  ✗ Erro: {e}")
                            # Registrar erro mas continuar
                            resultado = {
                                'experimento': experimento_count,
                                'config_nome': config_nome,
                                'dataset': dataset_nome,
                                'arquitetura': arquitetura,
                                'tipo_ruido': tipo_ruido,
                                'nivel_ruido': nivel_ruido,
                                'seed': seed,
                                'acuracia_teste': np.nan,
                                'erro': str(e)[:100]
                            }
                            resultados.append(resultado)
    
    # Finalizar
    tempo_total = time.time() - inicio_total
    
    logger.info("\n" + "=" * 80)
    logger.info("EXECUÇÃO COMPLETA")
    logger.info("=" * 80)
    logger.info(f"✓ Experimentos realizados: {len(resultados)}")
    logger.info(f"✓ Visualizações geradas: {total_visualizacoes}")
    logger.info(f"✓ Tempo total: {tempo_total/60:.1f} minutos")
    logger.info(f"✓ Tempo médio/experimento: {tempo_total/len(resultados):.1f}s")
    
    # Salvar resultados finais
    df_final = pd.DataFrame(resultados)
    csv_final = os.path.join(pasta_resultados, 'resultados_completos.csv')
    df_final.to_csv(csv_final, index=False)
    logger.info(f"✓ Resultados salvos: {csv_final}")
    
    # Estatísticas
    logger.info("\nESTATÍSTICAS:")
    df_valido = df_final[df_final['acuracia_teste'].notna()]
    logger.info(f"  Acurácia média: {df_valido['acuracia_teste'].mean():.4f}")
    logger.info(f"  Acurácia máxima: {df_valido['acuracia_teste'].max():.4f}")
    logger.info(f"  Acurácia mínima: {df_valido['acuracia_teste'].min():.4f}")
    logger.info(f"  Desvio padrão: {df_valido['acuracia_teste'].std():.4f}")
    
    # Comparar ruído benéfico
    logger.info("\nCOMPARAÇÃO RUÍDO:")
    for tipo_ruido in RUIDOS_SELECIONADOS:
        df_ruido = df_valido[df_valido['tipo_ruido'] == tipo_ruido]
        if len(df_ruido) > 0:
            logger.info(f"  {tipo_ruido}: {df_ruido['acuracia_teste'].mean():.4f} "
                       f"(±{df_ruido['acuracia_teste'].std():.4f})")
    
    # Melhor configuração
    if len(df_valido) > 0:
        idx_melhor = df_valido['acuracia_teste'].idxmax()
        melhor = df_valido.loc[idx_melhor]
        logger.info("\n🏆 MELHOR CONFIGURAÇÃO:")
        logger.info(f"  {melhor['config_nome']}")
        logger.info(f"  Acurácia: {melhor['acuracia_teste']:.4f}")
    
    # Listar visualizações geradas
    logger.info(f"\n📊 VISUALIZAÇÕES GERADAS:")
    viz_files = sorted(Path(pasta_viz).glob('*.png'))
    logger.info(f"  Total de arquivos: {len(viz_files)}")
    if viz_files:
        logger.info(f"  Pasta: {pasta_viz}")
        for i, f in enumerate(viz_files[:10], 1):
            logger.info(f"    {i}. {f.name}")
        if len(viz_files) > 10:
            logger.info(f"    ... e mais {len(viz_files) - 10} arquivos")
    
    logger.info("=" * 80)
    logger.info("✓ Execução completa! Verifique a pasta de resultados.")
    logger.info("=" * 80)
    
    return df_final, pasta_resultados


if __name__ == "__main__":
    try:
        logger.info("Iniciando execução otimizada de 2 horas...")
        df, pasta = executar_experimentos_com_visualizacoes()
        
        print(f"\n✓ Sucesso! Resultados em: {pasta}")
        print(f"✓ {len(df)} experimentos realizados")
        print(f"✓ Visualizações em: {pasta}/visualizacoes/")
        
        sys.exit(0)
    
    except KeyboardInterrupt:
        print("\n\n⚠️ Execução interrompida pelo usuário")
        sys.exit(1)
    
    except Exception as e:
        print(f"\n❌ Erro: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
