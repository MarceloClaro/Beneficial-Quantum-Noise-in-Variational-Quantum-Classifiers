"""
Exemplo de Uso Completo do Framework Qiskit v7.2
================================================

Este script demonstra como usar o framework Qiskit para replicar todos os
experimentos do artigo, incluindo:
- Treinamento de VQC com diferentes arquiteturas
- Análise de ruído quântico benéfico
- Visualizações de esfera de Bloch
- Diagramas de circuitos quânticos
- Gráficos 3D de estados quânticos
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path

# Adicionar diretório pai ao path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from framework_qiskit import (
    ClassificadorVQCQiskit,
    carregar_datasets,
    executar_experimento_qiskit,
    visualizar_bloch_sphere,
    visualizar_state_city_3d,
    visualizar_qsphere,
    ARQUITETURAS_QISKIT,
    MODELOS_RUIDO_QISKIT
)

import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def exemplo_1_experimento_basico():
    """Exemplo 1: Experimento básico com Qiskit."""
    logger.info("="*80)
    logger.info("EXEMPLO 1: Experimento Básico com Qiskit")
    logger.info("="*80)
    
    resultado = executar_experimento_qiskit(
        dataset_nome='moons',
        arquitetura='basico',
        tipo_ruido='depolarizante',
        nivel_ruido=0.01,
        n_epocas=15,
        pasta_resultados='resultados_exemplo_1',
        verbose=True
    )
    
    logger.info("\n✓ Experimento concluído!")
    logger.info(f"  Acurácia treino: {resultado['acuracia_treino']:.4f}")
    logger.info(f"  Acurácia teste: {resultado['acuracia_teste']:.4f}")
    logger.info(f"  Tempo de treino: {resultado['tempo_treino']:.2f}s")
    
    return resultado


def exemplo_2_comparar_arquiteturas():
    """Exemplo 2: Comparar diferentes arquiteturas de circuitos."""
    logger.info("\n" + "="*80)
    logger.info("EXEMPLO 2: Comparação de Arquiteturas")
    logger.info("="*80)
    
    arquiteturas = ['basico', 'strongly_entangling', 'hardware_efficient', 'brickwork']
    resultados = []
    
    datasets = carregar_datasets()
    dataset = datasets['moons']
    
    for arq in arquiteturas:
        logger.info(f"\nTestando arquitetura: {arq}")
        
        vqc = ClassificadorVQCQiskit(
            n_qubits=4,
            n_camadas=2,
            arquitetura=arq,
            tipo_ruido='sem_ruido',
            nivel_ruido=0.0,
            n_epocas=10,
            seed=42
        )
        
        vqc.fit(dataset['X_train'], dataset['y_train'])
        acuracia = vqc.score(dataset['X_test'], dataset['y_test'])
        
        resultados.append({
            'arquitetura': arq,
            'acuracia_teste': acuracia
        })
        
        logger.info(f"  Acurácia: {acuracia:.4f}")
    
    # Criar DataFrame
    df_resultados = pd.DataFrame(resultados)
    logger.info("\n" + "="*80)
    logger.info("RESUMO DAS ARQUITETURAS")
    logger.info("="*80)
    logger.info(f"\n{df_resultados.to_string(index=False)}")
    
    # Melhor arquitetura
    melhor = df_resultados.loc[df_resultados['acuracia_teste'].idxmax()]
    logger.info(f"\n🏆 Melhor arquitetura: {melhor['arquitetura']} ({melhor['acuracia_teste']:.4f})")
    
    return df_resultados


def exemplo_3_analise_ruido_benefico():
    """Exemplo 3: Análise de ruído quântico benéfico."""
    logger.info("\n" + "="*80)
    logger.info("EXEMPLO 3: Análise de Ruído Quântico Benéfico")
    logger.info("="*80)
    
    niveis_ruido = [0.0, 0.001, 0.005, 0.01, 0.02, 0.05]
    tipos_ruido = ['sem_ruido', 'depolarizante', 'amplitude_damping', 'phase_damping']
    
    resultados = []
    datasets = carregar_datasets()
    dataset = datasets['circles']
    
    for tipo in tipos_ruido:
        for nivel in niveis_ruido:
            if tipo == 'sem_ruido' and nivel > 0:
                continue
            
            logger.info(f"\nTestando: {tipo} - Nível: {nivel:.3f}")
            
            vqc = ClassificadorVQCQiskit(
                n_qubits=4,
                n_camadas=2,
                arquitetura='basico',
                tipo_ruido=tipo,
                nivel_ruido=nivel,
                n_epocas=12,
                seed=42
            )
            
            vqc.fit(dataset['X_train'], dataset['y_train'])
            acuracia_teste = vqc.score(dataset['X_test'], dataset['y_test'])
            acuracia_treino = vqc.score(dataset['X_train'], dataset['y_train'])
            
            resultados.append({
                'tipo_ruido': tipo,
                'nivel_ruido': nivel,
                'acuracia_teste': acuracia_teste,
                'acuracia_treino': acuracia_treino,
                'gap': acuracia_treino - acuracia_teste
            })
            
            logger.info(f"  Acurácia teste: {acuracia_teste:.4f}")
            logger.info(f"  Gap: {acuracia_treino - acuracia_teste:.4f}")
    
    # Criar DataFrame
    df_resultados = pd.DataFrame(resultados)
    
    # Salvar resultados
    os.makedirs('resultados_exemplo_3', exist_ok=True)
    csv_path = 'resultados_exemplo_3/analise_ruido_benefico.csv'
    df_resultados.to_csv(csv_path, index=False)
    
    logger.info("\n" + "="*80)
    logger.info("RESUMO DA ANÁLISE DE RUÍDO")
    logger.info("="*80)
    
    # Identificar melhor configuração
    melhor = df_resultados.loc[df_resultados['acuracia_teste'].idxmax()]
    logger.info(f"\n🏆 Melhor configuração:")
    logger.info(f"  Tipo de ruído: {melhor['tipo_ruido']}")
    logger.info(f"  Nível: {melhor['nivel_ruido']:.3f}")
    logger.info(f"  Acurácia teste: {melhor['acuracia_teste']:.4f}")
    
    # Verificar se há ruído benéfico
    sem_ruido_acc = df_resultados[df_resultados['tipo_ruido'] == 'sem_ruido']['acuracia_teste'].values[0]
    
    logger.info(f"\n📊 Comparação com baseline (sem ruído: {sem_ruido_acc:.4f}):")
    for tipo in ['depolarizante', 'amplitude_damping', 'phase_damping']:
        df_tipo = df_resultados[df_resultados['tipo_ruido'] == tipo]
        melhor_tipo = df_tipo.loc[df_tipo['acuracia_teste'].idxmax()]
        delta = melhor_tipo['acuracia_teste'] - sem_ruido_acc
        
        if delta > 0:
            logger.info(f"  ✓ {tipo}: {melhor_tipo['acuracia_teste']:.4f} (Δ=+{delta:.4f}) - BENÉFICO!")
        else:
            logger.info(f"  ✗ {tipo}: {melhor_tipo['acuracia_teste']:.4f} (Δ={delta:.4f})")
    
    logger.info(f"\n💾 Resultados salvos em: {csv_path}")
    
    return df_resultados


def exemplo_4_visualizacoes_completas():
    """Exemplo 4: Gerar todas as visualizações quânticas."""
    logger.info("\n" + "="*80)
    logger.info("EXEMPLO 4: Visualizações Completas (Bloch, 3D, Circuitos)")
    logger.info("="*80)
    
    # Criar diretório de saída
    pasta_viz = 'resultados_exemplo_4_visualizacoes'
    os.makedirs(pasta_viz, exist_ok=True)
    
    # Treinar modelo
    datasets = carregar_datasets()
    dataset = datasets['moons']
    
    logger.info("\nTreinando modelo...")
    vqc = ClassificadorVQCQiskit(
        n_qubits=4,
        n_camadas=2,
        arquitetura='strongly_entangling',
        tipo_ruido='phase_damping',
        nivel_ruido=0.005,
        n_epocas=15,
        seed=42
    )
    
    vqc.fit(dataset['X_train'], dataset['y_train'])
    acuracia = vqc.score(dataset['X_test'], dataset['y_test'])
    logger.info(f"✓ Modelo treinado - Acurácia: {acuracia:.4f}")
    
    # Selecionar exemplo
    x_exemplo = dataset['X_test'][0]
    logger.info(f"\nGerando visualizações para entrada: {x_exemplo}")
    
    # 1. Diagrama de circuito
    logger.info("\n1. Gerando diagrama de circuito...")
    circuit_path = os.path.join(pasta_viz, 'circuito_quantico_qiskit.png')
    vqc.get_circuit_diagram(circuit_path)
    logger.info(f"   ✓ Salvo: {circuit_path}")
    
    # 2. Esfera de Bloch
    logger.info("\n2. Gerando visualização da Esfera de Bloch...")
    bloch_path = os.path.join(pasta_viz, 'bloch_sphere_qiskit.png')
    visualizar_bloch_sphere(vqc, x_exemplo, bloch_path)
    logger.info(f"   ✓ Salvo: {bloch_path}")
    
    # 3. State City 3D
    logger.info("\n3. Gerando visualização 3D State City...")
    city_path = os.path.join(pasta_viz, 'state_city_3d_qiskit.png')
    visualizar_state_city_3d(vqc, x_exemplo, city_path)
    logger.info(f"   ✓ Salvo: {city_path}")
    
    # 4. Q-Sphere
    logger.info("\n4. Gerando visualização Q-Sphere...")
    qsphere_path = os.path.join(pasta_viz, 'qsphere_qiskit.png')
    visualizar_qsphere(vqc, x_exemplo, qsphere_path)
    logger.info(f"   ✓ Salvo: {qsphere_path}")
    
    logger.info("\n" + "="*80)
    logger.info("✓ TODAS AS VISUALIZAÇÕES GERADAS COM SUCESSO!")
    logger.info("="*80)
    logger.info(f"📁 Visualizações disponíveis em: {pasta_viz}/")
    logger.info("   - Diagrama de circuito quântico")
    logger.info("   - Esfera de Bloch (estados de qubit)")
    logger.info("   - State City 3D (densidade de probabilidade)")
    logger.info("   - Q-Sphere (representação esférica)")
    
    return pasta_viz


def exemplo_5_experimento_completo_multiplos_datasets():
    """Exemplo 5: Experimento completo com múltiplos datasets."""
    logger.info("\n" + "="*80)
    logger.info("EXEMPLO 5: Experimento Completo - Múltiplos Datasets")
    logger.info("="*80)
    
    datasets_nomes = ['moons', 'circles', 'iris']
    arquiteturas = ['basico', 'hardware_efficient', 'brickwork']
    tipos_ruido = ['sem_ruido', 'depolarizante', 'phase_damping']
    niveis_ruido = [0.0, 0.005, 0.01]
    
    resultados = []
    total_exp = len(datasets_nomes) * len(arquiteturas) * len(tipos_ruido) * len(niveis_ruido)
    exp_count = 0
    
    logger.info(f"\nExecutando {total_exp} experimentos...")
    logger.info("(Isso pode levar alguns minutos)\n")
    
    for ds_nome in datasets_nomes:
        datasets = carregar_datasets()
        dataset = datasets[ds_nome]
        
        for arq in arquiteturas:
            for tipo_ruido in tipos_ruido:
                for nivel in niveis_ruido:
                    if tipo_ruido == 'sem_ruido' and nivel > 0:
                        continue
                    
                    exp_count += 1
                    logger.info(f"[{exp_count}/{total_exp}] {ds_nome} | {arq} | {tipo_ruido} | γ={nivel:.3f}")
                    
                    try:
                        vqc = ClassificadorVQCQiskit(
                            n_qubits=4,
                            n_camadas=2,
                            arquitetura=arq,
                            tipo_ruido=tipo_ruido,
                            nivel_ruido=nivel,
                            n_epocas=10,
                            seed=42
                        )
                        
                        vqc.fit(dataset['X_train'], dataset['y_train'])
                        acuracia_teste = vqc.score(dataset['X_test'], dataset['y_test'])
                        
                        resultados.append({
                            'dataset': ds_nome,
                            'arquitetura': arq,
                            'tipo_ruido': tipo_ruido,
                            'nivel_ruido': nivel,
                            'acuracia_teste': acuracia_teste
                        })
                        
                        logger.info(f"     Acurácia: {acuracia_teste:.4f}")
                        
                    except Exception as e:
                        logger.warning(f"     Erro: {str(e)[:50]}")
    
    # Salvar resultados
    df_resultados = pd.DataFrame(resultados)
    pasta_out = 'resultados_exemplo_5_completo'
    os.makedirs(pasta_out, exist_ok=True)
    csv_path = os.path.join(pasta_out, 'resultados_completos_qiskit.csv')
    df_resultados.to_csv(csv_path, index=False)
    
    logger.info("\n" + "="*80)
    logger.info("ANÁLISE FINAL")
    logger.info("="*80)
    
    # Estatísticas gerais
    logger.info(f"\n📊 Estatísticas Gerais:")
    logger.info(f"   Total de experimentos: {len(df_resultados)}")
    logger.info(f"   Acurácia média: {df_resultados['acuracia_teste'].mean():.4f}")
    logger.info(f"   Acurácia máxima: {df_resultados['acuracia_teste'].max():.4f}")
    logger.info(f"   Acurácia mínima: {df_resultados['acuracia_teste'].min():.4f}")
    
    # Melhor configuração geral
    melhor = df_resultados.loc[df_resultados['acuracia_teste'].idxmax()]
    logger.info(f"\n🏆 Melhor Configuração Global:")
    logger.info(f"   Dataset: {melhor['dataset']}")
    logger.info(f"   Arquitetura: {melhor['arquitetura']}")
    logger.info(f"   Ruído: {melhor['tipo_ruido']}")
    logger.info(f"   Nível: {melhor['nivel_ruido']:.3f}")
    logger.info(f"   Acurácia: {melhor['acuracia_teste']:.4f}")
    
    # Análise por dataset
    logger.info(f"\n📈 Desempenho por Dataset:")
    for ds in datasets_nomes:
        df_ds = df_resultados[df_resultados['dataset'] == ds]
        media = df_ds['acuracia_teste'].mean()
        melhor_ds = df_ds['acuracia_teste'].max()
        logger.info(f"   {ds:10s}: média={media:.4f}, melhor={melhor_ds:.4f}")
    
    logger.info(f"\n💾 Resultados completos salvos em: {csv_path}")
    
    return df_resultados


def main():
    """Função principal - executa todos os exemplos."""
    logger.info("="*80)
    logger.info("FRAMEWORK QISKIT v7.2 - EXEMPLOS DE USO COMPLETO")
    logger.info("Ruído Quântico Benéfico em Classificadores Variacionais")
    logger.info("="*80)
    
    try:
        # Verificar disponibilidade do Qiskit
        from framework_qiskit import QISKIT_AVAILABLE
        if not QISKIT_AVAILABLE:
            logger.error("\n❌ Qiskit não está instalado!")
            logger.info("Instale com: pip install qiskit qiskit-aer qiskit-ibm-runtime")
            return
        
        logger.info("\n✓ Qiskit está disponível - iniciando exemplos...\n")
        
        # Menu interativo
        logger.info("Selecione qual exemplo executar:")
        logger.info("  1 - Experimento Básico")
        logger.info("  2 - Comparar Arquiteturas")
        logger.info("  3 - Análise de Ruído Benéfico")
        logger.info("  4 - Visualizações Completas (Bloch, 3D, Circuitos)")
        logger.info("  5 - Experimento Completo (Múltiplos Datasets)")
        logger.info("  0 - Executar TODOS os exemplos")
        
        escolha = input("\nEscolha (0-5): ").strip()
        
        if escolha == '1':
            exemplo_1_experimento_basico()
        elif escolha == '2':
            exemplo_2_comparar_arquiteturas()
        elif escolha == '3':
            exemplo_3_analise_ruido_benefico()
        elif escolha == '4':
            exemplo_4_visualizacoes_completas()
        elif escolha == '5':
            exemplo_5_experimento_completo_multiplos_datasets()
        elif escolha == '0':
            logger.info("\n🚀 Executando TODOS os exemplos...\n")
            exemplo_1_experimento_basico()
            exemplo_2_comparar_arquiteturas()
            exemplo_3_analise_ruido_benefico()
            exemplo_4_visualizacoes_completas()
            exemplo_5_experimento_completo_multiplos_datasets()
        else:
            logger.warning("Escolha inválida!")
            return
        
        logger.info("\n" + "="*80)
        logger.info("✓ EXECUÇÃO CONCLUÍDA COM SUCESSO!")
        logger.info("="*80)
        logger.info("\nTodos os resultados e visualizações foram salvos.")
        logger.info("Confira os diretórios 'resultados_exemplo_*' para os arquivos gerados.")
        
    except KeyboardInterrupt:
        logger.info("\n\n⚠️ Execução interrompida pelo usuário.")
    except Exception as e:
        logger.error(f"\n❌ Erro durante execução: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
