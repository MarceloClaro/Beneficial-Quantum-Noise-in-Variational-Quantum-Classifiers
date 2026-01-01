#!/usr/bin/env python3
# =============================================================================
# SCRIPT DE EXECUÇÃO: QAOA 100 QUBITS COM RUÍDO BENÉFICO
# =============================================================================
"""
Script para executar experimentos QAOA com 100 qubits e análise de ruído benéfico.

Uso:
    python executar_qaoa_100qubits.py --modo demo
    python executar_qaoa_100qubits.py --modo completo --n_qubits 100
    python executar_qaoa_100qubits.py --modo rapido --p_layers 2

Autor: Framework QAOA
Data: 2025-12-26
"""

import os
import sys
import time
import logging
from pathlib import Path

# Importar framework QAOA
from framework_qaoa_100qubits import (
    ConfigQAOA,
    ConstrutorCircuitoQAOA,
    OtimizadorQAOA,
    AnalisadorHiperparametrosQAOA,
    VisualizadorQAOA,
    demo_qaoa_100qubits,
    experimento_completo_ruido_benefico,
    QISKIT_AVAILABLE
)

logger = logging.getLogger(__name__)


def executar_demo_rapida():
    """Demonstração rápida com parâmetros reduzidos para teste."""
    print("\n" + "="*80)
    print("DEMONSTRAÇÃO RÁPIDA: QAOA COM RUÍDO BENÉFICO")
    print("="*80)
    print("\nConfiguração:")
    print("  - Qubits: 20 (reduzido para demo)")
    print("  - P-layers: 2")
    print("  - Densidade grafo: 0.2")
    print("  - Tipo ruído: depolarizing")
    print("  - Nível ruído: 0.001")
    print("\nExecutando...\n")
    
    # Criar grafo menor para demo rápida
    construtor = ConstrutorCircuitoQAOA(n_qubits=20, p_layers=2, seed=42)
    grafo = construtor.criar_grafo_aleatorio(densidade=0.2)
    
    # Testar sem ruído
    print("\n--- Experimento 1: SEM RUÍDO ---")
    config_sem_ruido = ConfigQAOA(
        n_qubits=20,
        p_layers=2,
        tipo_ruido='sem_ruido',
        nivel_ruido=0.0,
        max_iter=50,
        shots=512
    )
    
    otimizador_sem = OtimizadorQAOA(config_sem_ruido)
    resultado_sem = otimizador_sem.otimizar(grafo)
    
    print(f"Energia final: {resultado_sem.energia_final:.4f}")
    print(f"Tempo: {resultado_sem.tempo_execucao:.2f}s")
    print(f"Iterações: {resultado_sem.iteracoes}")
    
    # Testar com ruído
    print("\n--- Experimento 2: COM RUÍDO DEPOLARIZANTE (0.001) ---")
    config_com_ruido = ConfigQAOA(
        n_qubits=20,
        p_layers=2,
        tipo_ruido='depolarizing',
        nivel_ruido=0.001,
        max_iter=50,
        shots=512
    )
    
    otimizador_com = OtimizadorQAOA(config_com_ruido)
    resultado_com = otimizador_com.otimizar(grafo)
    
    print(f"Energia final: {resultado_com.energia_final:.4f}")
    print(f"Tempo: {resultado_com.tempo_execucao:.2f}s")
    print(f"Iterações: {resultado_com.iteracoes}")
    
    # Comparação
    print("\n" + "-"*80)
    print("ANÁLISE COMPARATIVA:")
    print("-"*80)
    
    melhoria = ((resultado_sem.energia_final - resultado_com.energia_final) / 
                resultado_sem.energia_final * 100)
    
    print(f"Energia sem ruído:  {resultado_sem.energia_final:.4f}")
    print(f"Energia com ruído:  {resultado_com.energia_final:.4f}")
    print(f"Diferença relativa: {melhoria:+.2f}%")
    
    if melhoria > 0:
        print("✅ RUÍDO BENÉFICO DETECTADO! Energia melhorou com ruído.")
    else:
        print("⚠️  Ruído prejudicou o resultado nesta configuração.")
    
    print("\n" + "="*80)
    print("DEMONSTRAÇÃO CONCLUÍDA")
    print("="*80 + "\n")
    
    return {
        'sem_ruido': resultado_sem,
        'com_ruido': resultado_com,
        'melhoria_percentual': melhoria
    }


def executar_grid_search_pequeno():
    """Grid search em escala reduzida para teste rápido."""
    print("\n" + "="*80)
    print("GRID SEARCH: ANÁLISE DE RUÍDO BENÉFICO")
    print("="*80)
    print("\nConfiguração:")
    print("  - Qubits: 30")
    print("  - P-layers: 2")
    print("  - Níveis ruído: [0.0, 0.0005, 0.001, 0.002]")
    print("  - Tipos ruído: [sem_ruido, depolarizing, phase_damping]")
    print("  - Repetições: 2")
    print("\nExecutando...\n")
    
    # Criar grafo
    construtor = ConstrutorCircuitoQAOA(n_qubits=30, p_layers=2, seed=42)
    grafo = construtor.criar_grafo_aleatorio(densidade=0.15)
    
    # Analisador
    analisador = AnalisadorHiperparametrosQAOA(
        pasta_resultados='resultados_qaoa_grid_pequeno'
    )
    
    # Grid search
    df = analisador.grid_search_ruido(
        grafo=grafo,
        niveis_ruido=[0.0, 0.0005, 0.001, 0.002],
        tipos_ruido=['sem_ruido', 'depolarizing', 'phase_damping'],
        p_layers=2,
        n_repeticoes=2
    )
    
    # Análise
    print("\n" + "-"*80)
    print("RESULTADOS POR TIPO DE RUÍDO (Energia média):")
    print("-"*80)
    
    resumo = df.groupby(['tipo_ruido', 'nivel_ruido'])['energia_final'].agg(['mean', 'std'])
    print(resumo)
    
    # Encontrar melhor configuração
    melhor_idx = df['energia_final'].idxmin()
    melhor = df.loc[melhor_idx]
    
    print("\n" + "-"*80)
    print("MELHOR CONFIGURAÇÃO:")
    print("-"*80)
    print(f"Tipo ruído:    {melhor['tipo_ruido']}")
    print(f"Nível ruído:   {melhor['nivel_ruido']:.4f}")
    print(f"Energia final: {melhor['energia_final']:.4f}")
    print(f"Tempo exec:    {melhor['tempo_execucao']:.2f}s")
    
    # Visualizar
    visualizador = VisualizadorQAOA()
    pasta_resultados = Path('resultados_qaoa_grid_pequeno')
    visualizador.plotar_comparacao_ruido(
        df,
        salvar=str(pasta_resultados / 'comparacao_ruido.html')
    )
    
    print(f"\n✅ Visualização salva em: {pasta_resultados / 'comparacao_ruido.html'}")
    
    print("\n" + "="*80)
    print("GRID SEARCH CONCLUÍDO")
    print("="*80 + "\n")
    
    return df


def executar_qaoa_100qubits_completo():
    """Execução completa com 100 qubits (mais demorado)."""
    print("\n" + "="*80)
    print("EXPERIMENTO COMPLETO: QAOA 100 QUBITS")
    print("="*80)
    print("\n⚠️  ATENÇÃO: Esta execução pode levar várias horas!")
    print("Recomenda-se começar com modos 'rapido' ou 'grid' para testes.\n")
    
    resposta = input("Deseja continuar? (s/N): ")
    if resposta.lower() != 's':
        print("Execução cancelada.")
        return None
    
    print("\nIniciando experimento completo...\n")
    
    resultados = experimento_completo_ruido_benefico(
        n_qubits=100,
        densidade_grafo=0.1,
        p_layers=3
    )
    
    print("\n" + "="*80)
    print("EXPERIMENTO COMPLETO CONCLUÍDO")
    print("="*80 + "\n")
    
    return resultados


def executar_demonstracao_niveis_ruido():
    """Demonstra impacto de diferentes níveis de ruído."""
    print("\n" + "="*80)
    print("DEMONSTRAÇÃO: IMPACTO DE DIFERENTES NÍVEIS DE RUÍDO")
    print("="*80)
    
    # Configuração
    n_qubits = 25
    p_layers = 2
    niveis = [0.0, 0.0001, 0.0005, 0.001, 0.002, 0.005]
    
    print(f"\nConfiguração:")
    print(f"  - Qubits: {n_qubits}")
    print(f"  - P-layers: {p_layers}")
    print(f"  - Tipo ruído: depolarizing")
    print(f"  - Níveis testados: {niveis}")
    print("\nExecutando...\n")
    
    # Criar grafo
    construtor = ConstrutorCircuitoQAOA(n_qubits=n_qubits, p_layers=p_layers)
    grafo = construtor.criar_grafo_aleatorio(densidade=0.2)
    
    resultados = []
    
    for nivel in niveis:
        print(f"\n--- Testando nível {nivel:.4f} ---")
        
        config = ConfigQAOA(
            n_qubits=n_qubits,
            p_layers=p_layers,
            tipo_ruido='depolarizing' if nivel > 0 else 'sem_ruido',
            nivel_ruido=nivel,
            max_iter=50,
            shots=512
        )
        
        otimizador = OtimizadorQAOA(config)
        resultado = otimizador.otimizar(grafo)
        
        resultados.append({
            'nivel_ruido': nivel,
            'energia': resultado.energia_final,
            'tempo': resultado.tempo_execucao,
            'iteracoes': resultado.iteracoes
        })
        
        print(f"Energia: {resultado.energia_final:.4f}")
        print(f"Tempo: {resultado.tempo_execucao:.2f}s")
    
    # Análise
    print("\n" + "="*80)
    print("RESUMO COMPARATIVO")
    print("="*80)
    print(f"\n{'Nível Ruído':<15} {'Energia':<12} {'Tempo (s)':<12} {'Iterações':<12}")
    print("-"*80)
    
    for r in resultados:
        print(f"{r['nivel_ruido']:<15.4f} {r['energia']:<12.4f} "
              f"{r['tempo']:<12.2f} {r['iteracoes']:<12}")
    
    # Encontrar região de ruído benéfico
    energias = [r['energia'] for r in resultados]
    melhor_idx = energias.index(min(energias))
    melhor_nivel = resultados[melhor_idx]['nivel_ruido']
    
    print("\n" + "-"*80)
    print(f"✅ Melhor nível de ruído: {melhor_nivel:.4f}")
    print(f"   Energia obtida: {resultados[melhor_idx]['energia']:.4f}")
    
    if melhor_nivel > 0:
        print("\n🎯 RUÍDO BENÉFICO CONFIRMADO!")
        print(f"   Um nível moderado de ruído ({melhor_nivel:.4f}) melhorou os resultados.")
    else:
        print("\n⚠️  Neste experimento, ruído não foi benéfico.")
    
    print("\n" + "="*80 + "\n")
    
    return resultados


def main():
    """Função principal com menu interativo."""
    if not QISKIT_AVAILABLE:
        print("\n❌ ERRO: Qiskit não está instalado!")
        print("\nPara instalar, execute:")
        print("  pip install qiskit qiskit-aer")
        print("\nOu use o arquivo requirements.txt:")
        print("  pip install -r requirements.txt\n")
        sys.exit(1)
    
    print("\n" + "="*80)
    print(" "*20 + "FRAMEWORK QAOA 100 QUBITS")
    print(" "*15 + "Análise de Ruído Quântico Benéfico")
    print("="*80)
    
    if len(sys.argv) > 1:
        # Modo CLI
        modo = sys.argv[1] if len(sys.argv) > 1 else 'rapido'
    else:
        # Menu interativo
        print("\nModos de execução disponíveis:\n")
        print("  1. rapido       - Demonstração rápida (20 qubits, ~2 min)")
        print("  2. grid         - Grid search reduzido (30 qubits, ~15 min)")
        print("  3. niveis       - Teste de níveis de ruído (25 qubits, ~10 min)")
        print("  4. completo     - Experimento completo 100 qubits (LONGO!)")
        print("  5. demo         - Demo original do framework")
        print("\n  0. Sair")
        
        escolha = input("\nEscolha um modo (1-5): ").strip()
        
        modos_map = {
            '1': 'rapido',
            '2': 'grid',
            '3': 'niveis',
            '4': 'completo',
            '5': 'demo',
            '0': 'sair'
        }
        
        modo = modos_map.get(escolha, 'rapido')
        
        if modo == 'sair':
            print("\nEncerrando...\n")
            sys.exit(0)
    
    # Executar modo selecionado
    inicio_total = time.time()
    
    try:
        if modo == 'rapido':
            resultados = executar_demo_rapida()
        
        elif modo == 'grid':
            resultados = executar_grid_search_pequeno()
        
        elif modo == 'niveis':
            resultados = executar_demonstracao_niveis_ruido()
        
        elif modo == 'completo':
            resultados = executar_qaoa_100qubits_completo()
        
        elif modo == 'demo':
            resultados = demo_qaoa_100qubits(
                densidade_grafo=0.1,
                p_layers=3,
                tipo_ruido='depolarizing',
                nivel_ruido=0.001
            )
        
        else:
            print(f"\n❌ Modo '{modo}' não reconhecido.")
            print("Modos válidos: rapido, grid, niveis, completo, demo\n")
            sys.exit(1)
        
        tempo_total = time.time() - inicio_total
        
        print("\n" + "="*80)
        print(f"EXECUÇÃO CONCLUÍDA EM {tempo_total:.1f} SEGUNDOS")
        print("="*80 + "\n")
        
        return resultados
    
    except KeyboardInterrupt:
        print("\n\n⚠️  Execução interrompida pelo usuário.\n")
        sys.exit(0)
    
    except Exception as e:
        print(f"\n\n❌ ERRO durante execução:")
        print(f"   {str(e)}\n")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
