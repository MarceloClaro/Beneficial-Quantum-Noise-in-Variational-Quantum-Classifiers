#!/usr/bin/env python3
"""
Exemplo Prático: QAOA com Análise de Ruído Benéfico

Este script demonstra como usar o framework QAOA para investigar
o fenômeno de ruído benéfico em algoritmos quânticos.

Experimentos:
1. Comparação direta: com vs. sem ruído
2. Varredura de níveis de ruído
3. Comparação entre tipos de ruído

Autor: Framework QAOA
Data: 2025-12-26
"""

import numpy as np
import pandas as pd
from pathlib import Path

# Verificar disponibilidade do framework
try:
    from framework_qaoa_100qubits import (
        ConfigQAOA,
        ConstrutorCircuitoQAOA,
        OtimizadorQAOA,
        QISKIT_AVAILABLE
    )
except ImportError as e:
    print(f"❌ Erro ao importar framework: {e}")
    print("\nVerifique se os arquivos estão no diretório correto:")
    print("  - framework_qaoa_100qubits.py")
    exit(1)

if not QISKIT_AVAILABLE:
    print("❌ Qiskit não está disponível!")
    print("\nPara instalar:")
    print("  pip install qiskit qiskit-aer")
    exit(1)


def exemplo_1_comparacao_basica():
    """
    Exemplo 1: Comparação básica entre QAOA com e sem ruído.
    
    Demonstra o fenômeno de ruído benéfico comparando a energia
    final obtida com e sem ruído quântico.
    """
    print("\n" + "="*80)
    print("EXEMPLO 1: COMPARAÇÃO BÁSICA - COM vs. SEM RUÍDO")
    print("="*80 + "\n")
    
    # Configuração
    n_qubits = 20
    p_layers = 2
    densidade = 0.3
    
    print(f"Configuração:")
    print(f"  • Qubits: {n_qubits}")
    print(f"  • P-layers: {p_layers}")
    print(f"  • Densidade grafo: {densidade}")
    print(f"  • Problema: MaxCut\n")
    
    # Criar grafo
    print("Criando grafo MaxCut...")
    construtor = ConstrutorCircuitoQAOA(n_qubits=n_qubits, p_layers=p_layers, seed=42)
    grafo = construtor.criar_grafo_aleatorio(densidade=densidade)
    n_arestas = np.sum(grafo > 0) // 2
    print(f"  ✓ Grafo criado: {n_arestas} arestas\n")
    
    # Experimento A: Sem ruído (baseline)
    print("--- Experimento A: SEM RUÍDO ---")
    config_sem = ConfigQAOA(
        n_qubits=n_qubits,
        p_layers=p_layers,
        tipo_ruido='sem_ruido',
        nivel_ruido=0.0,
        max_iter=50,
        shots=512,
        seed=42
    )
    
    otimizador_sem = OtimizadorQAOA(config_sem)
    resultado_sem = otimizador_sem.otimizar(grafo)
    
    print(f"  Energia final:     {resultado_sem.energia_final:.4f}")
    print(f"  Tempo execução:    {resultado_sem.tempo_execucao:.2f}s")
    print(f"  Iterações:         {resultado_sem.iteracoes}\n")
    
    # Experimento B: Com ruído depolarizante
    print("--- Experimento B: COM RUÍDO DEPOLARIZANTE (0.001) ---")
    config_com = ConfigQAOA(
        n_qubits=n_qubits,
        p_layers=p_layers,
        tipo_ruido='depolarizing',
        nivel_ruido=0.001,
        max_iter=50,
        shots=512,
        seed=42
    )
    
    otimizador_com = OtimizadorQAOA(config_com)
    resultado_com = otimizador_com.otimizar(grafo)
    
    print(f"  Energia final:     {resultado_com.energia_final:.4f}")
    print(f"  Tempo execução:    {resultado_com.tempo_execucao:.2f}s")
    print(f"  Iterações:         {resultado_com.iteracoes}\n")
    
    # Análise
    print("-"*80)
    print("ANÁLISE DOS RESULTADOS")
    print("-"*80 + "\n")
    
    diferenca_absoluta = resultado_com.energia_final - resultado_sem.energia_final
    diferenca_relativa = (diferenca_absoluta / abs(resultado_sem.energia_final)) * 100
    
    print(f"Energia sem ruído:      {resultado_sem.energia_final:.4f}")
    print(f"Energia com ruído:      {resultado_com.energia_final:.4f}")
    print(f"Diferença absoluta:     {diferenca_absoluta:+.4f}")
    print(f"Diferença relativa:     {diferenca_relativa:+.2f}%\n")
    
    if diferenca_relativa < -0.5:  # Melhoria > 0.5%
        print("✅ RUÍDO BENÉFICO DETECTADO!")
        print(f"   O ruído depolarizante (0.001) MELHOROU a energia em {abs(diferenca_relativa):.2f}%")
        print("   Isso demonstra que ruído moderado pode ajudar o algoritmo a")
        print("   escapar de mínimos locais e encontrar soluções melhores.\n")
    elif diferenca_relativa > 0.5:  # Piora > 0.5%
        print("⚠️  RUÍDO PREJUDICIAL neste caso")
        print(f"   O ruído piorou a energia em {abs(diferenca_relativa):.2f}%")
        print("   Nível de ruído pode ser muito alto para este problema.\n")
    else:
        print("➖ IMPACTO NEUTRO")
        print("   O ruído teve impacto mínimo (< 0.5%) neste caso.\n")
    
    return {
        'sem_ruido': resultado_sem,
        'com_ruido': resultado_com,
        'diferenca_relativa': diferenca_relativa
    }


def exemplo_2_varredura_niveis():
    """
    Exemplo 2: Varredura de níveis de ruído.
    
    Testa múltiplos níveis de ruído para encontrar a região
    de ruído benéfico ótima.
    """
    print("\n" + "="*80)
    print("EXEMPLO 2: VARREDURA DE NÍVEIS DE RUÍDO")
    print("="*80 + "\n")
    
    # Configuração
    n_qubits = 25
    p_layers = 2
    niveis = [0.0, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01]
    
    print(f"Configuração:")
    print(f"  • Qubits: {n_qubits}")
    print(f"  • P-layers: {p_layers}")
    print(f"  • Tipo ruído: depolarizing")
    print(f"  • Níveis testados: {len(niveis)}")
    print(f"  • Valores: {niveis}\n")
    
    # Criar grafo
    construtor = ConstrutorCircuitoQAOA(n_qubits=n_qubits, p_layers=p_layers, seed=42)
    grafo = construtor.criar_grafo_aleatorio(densidade=0.2)
    
    print("Executando varredura...\n")
    
    resultados = []
    
    for i, nivel in enumerate(niveis, 1):
        print(f"[{i}/{len(niveis)}] Testando nível {nivel:.4f}...")
        
        config = ConfigQAOA(
            n_qubits=n_qubits,
            p_layers=p_layers,
            tipo_ruido='depolarizing' if nivel > 0 else 'sem_ruido',
            nivel_ruido=nivel,
            max_iter=50,
            shots=512,
            seed=42
        )
        
        otimizador = OtimizadorQAOA(config)
        resultado = otimizador.otimizar(grafo)
        
        resultados.append({
            'nivel_ruido': nivel,
            'energia': resultado.energia_final,
            'tempo': resultado.tempo_execucao,
            'iteracoes': resultado.iteracoes
        })
        
        print(f"      Energia: {resultado.energia_final:.4f}\n")
    
    # Análise
    print("-"*80)
    print("RESULTADOS DA VARREDURA")
    print("-"*80 + "\n")
    
    df = pd.DataFrame(resultados)
    
    print(f"{'Nível Ruído':<15} {'Energia':<12} {'Tempo (s)':<12} {'Iterações':<12}")
    print("-"*80)
    for _, row in df.iterrows():
        print(f"{row['nivel_ruido']:<15.4f} {row['energia']:<12.4f} "
              f"{row['tempo']:<12.2f} {row['iteracoes']:<12.0f}")
    
    # Encontrar ótimo
    idx_melhor = df['energia'].idxmin()
    melhor = df.loc[idx_melhor]
    baseline = df.loc[0]  # Sem ruído
    
    print("\n" + "-"*80)
    print("ANÁLISE")
    print("-"*80 + "\n")
    
    print(f"Energia baseline (sem ruído): {baseline['energia']:.4f}")
    print(f"Melhor energia (com ruído):   {melhor['energia']:.4f}")
    print(f"Nível ótimo de ruído:         {melhor['nivel_ruido']:.4f}")
    
    melhoria = ((baseline['energia'] - melhor['energia']) / abs(baseline['energia'])) * 100
    print(f"Melhoria:                     {melhoria:+.2f}%\n")
    
    if melhoria > 0.5:
        print("✅ REGIÃO DE RUÍDO BENÉFICO IDENTIFICADA!")
        print(f"   O nível ótimo de ruído ({melhor['nivel_ruido']:.4f}) proporcionou")
        print(f"   uma melhoria de {melhoria:.2f}% em relação ao caso sem ruído.\n")
    else:
        print("⚠️  Ruído não foi benéfico neste experimento.")
        print("   Considere testar outros tipos de ruído ou configurações.\n")
    
    # Salvar resultados
    pasta = Path('resultados_exemplo_qaoa')
    pasta.mkdir(exist_ok=True)
    arquivo = pasta / 'varredura_niveis.csv'
    df.to_csv(arquivo, index=False)
    print(f"✓ Resultados salvos em: {arquivo}\n")
    
    return df


def exemplo_3_comparacao_tipos():
    """
    Exemplo 3: Comparação entre tipos de ruído.
    
    Compara depolarizing, amplitude damping e phase damping
    para identificar qual tipo é mais benéfico.
    """
    print("\n" + "="*80)
    print("EXEMPLO 3: COMPARAÇÃO ENTRE TIPOS DE RUÍDO")
    print("="*80 + "\n")
    
    # Configuração
    n_qubits = 20
    p_layers = 2
    nivel_fixo = 0.001
    tipos = ['sem_ruido', 'depolarizing', 'amplitude_damping', 'phase_damping']
    
    print(f"Configuração:")
    print(f"  • Qubits: {n_qubits}")
    print(f"  • P-layers: {p_layers}")
    print(f"  • Nível de ruído fixo: {nivel_fixo}")
    print(f"  • Tipos testados: {', '.join(tipos)}\n")
    
    # Criar grafo
    construtor = ConstrutorCircuitoQAOA(n_qubits=n_qubits, p_layers=p_layers, seed=42)
    grafo = construtor.criar_grafo_aleatorio(densidade=0.25)
    
    print("Executando comparação...\n")
    
    resultados = []
    
    for i, tipo in enumerate(tipos, 1):
        print(f"[{i}/{len(tipos)}] Testando {tipo}...")
        
        config = ConfigQAOA(
            n_qubits=n_qubits,
            p_layers=p_layers,
            tipo_ruido=tipo,
            nivel_ruido=nivel_fixo if tipo != 'sem_ruido' else 0.0,
            max_iter=50,
            shots=512,
            seed=42
        )
        
        otimizador = OtimizadorQAOA(config)
        resultado = otimizador.otimizar(grafo)
        
        resultados.append({
            'tipo_ruido': tipo,
            'energia': resultado.energia_final,
            'tempo': resultado.tempo_execucao,
            'iteracoes': resultado.iteracoes
        })
        
        print(f"      Energia: {resultado.energia_final:.4f}\n")
    
    # Análise
    print("-"*80)
    print("COMPARAÇÃO ENTRE TIPOS DE RUÍDO")
    print("-"*80 + "\n")
    
    df = pd.DataFrame(resultados)
    
    print(f"{'Tipo de Ruído':<25} {'Energia':<12} {'Tempo (s)':<12} {'Iterações':<12}")
    print("-"*80)
    for _, row in df.iterrows():
        print(f"{row['tipo_ruido']:<25} {row['energia']:<12.4f} "
              f"{row['tempo']:<12.2f} {row['iteracoes']:<12.0f}")
    
    # Análise comparativa
    baseline = df[df['tipo_ruido'] == 'sem_ruido'].iloc[0]
    df_com_ruido = df[df['tipo_ruido'] != 'sem_ruido'].copy()
    df_com_ruido['melhoria_pct'] = (
        (baseline['energia'] - df_com_ruido['energia']) / abs(baseline['energia']) * 100
    )
    
    print("\n" + "-"*80)
    print("ANÁLISE COMPARATIVA (relativo ao baseline sem ruído)")
    print("-"*80 + "\n")
    
    print(f"Energia baseline (sem ruído): {baseline['energia']:.4f}\n")
    
    print(f"{'Tipo de Ruído':<25} {'Energia':<12} {'Melhoria (%)':<15}")
    print("-"*80)
    for _, row in df_com_ruido.iterrows():
        simbolo = "✅" if row['melhoria_pct'] > 0.5 else ("⚠️" if row['melhoria_pct'] < -0.5 else "➖")
        print(f"{row['tipo_ruido']:<25} {row['energia']:<12.4f} "
              f"{row['melhoria_pct']:>+6.2f}%  {simbolo}")
    
    # Melhor tipo
    idx_melhor = df_com_ruido['energia'].idxmin()
    melhor = df_com_ruido.loc[idx_melhor]
    
    print("\n" + "-"*80)
    print("CONCLUSÃO")
    print("-"*80 + "\n")
    
    if melhor['melhoria_pct'] > 0.5:
        print(f"✅ MELHOR TIPO DE RUÍDO: {melhor['tipo_ruido']}")
        print(f"   Energia: {melhor['energia']:.4f}")
        print(f"   Melhoria: {melhor['melhoria_pct']:+.2f}% em relação ao baseline")
        print(f"\n   Este tipo de ruído foi o mais benéfico para este problema,")
        print(f"   sugerindo que suas características (taxa de decay, coerência, etc.)")
        print(f"   interagem favoravelmente com a dinâmica QAOA.\n")
    else:
        print("⚠️  Nenhum tipo de ruído foi significativamente benéfico.")
        print("   Considere ajustar o nível de ruído ou testar outros parâmetros.\n")
    
    # Salvar
    pasta = Path('resultados_exemplo_qaoa')
    pasta.mkdir(exist_ok=True)
    arquivo = pasta / 'comparacao_tipos.csv'
    df.to_csv(arquivo, index=False)
    print(f"✓ Resultados salvos em: {arquivo}\n")
    
    return df


def main():
    """Função principal - executa todos os exemplos."""
    print("\n" + "="*80)
    print(" "*20 + "EXEMPLOS PRÁTICOS: QAOA COM RUÍDO BENÉFICO")
    print("="*80 + "\n")
    
    print("Este script demonstra três experimentos fundamentais para")
    print("investigar o fenômeno de ruído benéfico em QAOA:\n")
    print("  1. Comparação básica: com vs. sem ruído")
    print("  2. Varredura de níveis de ruído")
    print("  3. Comparação entre tipos de ruído\n")
    
    escolha = input("Executar todos os exemplos? (s/N): ")
    
    if escolha.lower() == 's':
        # Executar todos
        exemplo_1_comparacao_basica()
        exemplo_2_varredura_niveis()
        exemplo_3_comparacao_tipos()
        
        print("\n" + "="*80)
        print("TODOS OS EXEMPLOS CONCLUÍDOS!")
        print("="*80 + "\n")
        print("Os resultados foram salvos em: resultados_exemplo_qaoa/\n")
    
    else:
        # Menu interativo
        while True:
            print("\nEscolha um exemplo para executar:\n")
            print("  1. Comparação básica (com vs. sem ruído)")
            print("  2. Varredura de níveis")
            print("  3. Comparação entre tipos")
            print("  4. Executar todos")
            print("\n  0. Sair\n")
            
            opcao = input("Opção: ").strip()
            
            if opcao == '1':
                exemplo_1_comparacao_basica()
            elif opcao == '2':
                exemplo_2_varredura_niveis()
            elif opcao == '3':
                exemplo_3_comparacao_tipos()
            elif opcao == '4':
                exemplo_1_comparacao_basica()
                exemplo_2_varredura_niveis()
                exemplo_3_comparacao_tipos()
                print("\n" + "="*80)
                print("TODOS OS EXEMPLOS CONCLUÍDOS!")
                print("="*80 + "\n")
                break
            elif opcao == '0':
                print("\nEncerrando...\n")
                break
            else:
                print("\n❌ Opção inválida. Tente novamente.\n")


if __name__ == "__main__":
    main()
