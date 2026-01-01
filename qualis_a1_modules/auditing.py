"""
Auditing and Traceability Module for Code→Method Mapping.

This module provides functions to generate explicit mappings between
scientific methods described in publications and their corresponding
code implementations, ensuring full auditability.

Purpose:
--------
Enable peer reviewers and readers to verify that published methods
are correctly implemented in code, satisfying Qualis A1 transparency requirements.

References:
-----------
Stodden, V. (2014). "The scientific method in practice: Reproducibility in the computational sciences."
    MIT Sloan Research Paper No. 4773-10.
"""

import os
import json
import csv
from datetime import datetime
from typing import List, Dict, Any, Optional
import logging

logger = logging.getLogger(__name__)


def gerar_tabela_codigo_metodo(
    pasta_resultados: str,
    formato: str = 'csv'
) -> str:
    """
    Gera tabela de rastreabilidade Código→Método para auditoria Qualis A1.
    
    A tabela mapeia explicitamente cada componente metodológico descrito
    no artigo científico para sua implementação em código, incluindo:
    - Arquivo, classe/função e número de linha
    - Parâmetros utilizados
    - Artefatos gerados
    - Referências bibliográficas
    
    Compliance com Padrões Qualis A1:
    ---------------------------------
    - Nature: "Code availability" section requirements
    - Quantum Journal: "Reproducibility statement"
    - Physical Review: "Supplementary material" standards
    
    Parameters:
    -----------
    pasta_resultados : str
        Diretório onde salvar a tabela
    formato : str, optional
        Formato de saída: 'csv', 'json' ou 'markdown' (padrão: 'csv')
    
    Returns:
    --------
    str
        Caminho do arquivo gerado
    
    Table Format:
    ------------
    | Componente Metodológico | Arquivo | Classe/Função | Linha | Parâmetros | Artefatos | Referência |
    |------------------------|---------|---------------|-------|------------|-----------|------------|
    
    Examples:
    ---------
    >>> tabela_path = gerar_tabela_codigo_metodo('./results', formato='csv')
    >>> # Tabela CSV gerada com mapeamento completo
    
    References:
    -----------
    Peng, R. D. (2011). "Reproducible research in computational science." Science, 334(6060), 1226-1227.
    """
    os.makedirs(pasta_resultados, exist_ok=True)
    
    # Mapeamento completo Método → Código
    mapeamento = [
        {
            'componente': 'Canal de Depolarização',
            'secao_artigo': 'Seção 2.1: Modelos de Ruído',
            'arquivo': 'framework_investigativo_completo.py',
            'classe_funcao': 'RuidoDepolarizante',
            'linha_inicio': 1523,
            'linha_fim': 1548,
            'parametros': 'nivel: float (probabilidade de erro p)',
            'equacao': 'ℰ(ρ) = (1-p)ρ + p/3(XρX + YρY + ZρZ)',
            'kraus_operators': 'K₀=√(1-p)I, K₁=√(p/3)X, K₂=√(p/3)Y, K₃=√(p/3)Z',
            'artefatos': 'Canal quântico aplicado via PennyLane',
            'referencia': 'Nielsen & Chuang (2010), Chapter 8',
            'validacao': 'tests/test_modelo_ruido.py::TestRuidoDepolarizante'
        },
        {
            'componente': 'Canal de Amplitude Damping',
            'secao_artigo': 'Seção 2.1: Modelos de Ruído',
            'arquivo': 'framework_investigativo_completo.py',
            'classe_funcao': 'RuidoAmplitudeDamping',
            'linha_inicio': 1551,
            'linha_fim': 1574,
            'parametros': 'nivel: float (taxa de decaimento γ)',
            'equacao': 'Relaxamento T1: |1⟩ → |0⟩',
            'kraus_operators': 'K₀=[[1,0],[0,√(1-γ)]], K₁=[[0,√γ],[0,0]]',
            'artefatos': 'Canal quântico aplicado via PennyLane',
            'referencia': 'Clerk et al. (2010), Rev. Mod. Phys.',
            'validacao': 'tests/test_modelo_ruido.py::TestRuidoAmplitudeDamping'
        },
        {
            'componente': 'Canal de Phase Damping',
            'secao_artigo': 'Seção 2.1: Modelos de Ruído',
            'arquivo': 'framework_investigativo_completo.py',
            'classe_funcao': 'RuidoPhaseDamping',
            'linha_inicio': 1577,
            'linha_fim': 1600,
            'parametros': 'nivel: float (taxa de perda de fase λ)',
            'equacao': 'Decoerência T2 (perda de fase)',
            'kraus_operators': 'K₀=[[1,0],[0,√(1-λ)]], K₁=[[0,0],[0,√λ]]',
            'artefatos': 'Canal quântico aplicado via PennyLane',
            'referencia': 'Schlosshauer (2007), Decoherence book',
            'validacao': 'tests/test_modelo_ruido.py::TestRuidoPhaseDamping'
        },
        {
            'componente': 'Canal de Crosstalk',
            'secao_artigo': 'Seção 2.1: Modelos de Ruído',
            'arquivo': 'framework_investigativo_completo.py',
            'classe_funcao': 'RuidoCrosstalk',
            'linha_inicio': 1604,
            'linha_fim': 1629,
            'parametros': 'nivel: float (intensidade de crosstalk)',
            'equacao': 'Erro correlacionado entre qubits vizinhos',
            'kraus_operators': 'Combinação de canais de depolarização + CNOT',
            'artefatos': 'Canal quântico aplicado via PennyLane',
            'referencia': 'Kandala et al. (2019), Nature',
            'validacao': 'tests/test_modelo_ruido.py::TestRuidoCrosstalk'
        },
        {
            'componente': 'Quantum Natural Gradient (QNG)',
            'secao_artigo': 'Seção 2.2: Otimizadores',
            'arquivo': 'framework_investigativo_completo.py',
            'classe_funcao': 'QNG',
            'linha_inicio': 904,
            'linha_fim': 919,
            'parametros': 'taxa_aprendizado: float, reg: float',
            'equacao': 'θ_{t+1} = θ_t - η g⁻¹(θ_t) ∇_θ L(θ_t)',
            'kraus_operators': 'N/A',
            'artefatos': 'Gradiente natural usando métrica de Fubini-Study',
            'referencia': 'Stokes et al. (2020), Quantum, 4, 269',
            'validacao': 'Validação via convergência em experimentos'
        },
        {
            'componente': 'Effect Size: Cohen\'s d',
            'secao_artigo': 'Seção 3.2: Análise Estatística',
            'arquivo': 'framework_investigativo_completo.py',
            'classe_funcao': 'TestesEstatisticosAvancados.cohen_d',
            'linha_inicio': 987,
            'linha_fim': 1012,
            'parametros': 'grupo1: array, grupo2: array',
            'equacao': 'd = (μ₁ - μ₂) / σ_pooled',
            'kraus_operators': 'N/A',
            'artefatos': 'Medida de tamanho de efeito',
            'referencia': 'Cohen (1988), Statistical Power Analysis',
            'validacao': 'Interpretação: |d|<0.2 pequeno, 0.2≤|d|<0.5 médio, |d|≥0.8 grande'
        },
        {
            'componente': 'Validação de Operadores de Kraus',
            'secao_artigo': 'Seção 2.1: Validação Matemática',
            'arquivo': 'qualis_a1_modules/validation.py',
            'classe_funcao': 'validar_operadores_kraus',
            'linha_inicio': 24,
            'linha_fim': 128,
            'parametros': 'operadores: List[ndarray], tol: float',
            'equacao': 'Σ_k K_k† K_k = I (condição de completude)',
            'kraus_operators': 'Validação de qualquer conjunto de Kraus',
            'artefatos': 'Verificação de preservação de traço',
            'referencia': 'Nielsen & Chuang (2010), Section 8.2.3',
            'validacao': 'qualis_a1_modules/validation.py (self-test)'
        },
        {
            'componente': 'Configuração de Seeds Reprodutíveis',
            'secao_artigo': 'Seção 3.1: Metodologia',
            'arquivo': 'qualis_a1_modules/reproducibility.py',
            'classe_funcao': 'configurar_seeds_reprodutiveis',
            'linha_inicio': 35,
            'linha_fim': 145,
            'parametros': 'seed: int',
            'equacao': 'N/A',
            'kraus_operators': 'N/A',
            'artefatos': 'Seeds para NumPy, random, PennyLane, Optuna, Qiskit',
            'referencia': 'Wilson et al. (2014), PLoS Biol',
            'validacao': 'Manifesto de execução documenta seeds'
        },
        {
            'componente': 'Correção de Bonferroni',
            'secao_artigo': 'Seção 3.2: Testes Post-Hoc',
            'arquivo': 'qualis_a1_modules/statistical_extensions.py',
            'classe_funcao': 'testes_post_hoc_com_correcao',
            'linha_inicio': 'TBD',
            'linha_fim': 'TBD',
            'parametros': 'df: DataFrame, metodo_correcao: str',
            'equacao': 'α_ajustado = α / m (m = número de comparações)',
            'kraus_operators': 'N/A',
            'artefatos': 'P-valores ajustados para múltiplas comparações',
            'referencia': 'Dunn (1961), JASA',
            'validacao': 'Testes estatísticos nos resultados'
        },
        {
            'componente': 'Análise de Poder Estatístico',
            'secao_artigo': 'Seção 3.2: Validação Estatística',
            'arquivo': 'qualis_a1_modules/statistical_extensions.py',
            'classe_funcao': 'analise_poder_estatistico',
            'linha_inicio': 'TBD',
            'linha_fim': 'TBD',
            'parametros': 'effect_size: float, n_per_group: int, alpha: float',
            'equacao': 'Poder = 1 - β = P(rejeitar H₀ | H₀ falsa)',
            'kraus_operators': 'N/A',
            'artefatos': 'Cálculo de poder (1-β)',
            'referencia': 'Cohen (1988), Statistical Power Analysis',
            'validacao': 'Reportado em resultados experimentais'
        },
        {
            'componente': 'Thermal Relaxation Noise',
            'secao_artigo': 'Seção 2.1: Modelos de Ruído',
            'arquivo': 'framework_investigativo_completo.py',
            'classe_funcao': 'RuidoThermal',
            'linha_inicio': 1444,
            'linha_fim': 1457,
            'parametros': 'nivel: float',
            'equacao': 'Combinação de Amplitude + Phase Damping (T1/T2)',
            'kraus_operators': 'Composição de K_amplitude e K_phase',
            'artefatos': 'Modelagem de relaxamento térmico',
            'referencia': 'Hardware typical T1~50μs, T2~20μs',
            'validacao': 'tests/test_modelo_ruido.py::TestRuidoThermal'
        },
        {
            'componente': 'Bit-Flip Noise',
            'secao_artigo': 'Seção 2.1: Modelos de Ruído',
            'arquivo': 'framework_investigativo_completo.py',
            'classe_funcao': 'RuidoBitFlip',
            'linha_inicio': 1459,
            'linha_fim': 1471,
            'parametros': 'nivel: float (probabilidade p)',
            'equacao': 'Aplicação de porta X com probabilidade p',
            'kraus_operators': 'K₀=√(1-p)I, K₁=√p X',
            'artefatos': 'Canal de bit-flip via PennyLane',
            'referencia': 'Nielsen & Chuang (2010)',
            'validacao': 'tests/test_modelo_ruido.py::TestRuidoBitFlip'
        },
        {
            'componente': 'Phase-Flip Noise',
            'secao_artigo': 'Seção 2.1: Modelos de Ruído',
            'arquivo': 'framework_investigativo_completo.py',
            'classe_funcao': 'RuidoPhaseFlip',
            'linha_inicio': 1473,
            'linha_fim': 1485,
            'parametros': 'nivel: float (probabilidade p)',
            'equacao': 'Aplicação de porta Z com probabilidade p',
            'kraus_operators': 'K₀=√(1-p)I, K₁=√p Z',
            'artefatos': 'Canal de phase-flip via PennyLane',
            'referencia': 'Nielsen & Chuang (2010)',
            'validacao': 'tests/test_modelo_ruido.py::TestRuidoPhaseFlip'
        },
        {
            'componente': 'Pink Noise (1/f)',
            'secao_artigo': 'Seção 2.1: Modelos de Ruído',
            'arquivo': 'framework_investigativo_completo.py',
            'classe_funcao': 'RuidoPinkNoise',
            'linha_inicio': 1487,
            'linha_fim': 1503,
            'parametros': 'nivel: float (escala Gaussiana)',
            'equacao': 'Phase damping com intensidade ~ N(0, σ²)',
            'kraus_operators': 'Variação estocástica por qubit',
            'artefatos': 'Modelagem de ruído de baixa frequência',
            'referencia': 'Característica de hardware superconductor',
            'validacao': 'tests/test_modelo_ruido.py::TestRuidoPinkNoise'
        },
        {
            'componente': 'Readout Error',
            'secao_artigo': 'Seção 2.1: Modelos de Ruído',
            'arquivo': 'framework_investigativo_completo.py',
            'classe_funcao': 'RuidoReadoutError',
            'linha_inicio': 1505,
            'linha_fim': 1520,
            'parametros': 'nivel: float (probabilidade de erro)',
            'equacao': 'Aproximação via bit-flip pós-medição',
            'kraus_operators': 'Similar a bit-flip',
            'artefatos': 'Modelagem de erros de medição',
            'referencia': 'Típico em hardware: 1-5% erro',
            'validacao': 'tests/test_modelo_ruido.py::TestRuidoReadoutError'
        },
        {
            'componente': 'Correlated Noise',
            'secao_artigo': 'Seção 2.1: Modelos de Ruído',
            'arquivo': 'framework_investigativo_completo.py',
            'classe_funcao': 'RuidoCorrelacionado',
            'linha_inicio': 1632,
            'linha_fim': 1654,
            'parametros': 'nivel: float',
            'equacao': 'Phase damping global (todos os qubits)',
            'kraus_operators': 'Aplicação coletiva de K_phase',
            'artefatos': 'Modelagem de ruído correlacionado',
            'referencia': 'Greenbaum (2015), Quantum Gate Set Tomography',
            'validacao': 'tests/test_modelo_ruido.py::TestRuidoCorrelacionado'
        },
    ]
    
    # Gerar arquivo no formato especificado
    if formato == 'csv':
        output_path = os.path.join(pasta_resultados, 'tabela_codigo_metodo.csv')
        with open(output_path, 'w', newline='', encoding='utf-8') as f:
            fieldnames = [
                'Componente Metodológico', 'Seção do Artigo', 'Arquivo',
                'Classe/Função', 'Linha Início', 'Linha Fim', 'Parâmetros',
                'Equação', 'Operadores de Kraus', 'Artefatos',
                'Referência', 'Validação'
            ]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for item in mapeamento:
                writer.writerow({
                    'Componente Metodológico': item['componente'],
                    'Seção do Artigo': item['secao_artigo'],
                    'Arquivo': item['arquivo'],
                    'Classe/Função': item['classe_funcao'],
                    'Linha Início': item['linha_inicio'],
                    'Linha Fim': item['linha_fim'],
                    'Parâmetros': item['parametros'],
                    'Equação': item['equacao'],
                    'Operadores de Kraus': item['kraus_operators'],
                    'Artefatos': item['artefatos'],
                    'Referência': item['referencia'],
                    'Validação': item['validacao']
                })
    
    elif formato == 'json':
        output_path = os.path.join(pasta_resultados, 'tabela_codigo_metodo.json')
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump({
                'metadata': {
                    'generated_at': datetime.now().isoformat(),
                    'purpose': 'Code-to-Method Traceability for Qualis A1 Compliance',
                    'total_mappings': len(mapeamento)
                },
                'mappings': mapeamento
            }, f, indent=2, ensure_ascii=False)
    
    elif formato == 'markdown':
        output_path = os.path.join(pasta_resultados, 'tabela_codigo_metodo.md')
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write("# Tabela de Rastreabilidade Código → Método\n\n")
            f.write("## Beneficial Quantum Noise in Variational Quantum Classifiers\n\n")
            f.write(f"**Gerado em:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write("---\n\n")
            
            for i, item in enumerate(mapeamento, 1):
                f.write(f"### {i}. {item['componente']}\n\n")
                f.write(f"- **Seção do Artigo:** {item['secao_artigo']}\n")
                f.write(f"- **Implementação:** `{item['arquivo']}::{item['classe_funcao']}` (linhas {item['linha_inicio']}-{item['linha_fim']})\n")
                f.write(f"- **Equação:** {item['equacao']}\n")
                if item['kraus_operators'] != 'N/A':
                    f.write(f"- **Operadores de Kraus:** {item['kraus_operators']}\n")
                f.write(f"- **Referência:** {item['referencia']}\n")
                f.write(f"- **Validação:** `{item['validacao']}`\n\n")
                f.write("---\n\n")
    
    else:
        raise ValueError(f"Formato não suportado: {formato}. Use 'csv', 'json' ou 'markdown'.")
    
    logger.info(f"✓ Tabela código→método gerada: {output_path}")
    logger.info(f"  Total de mapeamentos: {len(mapeamento)}")
    
    return output_path


if __name__ == "__main__":
    # Teste do módulo
    logging.basicConfig(level=logging.INFO)
    
    print("\n" + "=" * 80)
    print("GERAÇÃO DE TABELA DE RASTREABILIDADE CÓDIGO→MÉTODO")
    print("=" * 80 + "\n")
    
    # Gerar em todos os formatos
    for fmt in ['csv', 'json', 'markdown']:
        print(f"Gerando formato: {fmt}")
        path = gerar_tabela_codigo_metodo('/tmp/test_auditing', formato=fmt)
        print(f"  ✓ Arquivo: {path}\n")
    
    print("=" * 80)
    print("✓ TODOS OS FORMATOS GERADOS COM SUCESSO")
    print("=" * 80)
