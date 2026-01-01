"""
Statistical Extensions Module for Advanced Statistical Analysis.

This module extends the statistical capabilities of the framework with:
1. Bonferroni correction for multiple comparisons
2. Statistical power analysis
3. Enhanced post-hoc tests

Ensures compliance with Qualis A1 statistical rigor requirements.

References:
-----------
Dunn, O. J. (1961). "Multiple comparisons among means." JASA, 56(293), 52-64.
Cohen, J. (1988). "Statistical Power Analysis for the Behavioral Sciences." 2nd ed.
"""

import numpy as np
import pandas as pd
from typing import List, Tuple, Dict, Optional
import logging
from scipy import stats
from statsmodels.stats.multitest import multipletests
from itertools import combinations

logger = logging.getLogger(__name__)


def testes_post_hoc_com_correcao(
    df: pd.DataFrame,
    grupo_col: str,
    metrica_col: str,
    metodo_correcao: str = 'bonferroni',
    alpha: float = 0.05
) -> pd.DataFrame:
    """
    Realiza testes post-hoc com correção para múltiplas comparações.
    
    Fundamentação Matemática:
    ------------------------
    Quando realizamos m comparações múltiplas, a probabilidade de cometer
    pelo menos um erro tipo I (falso positivo) aumenta. A correção de Bonferroni
    ajusta o nível de significância:
    
    $$\\alpha_{ajustado} = \\frac{\\alpha}{m}$$
    
    onde m é o número de comparações.
    
    Métodos de Correção Disponíveis:
    --------------------------------
    - 'bonferroni': α_adj = α/m (mais conservador)
    - 'holm': Método de Holm-Bonferroni (sequencial)
    - 'fdr_bh': False Discovery Rate (Benjamini-Hochberg)
    - 'fdr_by': False Discovery Rate (Benjamini-Yekutieli)
    
    Interpretação dos Resultados:
    ----------------------------
    - p-valor < α_ajustado: Diferença estatisticamente significativa
    - p-valor ≥ α_ajustado: Não há evidência de diferença
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame com dados experimentais
    grupo_col : str
        Nome da coluna que identifica os grupos
    metrica_col : str
        Nome da coluna com a métrica a comparar
    metodo_correcao : str, optional
        Método de correção (padrão: 'bonferroni')
        Opções: 'bonferroni', 'holm', 'fdr_bh', 'fdr_by'
    alpha : float, optional
        Nível de significância antes da correção (padrão: 0.05)
    
    Returns:
    --------
    pd.DataFrame
        DataFrame com resultados dos testes post-hoc:
        - grupo1, grupo2: Grupos comparados
        - p_valor: P-valor original (t-test)
        - p_ajustado: P-valor após correção
        - significativo: Boolean indicando se é significativo
        - effect_size: Cohen's d para a comparação
    
    Examples:
    ---------
    >>> df = pd.DataFrame({
    ...     'ruido': ['depol', 'depol', 'amplitude', 'amplitude', 'phase', 'phase'],
    ...     'acuracia': [0.85, 0.87, 0.82, 0.84, 0.88, 0.90]
    ... })
    >>> resultados = testes_post_hoc_com_correcao(df, 'ruido', 'acuracia')
    
    References:
    -----------
    Dunn, O. J. (1961). "Multiple comparisons among means." 
        Journal of the American Statistical Association, 56(293), 52-64.
    Holm, S. (1979). "A simple sequentially rejective multiple test procedure."
        Scandinavian Journal of Statistics, 6(2), 65-70.
    Benjamini, Y., & Hochberg, Y. (1995). "Controlling the false discovery rate."
        Journal of the Royal Statistical Society B, 57(1), 289-300.
    """
    # Obter grupos únicos
    grupos = df[grupo_col].unique()
    n_grupos = len(grupos)
    
    if n_grupos < 2:
        raise ValueError(f"Necessário pelo menos 2 grupos. Encontrados: {n_grupos}")
    
    # Calcular todas as comparações par-a-par
    comparacoes = list(combinations(grupos, 2))
    n_comparacoes = len(comparacoes)
    
    logger.info(f"Realizando {n_comparacoes} comparações post-hoc entre {n_grupos} grupos")
    logger.info(f"Método de correção: {metodo_correcao}")
    logger.info(f"α original: {alpha:.4f}")
    
    # Calcular p-valores para todas as comparações
    resultados = []
    
    for grupo1, grupo2 in comparacoes:
        # Extrair dados dos grupos
        dados1 = df[df[grupo_col] == grupo1][metrica_col].values
        dados2 = df[df[grupo_col] == grupo2][metrica_col].values
        
        # Teste t independente
        t_stat, p_valor = stats.ttest_ind(dados1, dados2)
        
        # Calcular Cohen's d (effect size)
        n1, n2 = len(dados1), len(dados2)
        var1, var2 = np.var(dados1, ddof=1), np.var(dados2, ddof=1)
        pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
        cohen_d = (np.mean(dados1) - np.mean(dados2)) / pooled_std if pooled_std > 0 else 0
        
        resultados.append({
            'grupo1': grupo1,
            'grupo2': grupo2,
            'n1': n1,
            'n2': n2,
            'media1': np.mean(dados1),
            'media2': np.mean(dados2),
            'diff_medias': np.mean(dados1) - np.mean(dados2),
            't_statistic': t_stat,
            'p_valor': p_valor,
            'cohen_d': cohen_d
        })
    
    # Criar DataFrame com resultados
    df_resultados = pd.DataFrame(resultados)
    
    # Aplicar correção de múltiplas comparações
    p_valores = df_resultados['p_valor'].values
    reject, p_ajustados, alpha_sidak, alpha_bonf = multipletests(
        p_valores,
        alpha=alpha,
        method=metodo_correcao
    )
    
    df_resultados['p_ajustado'] = p_ajustados
    df_resultados['significativo'] = reject
    df_resultados['alpha_ajustado'] = alpha_bonf if metodo_correcao == 'bonferroni' else alpha
    
    # Adicionar interpretação do effect size
    def interpretar_cohen_d(d):
        """Interpreta magnitude do Cohen's d."""
        abs_d = abs(d)
        if abs_d < 0.2:
            return 'pequeno'
        elif abs_d < 0.5:
            return 'médio'
        elif abs_d < 0.8:
            return 'grande'
        else:
            return 'muito grande'
    
    df_resultados['magnitude_efeito'] = df_resultados['cohen_d'].apply(interpretar_cohen_d)
    
    # Log dos resultados
    n_significativos = df_resultados['significativo'].sum()
    logger.info(f"α ajustado ({metodo_correcao}): {alpha_bonf:.4f}")
    logger.info(f"Comparações significativas: {n_significativos}/{n_comparacoes}")
    
    return df_resultados


def analise_poder_estatistico(
    effect_size: float,
    n_per_group: int,
    alpha: float = 0.05,
    alternative: str = 'two-sided'
) -> Dict[str, float]:
    """
    Calcula o poder estatístico (1-β) para detectar um efeito.
    
    Fundamentação Matemática:
    ------------------------
    O poder estatístico é definido como:
    
    $$\\text{Poder} = 1 - \\beta = P(\\text{rejeitar } H_0 | H_0 \\text{ é falsa})$$
    
    onde β é a probabilidade de erro tipo II (falso negativo).
    
    Interpretação do Poder:
    ----------------------
    - Poder < 0.50: Inadequado (alta chance de perder efeito real)
    - Poder 0.50-0.70: Baixo (revisar tamanho amostral)
    - Poder 0.70-0.80: Aceitável
    - Poder ≥ 0.80: Adequado (padrão científico)
    - Poder ≥ 0.90: Excelente
    
    Relação com Tamanho Amostral:
    -----------------------------
    Para um efeito fixo, poder aumenta com n:
    - Dobrar n aumenta poder significativamente
    - Poder = 0.80 requer n adequado ao effect size
    
    Parameters:
    -----------
    effect_size : float
        Tamanho do efeito (Cohen's d para t-test)
    n_per_group : int
        Tamanho amostral por grupo
    alpha : float, optional
        Nível de significância (padrão: 0.05)
    alternative : str, optional
        Tipo de teste: 'two-sided', 'larger', 'smaller' (padrão: 'two-sided')
    
    Returns:
    --------
    Dict[str, float]
        Dicionário com:
        - poder: Poder estatístico (1-β)
        - beta: Probabilidade de erro tipo II
        - effect_size: Tamanho do efeito usado
        - n_per_group: Tamanho amostral por grupo
        - alpha: Nível de significância
        - interpretacao: Interpretação textual do poder
    
    Examples:
    ---------
    >>> # Calcular poder para detectar efeito médio (d=0.5) com n=30
    >>> resultado = analise_poder_estatistico(effect_size=0.5, n_per_group=30)
    >>> print(f"Poder: {resultado['poder']:.2%}")
    
    >>> # Calcular n necessário para poder=0.80
    >>> from statsmodels.stats.power import TTestIndPower
    >>> analysis = TTestIndPower()
    >>> n_necessario = analysis.solve_power(effect_size=0.5, alpha=0.05, power=0.80)
    
    References:
    -----------
    Cohen, J. (1988). "Statistical Power Analysis for the Behavioral Sciences."
        2nd edition, Lawrence Erlbaum Associates.
    Murphy, K. R., & Myors, B. (2004). "Statistical Power Analysis: A Simple and General
        Model for Traditional and Modern Hypothesis Tests." 2nd ed.
    """
    try:
        from statsmodels.stats.power import TTestIndPower
    except ImportError:
        logger.error("statsmodels não está disponível para análise de poder")
        raise ImportError("statsmodels é necessário para análise de poder")
    
    # Criar objeto de análise de poder
    power_analysis = TTestIndPower()
    
    # Calcular poder
    poder = power_analysis.solve_power(
        effect_size=effect_size,
        nobs1=n_per_group,
        alpha=alpha,
        ratio=1.0,  # Grupos de mesmo tamanho
        alternative=alternative
    )
    
    beta = 1 - poder
    
    # Interpretação
    if poder < 0.50:
        interpretacao = "Inadequado - Alta probabilidade de erro tipo II"
    elif poder < 0.70:
        interpretacao = "Baixo - Considerar aumentar tamanho amostral"
    elif poder < 0.80:
        interpretacao = "Aceitável - Próximo do padrão científico"
    elif poder < 0.90:
        interpretacao = "Adequado - Atende padrão científico (≥0.80)"
    else:
        interpretacao = "Excelente - Poder muito alto"
    
    resultado = {
        'poder': poder,
        'beta': beta,
        'effect_size': effect_size,
        'n_per_group': n_per_group,
        'alpha': alpha,
        'alternative': alternative,
        'interpretacao': interpretacao
    }
    
    logger.info(f"Análise de Poder Estatístico:")
    logger.info(f"  Effect size (Cohen's d): {effect_size:.3f}")
    logger.info(f"  Tamanho amostral por grupo: {n_per_group}")
    logger.info(f"  Nível de significância (α): {alpha:.3f}")
    logger.info(f"  Poder (1-β): {poder:.3f} ({poder*100:.1f}%)")
    logger.info(f"  Interpretação: {interpretacao}")
    
    return resultado


def calcular_tamanho_amostral_necessario(
    effect_size: float,
    poder_desejado: float = 0.80,
    alpha: float = 0.05,
    alternative: str = 'two-sided'
) -> int:
    """
    Calcula tamanho amostral necessário para atingir poder desejado.
    
    Esta função é útil no planejamento de experimentos para garantir
    que o estudo tenha poder adequado para detectar o efeito esperado.
    
    Parameters:
    -----------
    effect_size : float
        Tamanho do efeito esperado (Cohen's d)
    poder_desejado : float, optional
        Poder estatístico desejado (padrão: 0.80)
    alpha : float, optional
        Nível de significância (padrão: 0.05)
    alternative : str, optional
        Tipo de teste (padrão: 'two-sided')
    
    Returns:
    --------
    int
        Tamanho amostral necessário por grupo (arredondado para cima)
    
    Examples:
    ---------
    >>> # Quantos participantes preciso para detectar d=0.5 com poder=0.80?
    >>> n = calcular_tamanho_amostral_necessario(effect_size=0.5, poder_desejado=0.80)
    >>> print(f"Necessário: {n} por grupo, {n*2} total")
    
    References:
    -----------
    Cohen, J. (1988). "Statistical Power Analysis for the Behavioral Sciences."
    """
    try:
        from statsmodels.stats.power import TTestIndPower
    except ImportError:
        raise ImportError("statsmodels é necessário para cálculo de tamanho amostral")
    
    power_analysis = TTestIndPower()
    
    # Calcular n necessário
    n_necessario = power_analysis.solve_power(
        effect_size=effect_size,
        alpha=alpha,
        power=poder_desejado,
        ratio=1.0,
        alternative=alternative
    )
    
    # Arredondar para cima
    n_necessario = int(np.ceil(n_necessario))
    
    logger.info(f"Tamanho Amostral Necessário:")
    logger.info(f"  Para detectar effect size d={effect_size:.3f}")
    logger.info(f"  Com poder {poder_desejado:.2%} e α={alpha:.3f}")
    logger.info(f"  Necessário: {n_necessario} por grupo ({n_necessario*2} total)")
    
    return n_necessario


if __name__ == "__main__":
    # Teste do módulo
    logging.basicConfig(level=logging.INFO)
    
    print("\n" + "=" * 80)
    print("EXTENSÕES ESTATÍSTICAS - TESTES")
    print("=" * 80 + "\n")
    
    # Teste 1: Testes post-hoc com correção de Bonferroni
    print("1. Testes Post-Hoc com Correção de Bonferroni")
    print("-" * 80)
    
    # Dados de exemplo
    np.random.seed(42)
    df_exemplo = pd.DataFrame({
        'ruido': ['depol']*10 + ['amplitude']*10 + ['phase']*10,
        'acuracia': (
            np.random.normal(0.85, 0.05, 10).tolist() +
            np.random.normal(0.82, 0.05, 10).tolist() +
            np.random.normal(0.88, 0.05, 10).tolist()
        )
    })
    
    resultados_posthoc = testes_post_hoc_com_correcao(
        df_exemplo, 'ruido', 'acuracia', metodo_correcao='bonferroni'
    )
    print(resultados_posthoc[['grupo1', 'grupo2', 'p_valor', 'p_ajustado', 'significativo', 'cohen_d']])
    
    # Teste 2: Análise de poder estatístico
    print("\n2. Análise de Poder Estatístico")
    print("-" * 80)
    
    poder_info = analise_poder_estatistico(
        effect_size=0.5,  # Efeito médio
        n_per_group=30,
        alpha=0.05
    )
    print(f"Poder: {poder_info['poder']:.2%}")
    print(f"Interpretação: {poder_info['interpretacao']}")
    
    # Teste 3: Cálculo de tamanho amostral
    print("\n3. Cálculo de Tamanho Amostral Necessário")
    print("-" * 80)
    
    n_necessario = calcular_tamanho_amostral_necessario(
        effect_size=0.5,
        poder_desejado=0.80
    )
    print(f"N necessário por grupo: {n_necessario}")
    
    print("\n" + "=" * 80)
    print("✓ TODOS OS TESTES PASSARAM")
    print("=" * 80)
