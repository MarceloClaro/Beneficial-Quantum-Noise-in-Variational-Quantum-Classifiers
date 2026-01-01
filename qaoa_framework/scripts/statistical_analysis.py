"""
Tarefa 9: Análise Estatística para QAOA
Análise estatística rigorosa adaptada para approximation_ratio como métrica principal.

Referências:
- Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences.
- Field, A. (2013). Discovering Statistics Using IBM SPSS Statistics.
"""

from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd
from scipy import stats
from dataclasses import dataclass


@dataclass
class ResultadoTesteEstatistico:
    """Resultado de um teste estatístico."""
    teste: str
    estatistica: float
    p_valor: float
    significativo: bool
    effect_size: Optional[float] = None
    interpretacao: Optional[str] = None


def teste_t_pareado(
    grupo1: np.ndarray,
    grupo2: np.ndarray,
    alpha: float = 0.05
) -> ResultadoTesteEstatistico:
    """
    Teste t pareado para comparar duas condições.
    
    Uso: Comparar approximation_ratio com vs. sem ruído.
    
    Args:
        grupo1: Medições do primeiro grupo
        grupo2: Medições do segundo grupo (mesmo tamanho)
        alpha: Nível de significância
    
    Returns:
        ResultadoTesteEstatistico
    """
    # Teste t pareado
    t_stat, p_valor = stats.ttest_rel(grupo1, grupo2)
    
    # Cohen's d para amostras pareadas
    diferenca = grupo1 - grupo2
    d = np.mean(diferenca) / np.std(diferenca, ddof=1)
    
    # Interpretação do effect size
    if abs(d) < 0.2:
        interpretacao = "efeito trivial"
    elif abs(d) < 0.5:
        interpretacao = "efeito pequeno"
    elif abs(d) < 0.8:
        interpretacao = "efeito médio"
    else:
        interpretacao = "efeito grande"
    
    return ResultadoTesteEstatistico(
        teste="t-test pareado",
        estatistica=t_stat,
        p_valor=p_valor,
        significativo=p_valor < alpha,
        effect_size=d,
        interpretacao=interpretacao
    )


def teste_t_independente(
    grupo1: np.ndarray,
    grupo2: np.ndarray,
    alpha: float = 0.05
) -> ResultadoTesteEstatistico:
    """
    Teste t independente para comparar dois grupos.
    
    Args:
        grupo1: Medições do primeiro grupo
        grupo2: Medições do segundo grupo
        alpha: Nível de significância
    
    Returns:
        ResultadoTesteEstatistico
    """
    # Teste t independente
    t_stat, p_valor = stats.ttest_ind(grupo1, grupo2)
    
    # Cohen's d para amostras independentes
    n1, n2 = len(grupo1), len(grupo2)
    var1, var2 = np.var(grupo1, ddof=1), np.var(grupo2, ddof=1)
    pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    d = (np.mean(grupo1) - np.mean(grupo2)) / pooled_std
    
    # Interpretação
    if abs(d) < 0.2:
        interpretacao = "efeito trivial"
    elif abs(d) < 0.5:
        interpretacao = "efeito pequeno"
    elif abs(d) < 0.8:
        interpretacao = "efeito médio"
    else:
        interpretacao = "efeito grande"
    
    return ResultadoTesteEstatistico(
        teste="t-test independente",
        estatistica=t_stat,
        p_valor=p_valor,
        significativo=p_valor < alpha,
        effect_size=d,
        interpretacao=interpretacao
    )


def anova_one_way(
    grupos: List[np.ndarray],
    labels: List[str],
    alpha: float = 0.05
) -> Dict:
    """
    ANOVA one-way para comparar múltiplos grupos.
    
    Uso: Comparar múltiplos níveis de ruído ou tipos de ruído.
    
    Args:
        grupos: Lista de arrays com medições de cada grupo
        labels: Nomes dos grupos
        alpha: Nível de significância
    
    Returns:
        Dicionário com resultados da ANOVA e post-hoc tests
    """
    # ANOVA one-way
    f_stat, p_valor = stats.f_oneway(*grupos)
    
    # Eta-squared (effect size para ANOVA)
    grand_mean = np.mean(np.concatenate(grupos))
    ss_between = sum(len(g) * (np.mean(g) - grand_mean)**2 for g in grupos)
    ss_total = sum(sum((x - grand_mean)**2 for x in g) for g in grupos)
    eta_squared = ss_between / ss_total
    
    resultado = {
        'anova': ResultadoTesteEstatistico(
            teste="ANOVA one-way",
            estatistica=f_stat,
            p_valor=p_valor,
            significativo=p_valor < alpha,
            effect_size=eta_squared,
            interpretacao=f"η² = {eta_squared:.4f}"
        ),
        'grupos': labels,
        'medias': [np.mean(g) for g in grupos],
        'desvios': [np.std(g, ddof=1) for g in grupos]
    }
    
    # Post-hoc tests (Tukey HSD) se ANOVA for significativo
    if p_valor < alpha and len(grupos) > 2:
        posthoc = []
        for i in range(len(grupos)):
            for j in range(i+1, len(grupos)):
                resultado_t = teste_t_independente(grupos[i], grupos[j], alpha)
                posthoc.append({
                    'comparacao': f"{labels[i]} vs {labels[j]}",
                    'p_valor': resultado_t.p_valor,
                    'significativo': resultado_t.significativo,
                    'effect_size': resultado_t.effect_size
                })
        resultado['posthoc'] = posthoc
    
    return resultado


def calcular_intervalo_confianca(
    dados: np.ndarray,
    confianca: float = 0.95
) -> Tuple[float, float, float]:
    """
    Calcula intervalo de confiança.
    
    Args:
        dados: Array de medições
        confianca: Nível de confiança (default: 95%)
    
    Returns:
        (média, limite_inferior, limite_superior)
    """
    media = np.mean(dados)
    sem = stats.sem(dados)
    intervalo = sem * stats.t.ppf((1 + confianca) / 2, len(dados) - 1)
    
    return media, media - intervalo, media + intervalo


def analise_poder_estatistico(
    effect_size: float,
    n_amostras: int,
    alpha: float = 0.05
) -> Dict:
    """
    Análise de poder estatístico.
    
    Args:
        effect_size: Tamanho do efeito esperado (Cohen's d)
        n_amostras: Número de amostras por grupo
        alpha: Nível de significância
    
    Returns:
        Dicionário com poder e tamanho de amostra recomendado
    """
    from scipy.stats import nct
    
    # Calcular poder (1 - β)
    # Para teste t independente
    ncp = effect_size * np.sqrt(n_amostras / 2)
    t_crit = stats.t.ppf(1 - alpha/2, 2*n_amostras - 2)
    poder = 1 - nct.cdf(t_crit, 2*n_amostras - 2, ncp)
    
    # Tamanho de amostra recomendado para poder = 0.8
    poder_desejado = 0.8
    n_recomendado = 16  # Aproximação inicial
    
    # Busca binária para encontrar n
    for n in range(5, 1000):
        ncp_test = effect_size * np.sqrt(n / 2)
        t_crit_test = stats.t.ppf(1 - alpha/2, 2*n - 2)
        poder_test = 1 - nct.cdf(t_crit_test, 2*n - 2, ncp_test)
        if poder_test >= poder_desejado:
            n_recomendado = n
            break
    
    return {
        'poder_atual': poder,
        'n_amostras_atual': n_amostras,
        'n_recomendado_80': n_recomendado,
        'suficiente': poder >= 0.8,
        'interpretacao': (
            f"Poder atual: {poder:.2%}. "
            f"{'Suficiente' if poder >= 0.8 else 'Insuficiente'} para detectar efeito."
        )
    }


def comparar_com_baseline(
    baseline: np.ndarray,
    experimental: np.ndarray,
    tipo_teste: str = 'pareado',
    alpha: float = 0.05
) -> Dict:
    """
    Compara configuração experimental com baseline.
    
    Uso típico: Comparar QAOA com ruído vs. sem ruído.
    
    Args:
        baseline: Medições baseline (ex: sem ruído)
        experimental: Medições experimentais (ex: com ruído)
        tipo_teste: 'pareado' ou 'independente'
        alpha: Nível de significância
    
    Returns:
        Dicionário com análise completa
    """
    # Teste estatístico
    if tipo_teste == 'pareado':
        resultado_teste = teste_t_pareado(experimental, baseline, alpha)
    else:
        resultado_teste = teste_t_independente(experimental, baseline, alpha)
    
    # Intervalos de confiança
    media_base, ic_base_inf, ic_base_sup = calcular_intervalo_confianca(baseline)
    media_exp, ic_exp_inf, ic_exp_sup = calcular_intervalo_confianca(experimental)
    
    # Melhoria relativa
    melhoria_absoluta = media_exp - media_base
    melhoria_relativa = (melhoria_absoluta / media_base) * 100
    
    # Análise de poder
    poder = analise_poder_estatistico(
        resultado_teste.effect_size,
        len(baseline),
        alpha
    )
    
    # Interpretação
    if resultado_teste.significativo:
        if melhoria_absoluta > 0:
            conclusao = "✅ RUÍDO BENÉFICO CONFIRMADO estatisticamente"
        else:
            conclusao = "⚠️  RUÍDO PREJUDICIAL confirmado estatisticamente"
    else:
        conclusao = "➖ SEM DIFERENÇA SIGNIFICATIVA entre configurações"
    
    return {
        'teste': resultado_teste,
        'baseline': {
            'media': media_base,
            'ic_95': (ic_base_inf, ic_base_sup),
            'std': np.std(baseline, ddof=1)
        },
        'experimental': {
            'media': media_exp,
            'ic_95': (ic_exp_inf, ic_exp_sup),
            'std': np.std(experimental, ddof=1)
        },
        'melhoria': {
            'absoluta': melhoria_absoluta,
            'relativa_pct': melhoria_relativa
        },
        'poder': poder,
        'conclusao': conclusao
    }


def gerar_relatorio_estatistico(
    df_resultados: pd.DataFrame,
    coluna_metrica: str = 'approx_ratio',
    coluna_grupo: str = 'tipo_ruido',
    baseline: str = 'sem_ruido',
    output_file: Optional[str] = None
) -> str:
    """
    Gera relatório estatístico completo.
    
    Args:
        df_resultados: DataFrame com resultados
        coluna_metrica: Nome da coluna com métrica
        coluna_grupo: Nome da coluna com grupos
        baseline: Nome do grupo baseline
        output_file: Arquivo para salvar relatório
    
    Returns:
        String com relatório formatado
    """
    linhas = []
    linhas.append("="*80)
    linhas.append("RELATÓRIO DE ANÁLISE ESTATÍSTICA - QAOA")
    linhas.append("="*80)
    linhas.append("")
    
    # Estatísticas descritivas
    linhas.append("1. ESTATÍSTICAS DESCRITIVAS")
    linhas.append("-"*80)
    
    grupos = df_resultados.groupby(coluna_grupo)[coluna_metrica]
    for nome, dados in grupos:
        media, ic_inf, ic_sup = calcular_intervalo_confianca(dados.values)
        linhas.append(f"\n{nome}:")
        linhas.append(f"  Média: {media:.4f}")
        linhas.append(f"  IC 95%: [{ic_inf:.4f}, {ic_sup:.4f}]")
        linhas.append(f"  Desvio padrão: {np.std(dados, ddof=1):.4f}")
        linhas.append(f"  N: {len(dados)}")
    
    # ANOVA
    linhas.append("\n\n2. ANOVA (COMPARAÇÃO DE MÚLTIPLOS GRUPOS)")
    linhas.append("-"*80)
    
    grupos_list = [dados.values for _, dados in grupos]
    labels_list = [nome for nome, _ in grupos]
    
    if len(grupos_list) > 2:
        resultado_anova = anova_one_way(grupos_list, labels_list)
        
        linhas.append(f"\nF-statistic: {resultado_anova['anova'].estatistica:.4f}")
        linhas.append(f"P-valor: {resultado_anova['anova'].p_valor:.4e}")
        linhas.append(f"Effect size (η²): {resultado_anova['anova'].effect_size:.4f}")
        linhas.append(f"Significativo: {'SIM' if resultado_anova['anova'].significativo else 'NÃO'}")
        
        if 'posthoc' in resultado_anova:
            linhas.append("\nTestes Post-hoc (comparações pareadas):")
            for ph in resultado_anova['posthoc']:
                linhas.append(f"  {ph['comparacao']}:")
                linhas.append(f"    p-valor: {ph['p_valor']:.4e}")
                linhas.append(f"    Cohen's d: {ph['effect_size']:.4f}")
                linhas.append(f"    Significativo: {'SIM' if ph['significativo'] else 'NÃO'}")
    
    # Comparações com baseline
    linhas.append("\n\n3. COMPARAÇÕES COM BASELINE")
    linhas.append("-"*80)
    
    if baseline in labels_list:
        baseline_dados = df_resultados[
            df_resultados[coluna_grupo] == baseline
        ][coluna_metrica].values
        
        for nome, dados in grupos:
            if nome != baseline:
                comp = comparar_com_baseline(
                    baseline_dados,
                    dados.values,
                    tipo_teste='independente'
                )
                
                linhas.append(f"\n{baseline} vs {nome}:")
                linhas.append(f"  {comp['conclusao']}")
                linhas.append(f"  Melhoria: {comp['melhoria']['relativa_pct']:+.2f}%")
                linhas.append(f"  P-valor: {comp['teste'].p_valor:.4e}")
                linhas.append(f"  Effect size: {comp['teste'].effect_size:.4f} ({comp['teste'].interpretacao})")
                linhas.append(f"  Poder: {comp['poder']['poder_atual']:.2%}")
    
    linhas.append("\n" + "="*80)
    
    relatorio = "\n".join(linhas)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(relatorio)
        print(f"✓ Relatório salvo em: {output_file}")
    
    return relatorio


if __name__ == "__main__":
    print("Testando módulo de análise estatística...")
    
    # Teste 1: Teste t pareado
    print("\n1. Teste t pareado:")
    sem_ruido = np.random.normal(0.85, 0.05, 20)
    com_ruido = sem_ruido + np.random.normal(0.03, 0.02, 20)
    
    resultado = teste_t_pareado(com_ruido, sem_ruido)
    print(f"   Estatística t: {resultado.estatistica:.4f}")
    print(f"   P-valor: {resultado.p_valor:.4e}")
    print(f"   Cohen's d: {resultado.effect_size:.4f}")
    print(f"   Significativo: {resultado.significativo}")
    print(f"   Interpretação: {resultado.interpretacao}")
    
    # Teste 2: ANOVA
    print("\n2. ANOVA one-way:")
    grupo1 = np.random.normal(0.80, 0.05, 15)
    grupo2 = np.random.normal(0.85, 0.05, 15)
    grupo3 = np.random.normal(0.88, 0.05, 15)
    
    resultado_anova = anova_one_way(
        [grupo1, grupo2, grupo3],
        ['sem_ruido', 'ruido_baixo', 'ruido_medio']
    )
    
    print(f"   F-statistic: {resultado_anova['anova'].estatistica:.4f}")
    print(f"   P-valor: {resultado_anova['anova'].p_valor:.4e}")
    print(f"   Significativo: {resultado_anova['anova'].significativo}")
    
    # Teste 3: Análise de poder
    print("\n3. Análise de poder estatístico:")
    poder = analise_poder_estatistico(effect_size=0.5, n_amostras=20)
    print(f"   Poder atual: {poder['poder_atual']:.2%}")
    print(f"   N recomendado: {poder['n_recomendado_80']}")
    print(f"   {poder['interpretacao']}")
    
    print("\n✓ Todos os testes concluídos")
