"""
Tarefa 8: Visualizações QAOA
Gera visualizações específicas para análise de resultados QAOA.

Referências:
- Plotly Documentation: https://plotly.com/python/
- Matplotlib Documentation: https://matplotlib.org/
"""

from typing import List, Dict, Optional, Any
import numpy as np
import pandas as pd
from pathlib import Path

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("⚠️  Plotly não disponível")

try:
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    print("⚠️  Matplotlib não disponível")


def visualizar_convergencia_qaoa(
    historico_energia: List[float],
    titulo: str = "Convergência QAOA",
    salvar: Optional[str] = None
) -> Optional[go.Figure]:
    """
    Visualiza convergência da energia durante otimização QAOA.
    
    Args:
        historico_energia: Lista de energias por iteração
        titulo: Título do gráfico
        salvar: Caminho para salvar figura (HTML ou PNG)
    
    Returns:
        Figura Plotly (se disponível)
    """
    if not PLOTLY_AVAILABLE:
        print("Plotly não disponível. Não é possível gerar visualização.")
        return None
    
    fig = go.Figure()
    
    # Linha de convergência
    fig.add_trace(go.Scatter(
        y=historico_energia,
        mode='lines+markers',
        name='Energia',
        line=dict(color='blue', width=2),
        marker=dict(size=6, symbol='circle')
    ))
    
    # Melhor energia (linha horizontal)
    melhor_energia = min(historico_energia)
    fig.add_hline(
        y=melhor_energia,
        line_dash="dash",
        line_color="red",
        annotation_text=f"Mínimo: {melhor_energia:.4f}",
        annotation_position="right"
    )
    
    fig.update_layout(
        title=titulo,
        xaxis_title="Iteração",
        yaxis_title="Energia",
        template='plotly_white',
        width=900,
        height=500,
        font=dict(size=12)
    )
    
    if salvar:
        if salvar.endswith('.html'):
            fig.write_html(salvar)
        elif salvar.endswith('.png'):
            fig.write_image(salvar, width=900, height=500, scale=2)
        print(f"✓ Gráfico salvo: {salvar}")
    
    return fig


def visualizar_histograma_solucoes(
    contagens: Dict[str, int],
    top_n: int = 20,
    titulo: str = "Distribuição de Soluções QAOA",
    salvar: Optional[str] = None
) -> Optional[go.Figure]:
    """
    Visualiza histograma das soluções mais frequentes.
    
    Args:
        contagens: Dicionário {bitstring: count}
        top_n: Número de soluções a mostrar
        titulo: Título do gráfico
        salvar: Caminho para salvar
    
    Returns:
        Figura Plotly
    """
    if not PLOTLY_AVAILABLE:
        return None
    
    # Ordenar e pegar top N
    contagens_ordenadas = sorted(contagens.items(), key=lambda x: x[1], reverse=True)[:top_n]
    bitstrings, counts = zip(*contagens_ordenadas) if contagens_ordenadas else ([], [])
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=list(bitstrings),
        y=list(counts),
        marker=dict(
            color=list(counts),
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title="Contagem")
        )
    ))
    
    fig.update_layout(
        title=titulo,
        xaxis_title="Bitstring (solução)",
        yaxis_title="Frequência",
        template='plotly_white',
        width=1000,
        height=500,
        xaxis=dict(tickangle=-45)
    )
    
    if salvar:
        fig.write_html(salvar) if salvar.endswith('.html') else fig.write_image(salvar)
        print(f"✓ Histograma salvo: {salvar}")
    
    return fig


def visualizar_comparacao_ruido(
    df_resultados: pd.DataFrame,
    metrica: str = 'approx_ratio',
    titulo: str = "Impacto do Ruído no QAOA",
    salvar: Optional[str] = None
) -> Optional[go.Figure]:
    """
    Visualiza comparação de diferentes níveis e tipos de ruído.
    
    Args:
        df_resultados: DataFrame com colunas [tipo_ruido, nivel_ruido, metrica]
        metrica: Nome da coluna com métrica (approx_ratio, energia_final)
        titulo: Título do gráfico
        salvar: Caminho para salvar
    
    Returns:
        Figura Plotly
    """
    if not PLOTLY_AVAILABLE:
        return None
    
    fig = go.Figure()
    
    # Para cada tipo de ruído, criar uma linha
    for tipo in df_resultados['tipo_ruido'].unique():
        df_tipo = df_resultados[df_resultados['tipo_ruido'] == tipo]
        
        # Agrupar por nível de ruído
        grouped = df_tipo.groupby('nivel_ruido').agg({
            metrica: ['mean', 'std']
        }).reset_index()
        
        fig.add_trace(go.Scatter(
            x=grouped['nivel_ruido'],
            y=grouped[metrica]['mean'],
            error_y=dict(
                type='data',
                array=grouped[metrica]['std'],
                visible=True
            ),
            mode='lines+markers',
            name=tipo,
            line=dict(width=2),
            marker=dict(size=8)
        ))
    
    fig.update_layout(
        title=titulo,
        xaxis_title="Nível de Ruído",
        yaxis_title=metrica.replace('_', ' ').title(),
        template='plotly_white',
        width=1000,
        height=600,
        legend=dict(x=0.02, y=0.98),
        xaxis_type='log' if df_resultados['nivel_ruido'].min() > 0 else 'linear'
    )
    
    if salvar:
        fig.write_html(salvar) if salvar.endswith('.html') else fig.write_image(salvar)
        print(f"✓ Comparação salva: {salvar}")
    
    return fig


def visualizar_heatmap_hiperparametros(
    df_optuna: pd.DataFrame,
    param_x: str,
    param_y: str,
    metrica: str = 'value',
    titulo: str = "Mapa de Hiperparâmetros",
    salvar: Optional[str] = None
) -> Optional[go.Figure]:
    """
    Cria heatmap de 2 hiperparâmetros vs métrica.
    
    Args:
        df_optuna: DataFrame de trials do Optuna
        param_x: Nome do parâmetro no eixo X
        param_y: Nome do parâmetro no eixo Y
        metrica: Coluna com métrica
        titulo: Título do gráfico
        salvar: Caminho para salvar
    
    Returns:
        Figura Plotly
    """
    if not PLOTLY_AVAILABLE:
        return None
    
    # Criar pivot table
    pivot = df_optuna.pivot_table(
        values=metrica,
        index=param_y,
        columns=param_x,
        aggfunc='mean'
    )
    
    fig = go.Figure(data=go.Heatmap(
        z=pivot.values,
        x=pivot.columns,
        y=pivot.index,
        colorscale='Viridis',
        colorbar=dict(title=metrica.replace('_', ' ').title())
    ))
    
    fig.update_layout(
        title=titulo,
        xaxis_title=param_x.replace('_', ' ').title(),
        yaxis_title=param_y.replace('_', ' ').title(),
        template='plotly_white',
        width=800,
        height=600
    )
    
    if salvar:
        fig.write_html(salvar) if salvar.endswith('.html') else fig.write_image(salvar)
        print(f"✓ Heatmap salvo: {salvar}")
    
    return fig


def visualizar_painel_qaoa(
    resultados: Dict[str, Any],
    pasta_saida: Optional[str] = None
) -> Dict[str, go.Figure]:
    """
    Cria painel completo de visualizações QAOA.
    
    Args:
        resultados: Dicionário com todos os resultados do experimento
        pasta_saida: Pasta para salvar visualizações
    
    Returns:
        Dicionário com todas as figuras geradas
    """
    if pasta_saida:
        Path(pasta_saida).mkdir(parents=True, exist_ok=True)
    
    figuras = {}
    
    # 1. Convergência
    if 'historico_energia' in resultados:
        fig = visualizar_convergencia_qaoa(
            resultados['historico_energia'],
            salvar=f"{pasta_saida}/convergencia.html" if pasta_saida else None
        )
        figuras['convergencia'] = fig
    
    # 2. Histograma de soluções
    if 'contagens' in resultados:
        fig = visualizar_histograma_solucoes(
            resultados['contagens'],
            salvar=f"{pasta_saida}/histograma_solucoes.html" if pasta_saida else None
        )
        figuras['histograma'] = fig
    
    # 3. Comparação de ruído
    if 'df_resultados' in resultados:
        fig = visualizar_comparacao_ruido(
            resultados['df_resultados'],
            salvar=f"{pasta_saida}/comparacao_ruido.html" if pasta_saida else None
        )
        figuras['comparacao_ruido'] = fig
    
    print(f"\n✓ Painel completo gerado: {len(figuras)} visualizações")
    if pasta_saida:
        print(f"  Salvo em: {pasta_saida}")
    
    return figuras


def gerar_relatorio_visual_html(
    resultados: Dict[str, Any],
    arquivo_saida: str
):
    """
    Gera relatório HTML completo com todas as visualizações.
    
    Args:
        resultados: Dicionário com resultados
        arquivo_saida: Caminho do arquivo HTML de saída
    """
    html_parts = []
    
    # Cabeçalho
    html_parts.append("""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Relatório QAOA - Análise de Ruído Benéfico</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            h1 { color: #2c3e50; }
            h2 { color: #3498db; margin-top: 40px; }
            .summary { background: #ecf0f1; padding: 20px; border-radius: 5px; }
            .metric { font-size: 24px; font-weight: bold; color: #27ae60; }
        </style>
    </head>
    <body>
        <h1>📊 Relatório QAOA - Análise de Ruído Benéfico</h1>
    """)
    
    # Resumo executivo
    html_parts.append("""
        <div class="summary">
            <h2>Resumo Executivo</h2>
    """)
    
    if 'best_approx_ratio' in resultados:
        html_parts.append(f"""
            <p>Melhor Approximation Ratio: <span class="metric">{resultados['best_approx_ratio']:.4f}</span></p>
        """)
    
    if 'best_config' in resultados:
        html_parts.append("<p><strong>Melhor configuração:</strong></p><ul>")
        for k, v in resultados['best_config'].items():
            html_parts.append(f"<li>{k}: {v}</li>")
        html_parts.append("</ul>")
    
    html_parts.append("</div>")
    
    # Visualizações
    figuras = visualizar_painel_qaoa(resultados)
    
    for nome, fig in figuras.items():
        if fig:
            html_parts.append(f"<h2>{nome.replace('_', ' ').title()}</h2>")
            html_parts.append(fig.to_html(full_html=False, include_plotlyjs='cdn'))
    
    # Rodapé
    html_parts.append("""
        <hr>
        <p><em>Gerado automaticamente pelo Framework QAOA</em></p>
    </body>
    </html>
    """)
    
    # Salvar
    with open(arquivo_saida, 'w', encoding='utf-8') as f:
        f.write('\n'.join(html_parts))
    
    print(f"✓ Relatório HTML gerado: {arquivo_saida}")


if __name__ == "__main__":
    print("Testando módulo de visualizações QAOA...")
    
    if PLOTLY_AVAILABLE:
        print("✓ Plotly disponível")
        
        # Teste 1: Convergência
        print("\n1. Testando visualização de convergência...")
        historico = [-10.0, -12.5, -13.8, -14.2, -14.5, -14.6, -14.6]
        fig = visualizar_convergencia_qaoa(historico)
        if fig:
            print("   ✓ Gráfico de convergência criado")
        
        # Teste 2: Histograma
        print("\n2. Testando histograma de soluções...")
        contagens = {
            '0101': 150,
            '1010': 145,
            '0011': 90,
            '1100': 85,
            '0000': 30
        }
        fig = visualizar_histograma_solucoes(contagens)
        if fig:
            print("   ✓ Histograma criado")
        
        # Teste 3: Comparação de ruído
        print("\n3. Testando comparação de ruído...")
        df_test = pd.DataFrame({
            'tipo_ruido': ['sem_ruido']*3 + ['depolarizing']*3,
            'nivel_ruido': [0.0, 0.0, 0.0, 0.001, 0.005, 0.01],
            'approx_ratio': [0.85, 0.87, 0.86, 0.88, 0.90, 0.87]
        })
        fig = visualizar_comparacao_ruido(df_test)
        if fig:
            print("   ✓ Gráfico de comparação criado")
        
        print("\n✓ Todos os testes de visualização concluídos")
    
    else:
        print("❌ Plotly não disponível para testes")
