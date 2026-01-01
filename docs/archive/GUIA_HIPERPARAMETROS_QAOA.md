# Guia de Otimização de Hiperparâmetros QAOA

## 🎯 Objetivo

Este guia mostra como encontrar os hiperparâmetros ótimos (níveis de ruído, tipos de ruído, profundidade p, etc.) para QAOA usando o framework desenvolvido.

---


## 📊 Hiperparâmetros do QAOA

### Hiperparâmetros Principais

| Parâmetro | Descrição | Faixa Típica | Impacto |
|-----------|-----------|--------------|---------|
| **n_qubits** | Número de qubits | 10-100 | Tamanho do problema |
| **p_layers** | Profundidade QAOA | 1-10 | Qualidade da solução vs. profundidade |
| **tipo_ruido** | Tipo de canal de ruído | 4 tipos | Mecanismo de interação |
| **nivel_ruido** | Taxa de erro | 0.0-0.05 | Intensidade do efeito |
| **shots** | Medições por execução | 512-8192 | Precisão estatística |
| **max_iter** | Iterações máximas | 50-200 | Convergência |
| **otimizador** | Otimizador clássico | COBYLA, SLSQP, etc. | Estratégia de busca |

### Hiperparâmetros de Ruído

Tipos disponíveis:

- `'sem_ruido'`: Baseline (ideal)
- `'depolarizing'`: Erro genérico (todas direções)
- `'amplitude_damping'`: Relaxação de energia (T1)
- `'phase_damping'`: Perda de coerência (T2)
- `'thermal'`: Modelo realista (T1 + T2)


---


## 🔍 Método 1: Grid Search

### Uso Básico

```python
from framework_qaoa_100qubits import (
    ConstrutorCircuitoQAOA,
    AnalisadorHiperparametrosQAOA
)

# 1. Criar problema
construtor = ConstrutorCircuitoQAOA(n_qubits=30, p_layers=3)
grafo = construtor.criar_grafo_aleatorio(densidade=0.2)

# 2. Configurar analisador
analisador = AnalisadorHiperparametrosQAOA(
    pasta_resultados='resultados_grid_search'
)

# 3. Definir espaço de busca
niveis_ruido = [0.0, 0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01]
tipos_ruido = ['sem_ruido', 'depolarizing', 'amplitude_damping', 'phase_damping']

# 4. Executar grid search
df_resultados = analisador.grid_search_ruido(
    grafo=grafo,
    niveis_ruido=niveis_ruido,
    tipos_ruido=tipos_ruido,
    p_layers=3,
    n_repeticoes=5  # Repetições para análise estatística
)

# 5. Analisar
print(df_resultados.groupby(['tipo_ruido', 'nivel_ruido']).agg({
    'energia_final': ['mean', 'std', 'min']
}))

```text

### Encontrar Melhor Configuração

```python

# Melhor configuração global
melhor_idx = df_resultados['energia_final'].idxmin()
melhor = df_resultados.loc[melhor_idx]

print(f"Melhor configuração:")
print(f"  Tipo ruído:  {melhor['tipo_ruido']}")
print(f"  Nível ruído: {melhor['nivel_ruido']:.4f}")
print(f"  Energia:     {melhor['energia_final']:.4f}")

# Melhor para cada tipo de ruído
for tipo in tipos_ruido:
    df_tipo = df_resultados[df_resultados['tipo_ruido'] == tipo]
    melhor_tipo = df_tipo.loc[df_tipo['energia_final'].idxmin()]
    print(f"\n{tipo}:")
    print(f"  Melhor nível: {melhor_tipo['nivel_ruido']:.4f}")
    print(f"  Energia:      {melhor_tipo['energia_final']:.4f}")

```text

### Visualizar Resultados

```python
from framework_qaoa_100qubits import VisualizadorQAOA

visualizador = VisualizadorQAOA()
visualizador.plotar_comparacao_ruido(
    df_resultados,
    salvar='resultados_grid_search/comparacao.html'
)

```text

---


## 🎲 Método 2: Otimização Bayesiana (Optuna)

### Uso Básico

```python
from framework_qaoa_100qubits import AnalisadorHiperparametrosQAOA

analisador = AnalisadorHiperparametrosQAOA(
    pasta_resultados='resultados_bayesiana'
)

# Otimização Bayesiana
resultado_bayes = analisador.otimizacao_bayesiana(
    grafo=grafo,
    n_trials=50  # Número de trials
)

# Melhores parâmetros
print("Melhores hiperparâmetros encontrados:")
print(f"  Tipo ruído:  {resultado_bayes['best_params']['tipo_ruido']}")
print(f"  Nível ruído: {resultado_bayes['best_params']['nivel_ruido']:.4f}")
print(f"  P-layers:    {resultado_bayes['best_params']['p_layers']}")
print(f"\nMelhor energia: {resultado_bayes['best_value']:.4f}")

```text

### Customizar Espaço de Busca

Para controle mais fino, edite a função `otimizacao_bayesiana` em `framework_qaoa_100qubits.py`:

```python
def objetivo_optuna(trial):

    # Sugerir hiperparâmetros
    tipo_ruido = trial.suggest_categorical(
        'tipo_ruido',
        ['depolarizing', 'amplitude_damping', 'phase_damping', 'sem_ruido']
    )
    
    if tipo_ruido != 'sem_ruido':

        # Log-scale para nível de ruído
        nivel_ruido = trial.suggest_float('nivel_ruido', 0.0001, 0.01, log=True)
    else:
        nivel_ruido = 0.0
    
    # Profundidade QAOA
    p_layers = trial.suggest_int('p_layers', 1, 5)
    
    # Shots (opcional)
    shots = trial.suggest_categorical('shots', [512, 1024, 2048])
    
    # ... resto do código

```text

### Analisar Trials

```python
study = resultado_bayes['study']

# Importância dos parâmetros
import optuna.visualization as vis

fig = vis.plot_param_importances(study)
fig.show()

# Histórico de otimização
fig = vis.plot_optimization_history(study)
fig.show()

# Slice plot
fig = vis.plot_slice(study)
fig.show()

```text

---


## 📐 Estratégias de Busca

### 1. Busca Grosseira → Fina (Coarse-to-Fine)

```python

# Fase 1: Grid grosseiro
niveis_fase1 = [0.0, 0.001, 0.01]
df_fase1 = analisador.grid_search_ruido(
    grafo=grafo,
    niveis_ruido=niveis_fase1,
    tipos_ruido=['sem_ruido', 'depolarizing'],
    p_layers=2,
    n_repeticoes=3
)

# Identificar região promissora
melhor_fase1 = df_fase1.loc[df_fase1['energia_final'].idxmin()]
nivel_central = melhor_fase1['nivel_ruido']

# Fase 2: Grid fino ao redor do ótimo
niveis_fase2 = [
    nivel_central * 0.5,
    nivel_central * 0.75,
    nivel_central,
    nivel_central * 1.25,
    nivel_central * 1.5
]
df_fase2 = analisador.grid_search_ruido(
    grafo=grafo,
    niveis_ruido=niveis_fase2,
    tipos_ruido=[melhor_fase1['tipo_ruido']],
    p_layers=3,
    n_repeticoes=5
)

```text

### 2. Busca Multi-Objetivo

Otimizar múltiplos objetivos simultaneamente:

```python
def objetivo_multi(params, grafo):
    resultado = otimizador.otimizar(grafo)
    
    # Objetivo 1: Energia (minimizar)
    energia = resultado.energia_final
    
    # Objetivo 2: Tempo (minimizar)
    tempo = resultado.tempo_execucao
    
    # Objetivo 3: Convergência (maximizar ~ minimizar iterações)
    convergencia = resultado.iteracoes
    
    # Combinar (ponderado)
    score = energia + 0.1 * tempo + 0.01 * convergencia
    return score

```text

### 3. Busca Sequencial por Tipo

```python
resultados_por_tipo = {}

for tipo in ['depolarizing', 'amplitude_damping', 'phase_damping']:
    print(f"\nOtimizando {tipo}...")
    
    # Grid específico para este tipo
    df = analisador.grid_search_ruido(
        grafo=grafo,
        niveis_ruido=np.logspace(-4, -2, 10),  # 10 níveis em log-scale
        tipos_ruido=[tipo],
        p_layers=3,
        n_repeticoes=5
    )
    
    resultados_por_tipo[tipo] = df

# Comparar melhores de cada tipo
for tipo, df in resultados_por_tipo.items():
    melhor = df.loc[df['energia_final'].idxmin()]
    print(f"{tipo}: {melhor['energia_final']:.4f} @ {melhor['nivel_ruido']:.4f}")

```text

---


## 📊 Análise Estatística

### Teste de Significância

```python
from scipy.stats import ttest_ind

# Comparar com vs. sem ruído
energias_sem = df_resultados[
    df_resultados['tipo_ruido'] == 'sem_ruido'
]['energia_final'].values

energias_com = df_resultados[
    df_resultados['tipo_ruido'] == 'depolarizing'
]['energia_final'].values

stat, p_value = ttest_ind(energias_sem, energias_com)

print(f"t-statistic: {stat:.4f}")
print(f"p-value: {p_value:.4f}")

if p_value < 0.05:
    print("✅ Diferença estatisticamente significativa!")
else:
    print("⚠️  Diferença não significativa")

```text

### Effect Size (Cohen's d)

```python
import numpy as np

def cohens_d(group1, group2):
    mean1, mean2 = np.mean(group1), np.mean(group2)
    std1, std2 = np.std(group1, ddof=1), np.std(group2, ddof=1)
    n1, n2 = len(group1), len(group2)
    
    pooled_std = np.sqrt(((n1-1)*std1**2 + (n2-1)*std2**2) / (n1+n2-2))
    d = (mean1 - mean2) / pooled_std
    return d

d = cohens_d(energias_sem, energias_com)
print(f"Cohen's d: {d:.4f}")

if abs(d) < 0.2:
    print("Efeito pequeno")
elif abs(d) < 0.5:
    print("Efeito médio")
else:
    print("Efeito grande")

```text

---


## 🎯 Receitas Práticas

### Receita 1: Encontrar Ruído Benéfico Rapidamente

```python

# 1. Teste rápido com poucos qubits
construtor = ConstrutorCircuitoQAOA(n_qubits=20, p_layers=2)
grafo = construtor.criar_grafo_aleatorio(densidade=0.3)

# 2. Grid pequeno
niveis = [0.0, 0.0005, 0.001, 0.002]
tipos = ['sem_ruido', 'depolarizing', 'phase_damping']

analisador = AnalisadorHiperparametrosQAOA()
df = analisador.grid_search_ruido(
    grafo=grafo,
    niveis_ruido=niveis,
    tipos_ruido=tipos,
    p_layers=2,
    n_repeticoes=3
)

# 3. Verificar se há benefício
baseline = df[df['tipo_ruido'] == 'sem_ruido']['energia_final'].mean()
melhor = df[df['tipo_ruido'] != 'sem_ruido']['energia_final'].min()

if melhor < baseline:
    print(f"✅ Ruído benéfico encontrado!")
    print(f"   Baseline: {baseline:.4f}")
    print(f"   Melhor:   {melhor:.4f}")
else:
    print("⚠️  Ruído não foi benéfico neste teste")

```text

### Receita 2: Otimização Completa para Publicação

```python

# 1. Problema real (mais qubits)
construtor = ConstrutorCircuitoQAOA(n_qubits=50, p_layers=4)
grafo = construtor.criar_grafo_aleatorio(densidade=0.15)

# 2. Grid search detalhado
niveis = np.logspace(-4, -2, 15)  # 15 pontos em log-scale
tipos = ['sem_ruido', 'depolarizing', 'amplitude_damping', 'phase_damping']

analisador = AnalisadorHiperparametrosQAOA(
    pasta_resultados='resultados_publicacao'
)

df = analisador.grid_search_ruido(
    grafo=grafo,
    niveis_ruido=niveis,
    tipos_ruido=tipos,
    p_layers=4,
    n_repeticoes=10  # Muitas repetições para estatística
)

# 3. Análise estatística completa
# ... ANOVA, post-hoc, effect sizes ...

# 4. Visualizações de alta qualidade
visualizador = VisualizadorQAOA()
visualizador.plotar_comparacao_ruido(
    df,
    salvar='resultados_publicacao/figura_principal.html'
)

# 5. Salvar tabela para paper
resumo = df.groupby(['tipo_ruido', 'nivel_ruido']).agg({
    'energia_final': ['mean', 'std', 'min', 'max']
})
resumo.to_csv('resultados_publicacao/tabela_resultados.csv')

```text

### Receita 3: Validação Cruzada

```python

# Testar em múltiplos grafos para generalização
n_grafos = 5
resultados_todos = []

for i in range(n_grafos):
    print(f"\nGrafo {i+1}/{n_grafos}")
    
    # Grafo diferente
    construtor = ConstrutorCircuitoQAOA(n_qubits=30, p_layers=3, seed=42+i)
    grafo = construtor.criar_grafo_aleatorio(densidade=0.2)
    
    # Grid search
    df = analisador.grid_search_ruido(
        grafo=grafo,
        niveis_ruido=[0.0, 0.0005, 0.001, 0.002],
        tipos_ruido=['sem_ruido', 'depolarizing'],
        p_layers=3,
        n_repeticoes=3
    )
    
    df['grafo_id'] = i
    resultados_todos.append(df)

# Concatenar
df_total = pd.concat(resultados_todos, ignore_index=True)

# Analisar consistência
print("\nConsistência entre grafos:")
print(df_total.groupby(['tipo_ruido', 'nivel_ruido', 'grafo_id'])['energia_final'].mean())

```

---


## 💡 Dicas e Boas Práticas

### Performance

1. **Comece pequeno**: Teste com 10-20 qubits primeiro
2. **Use p-layers = 2-3**: Profundidade maior nem sempre é melhor
3. **Log-scale para ruído**: `np.logspace(-4, -2, 10)` ao invés de `np.linspace`
4. **Shots moderados**: 512-1024 é suficiente para exploração


### Estatística

1. **Sempre use n_repeticoes ≥ 3**: Mínimo para estatística básica
2. **n_repeticoes ≥ 10 para publicação**: Análise robusta
3. **Calcule intervalos de confiança**: Use `df.groupby(...).agg(['mean', 'std'])`
4. **Teste significância**: t-test, ANOVA conforme apropriado


### Reprodutibilidade

1. **Fixe seeds**: `seed=42` em todos os lugares
2. **Documente configuração**: Salve ConfigQAOA como JSON
3. **Versione resultados**: Use timestamps nos nomes de arquivos
4. **Salve grafos**: `np.save('grafo.npy', grafo)` para reuso


---


## 📚 Referências

1. **Optuna**: <https://optuna.org/>
2. **Scipy Optimize**: <https://docs.scipy.org/doc/scipy/reference/optimize.html>
3. **QAOA Hyperparameter Optimization**: Zhou et al. (2020), PRX Quantum


---


**Última atualização**: 2025-12-26

