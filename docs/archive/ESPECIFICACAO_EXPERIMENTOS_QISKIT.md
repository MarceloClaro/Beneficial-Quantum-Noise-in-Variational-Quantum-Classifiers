# Framework Qiskit - Especificação Completa de Experimentos

## 📋 Parâmetros Experimentais

Este documento descreve a configuração **exata** dos experimentos do framework Qiskit, replicando o framework PennyLane.

### Grid Search Completo

O framework Qiskit suporta o mesmo grid search extensivo do PennyLane:

**Total de combinações possíveis**: 5 × 9 × 4 × 6 × 9 × 5 = **48,600 experimentos**


## 📊 Datasets (5)

| Dataset | Amostras | Features | Classes | Descrição |
|---------|----------|----------|---------|-----------|
| **moons** | 200 | 2 | 2 | Duas luas entrelaçadas (sintético) |
| **circles** | 200 | 2 | 2 | Círculos concêntricos (sintético) |
| **iris** | 150 | 2 | 2 | Iris dataset (binário: setosa vs resto) |
| **breast_cancer** | 569 | 4 | 2 | Wisconsin Breast Cancer (PCA) |
| **wine** | 178 | 4 | 2 | Wine dataset (binário: classe 0 vs resto) |

## 🏗️ Arquiteturas de Circuitos (9)

| ID | Nome | Parâmetros | Entanglement | Uso |
|----|------|------------|--------------|-----|
| 1 | **basico** | 3n×l | Cadeia linear | Baseline |
| 2 | **strongly_entangling** | 3n×l | All-to-all | Máximo entrelaçamento |
| 3 | **hardware_efficient** | n×l | Linear | Otimizado para hardware |
| 4 | **alternating_layers** | 2n×l | Alternado | Camadas pares/ímpares |
| 5 | **brickwork** | 2n×l | Padrão brickwork | Médio entrelaçamento |
| 6 | **random_entangling** | 2n×l | Aleatório | Exploração |
| 7 | **tree** | 2n×l | Árvore binária | Hierárquico |
| 8 | **star_entanglement** | 2n×l | Estrela (hub central) | Hub-and-spoke |
| 9 | **qaoa** | 2n×l | QAOA-like | Mixing + Problem |

**Onde**: n = n_qubits, l = n_camadas


## 🎲 Estratégias de Inicialização (4)

| Estratégia | Descrição | Base |
|------------|-----------|------|
| **matematico** | Constantes matemáticas | π, e, φ, √2, ln2, γ |
| **quantico** | Constantes físicas | ℏ, α (fine-structure), R∞ |
| **aleatorio** | Uniforme [-π, π] | Baseline |
| **fibonacci_spiral** | Espiral de Fibonacci | Razão áurea |

## 🔬 Tipos de Ruído Quântico (6)

| Tipo | Código | Descrição | Modelo Físico |
|------|--------|-----------|---------------|
| **Sem ruído** | `sem_ruido` | Simulação ideal | Baseline |
| **Depolarizante** | `depolarizante` | Ruído isotrópico | Erro simétrico X, Y, Z |
| **Amplitude Damping** | `amplitude_damping` | Perda de energia | Relaxação T1 (|1⟩→|0⟩) |
| **Phase Damping** | `phase_damping` | Perda de coerência | Decoerência T2 |
| **Crosstalk** | `crosstalk` | Interferência entre qubits | Erro aumentado em 2-qubit gates |
| **Correlacionado** | `correlacionado` | Erro correlacionado | Combinação dep+phase correlacionados |

### Implementação dos Modelos de Ruído

```python

# Depolarizante: X, Y, Z com mesma probabilidade
error_1q = depolarizing_error(γ, 1)
error_2q = depolarizing_error(2γ, 2)

# Amplitude Damping: |1⟩ → |0⟩
error = amplitude_damping_error(γ)

# Phase Damping: Perda de fase
error = phase_damping_error(γ)

# Crosstalk: Maior erro em 2-qubit
error_1q = depolarizing_error(0.3γ, 1)
error_2q = depolarizing_error(1.5γ, 2)

# Correlacionado: Múltiplos erros combinados
error_dep = depolarizing_error(0.6γ, 1)
error_phase = phase_damping_error(0.4γ)

```text

## 📉 Níveis de Ruído (9)

| Nível | γ | Regime | Aplicação |
|-------|---|--------|-----------|
| 1 | **0.0000** | Sem ruído | Baseline ideal |
| 2 | **0.0025** | Muito baixo | NISQ alta qualidade |
| 3 | **0.0050** | Baixo | **Regime benéfico típico** |
| 4 | **0.0075** | Médio-baixo | NISQ média qualidade |
| 5 | **0.0100** | Médio | Limite benéfico |
| 6 | **0.0125** | Médio-alto | Transição prejudicial |
| 7 | **0.0150** | Alto | Prejudicial |
| 8 | **0.0175** | Muito alto | Fortemente prejudicial |
| 9 | **0.0200** | Crítico | Perda significativa |

**Região de interesse**: γ ∈ [0.001, 0.01] (ruído potencialmente benéfico)


## 🎯 Seeds Aleatórias (5)

| Seed | Uso |
|------|-----|
| **42** | Principal (padrão) |
| **43** | Replicação 1 |
| **44** | Replicação 2 |
| **45** | Replicação 3 |
| **46** | Replicação 4 |

**Propósito**: Validação estatística e cálculo de intervalos de confiança


## ⚙️ Hiperparâmetros Fixos

| Parâmetro | Valor | Justificativa |
|-----------|-------|---------------|
| **n_qubits** | 4 | Balanceamento simulação/expressividade |
| **n_camadas** | 2 | Profundidade moderada |
| **n_epocas** | 10 | Rápido (PennyLane: 20) |
| **batch_size** | 32 | Mini-batch padrão |
| **taxa_aprendizado** | 0.01 | Learning rate conservador |
| **shots** | 1024 | Precisão estatística |
| **otimizador** | Adam | Adaptativo |

## 🚀 Modos de Execução

### Modo Rápido (Demonstração)

```bash
python executar_grid_search_qiskit.py --modo rapido

```text

- 1 dataset (moons)
- 2 arquiteturas (basico, strongly_entangling)
- 2 init (aleatorio, matematico)
- 2 ruídos (sem_ruido, depolarizante)
- 3 níveis (0.0, 0.005, 0.01)
- 2 seeds (42, 43)


**Total**: ~24 experimentos (~1 hora)


### Modo Médio (Validação)

```bash
python executar_grid_search_qiskit.py --modo medio

```text

- 3 datasets (moons, circles, iris)
- 5 arquiteturas (primeiras 5)
- 4 init (todas)
- 4 ruídos (primeiros 4)
- 5 níveis (pares)
- 3 seeds (42, 43, 44)


**Total**: ~1,800 experimentos (~6-8 horas)


### Modo Completo (Produção)

```bash
python executar_grid_search_qiskit.py --modo completo

```text

- 5 datasets (todos)
- 9 arquiteturas (todas)
- 4 init (todas)
- 6 ruídos (todos)
- 9 níveis (todos)
- 5 seeds (todas)


**Total**: ~48,600 experimentos (~5-7 dias)


## 📊 Formato de Saída

### CSV de Resultados

Cada linha contém:

```

dataset, arquitetura, estrategia_init, tipo_ruido, nivel_ruido, seed,
n_qubits, n_camadas, n_epocas, acuracia_treino, acuracia_teste,
tempo_treino, framework, shots

```text

### Exemplo

```csv
dataset,arquitetura,estrategia_init,tipo_ruido,nivel_ruido,seed,acuracia_teste
moons,strongly_entangling,matematico,phase_damping,0.005,42,0.6333
moons,strongly_entangling,matematico,phase_damping,0.005,43,0.6500
moons,strongly_entangling,matematico,phase_damping,0.005,44,0.6167

```text

## 🔍 Análises Estatísticas

Com 5 seeds por configuração, é possível calcular:

1. **Média e desvio padrão** de acurácia
2. **Intervalos de confiança 95%** (t-Student)
3. **Testes de hipótese** (ANOVA, t-test)
4. **Effect sizes** (Cohen's d, Hedges' g)


## 🎯 Comparação com PennyLane

| Aspecto | PennyLane | Qiskit |
|---------|-----------|--------|
| **Datasets** | 5 ✅ | 5 ✅ |
| **Arquiteturas** | 9 ✅ | 9 ✅ |
| **Init** | 4 ✅ | 4 ✅ |
| **Ruído** | 6 ✅ | 6 ✅ |
| **Níveis** | 9 ✅ | 9 ✅ |
| **Seeds** | 5 ✅ | 5 ✅ |
| **Hiperparâmetros** | Idênticos ✅ | Idênticos ✅ |
| **Total experimentos** | 48,600 | 48,600 |

**Status**: ✅ **PARIDADE COMPLETA ATINGIDA**


## 📝 Uso Programático

### Experimento Individual

```python
from framework_qiskit import ClassificadorVQCQiskit, carregar_datasets

datasets = carregar_datasets()

vqc = ClassificadorVQCQiskit(
    n_qubits=4,
    n_camadas=2,
    arquitetura='strongly_entangling',
    estrategia_init='matematico',
    tipo_ruido='phase_damping',
    nivel_ruido=0.005,
    n_epocas=10,
    seed=42,
    shots=1024
)

vqc.fit(datasets['moons']['X_train'], datasets['moons']['y_train'])
acc = vqc.score(datasets['moons']['X_test'], datasets['moons']['y_test'])

```text

### Grid Search Customizado

```python
from executar_grid_search_qiskit import executar_experimento_completo

# Executar subset customizado
df = executar_experimento_completo(
    modo='medio',
    max_experimentos=100,
    pasta_resultados='resultados_custom',
    verbose=True
)

# Analisar resultados
print(df.groupby('tipo_ruido')['acuracia_teste'].mean())

```

## 🔗 Referências

- Framework PennyLane: `framework_investigativo_completo.py`
- Documentação Qiskit: `docs/GUIA_QISKIT.md`
- Comparação: `docs/COMPARACAO_PENNYLANE_QISKIT.md`


---


**Última Atualização**: 24/12/2025  
**Versão**: Framework Qiskit v7.2  
**Status**: Paridade completa com PennyLane ✅

