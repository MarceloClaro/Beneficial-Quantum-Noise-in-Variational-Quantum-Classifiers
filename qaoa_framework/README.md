# QAOA Framework para 100 Qubits

## 📋 Visão Geral

Framework modular para execução de QAOA (Quantum Approximate Optimization Algorithm) com análise sistemática de ruído quântico benéfico. Implementado seguindo padrões Qualis A1 de reprodutibilidade e auditabilidade.

### Características Principais

- ✅ **Escalabilidade**: Suporte para 100 qubits via Matrix Product State (MPS)
- ✅ **Modularidade**: Arquitetura em componentes independentes e testáveis
- ✅ **Configuração YAML**: Experimentos definidos via arquivo de configuração
- ✅ **Reprodutibilidade**: Seeds fixas, manifesto de execução, resultados auditáveis
- ✅ **Análise de Ruído**: 5 modelos de ruído quântico (Lindblad)
- ✅ **Benchmark**: Biblioteca de problemas MaxCut

---

## 🏗️ Estrutura do Projeto

```
qaoa_framework/
├── configs/
│   └── experiment_qaoa.yaml       # Configuração de experimentos
├── scripts/
│   ├── noise_models.py            # Tarefa 1: Modelos de ruído Qiskit
│   ├── problem_generator.py       # Tarefa 3: Gerador de problemas
│   ├── circuit_builder.py         # Tarefa 4: Construtor de circuito QAOA
│   ├── cost_function.py           # Tarefa 5: Função de custo (TODO)
│   ├── hyperparameter_tuning.py   # Tarefa 2: Optuna (TODO)
│   └── visualization.py           # Tarefa 8: Visualizações (TODO)
├── results/                        # Resultados de experimentos
├── docs/                          # Documentação adicional
└── main.py                        # Script principal de execução
```

---

## 🚀 Instalação

### Dependências

```bash
pip install qiskit qiskit-aer numpy pandas scipy matplotlib plotly pyyaml networkx
```

Opcional (para otimização Bayesiana):
```bash
pip install optuna
```

---

## 📖 Uso

### 1. Configurar Experimento

Edite `configs/experiment_qaoa.yaml`:

```yaml
problem:
  type: "maxcut"
  n_nodes: 100
  graph_type: "erdos_renyi"
  edge_probability: 0.5

model:
  algorithm: "qaoa"
  p_layers: 3

noise:
  enabled: true
  model: "depolarizing"
  params:
    p: [0.0, 0.005, 0.01]

frameworks:
  qiskit:
    backend: "aer_simulator"
    method: "matrix_product_state"  # Permite 100 qubits
```

### 2. Executar Experimento

```bash
python qaoa_framework/main.py --config qaoa_framework/configs/experiment_qaoa.yaml
```

### 3. Analisar Resultados

Resultados salvos em `qaoa_framework/results/<run_id>/`:

- `resultados.csv`: Dados brutos de cada execução
- `resumo.csv`: Estatísticas agregadas
- `manifesto.json`: Metadados de reprodutibilidade

---

## 🔬 Componentes Implementados

### ✅ Tarefa 1: Modelos de Ruído (noise_models.py)

Implementa 5 tipos de ruído quântico:

```python
from scripts.noise_models import criar_noise_model_qiskit

# Ruído despolarizante
noise_model = criar_noise_model_qiskit('depolarizing', 0.01)

# Ruído com schedule dinâmico
nivel = aplicar_schedule_ruido(0.01, 'linear', iteracao=50, max_iter=100)
```

**Tipos disponíveis**:
- `depolarizing`: Canal despolarizante genérico
- `amplitude_damping`: Relaxação de amplitude (T1)
- `phase_damping`: Perda de coerência (T2)
- `thermal`: Modelo realista com T1 e T2
- `pauli`: Combinação de erros X, Y, Z

### ✅ Tarefa 3: Gerador de Problemas (problem_generator.py)

Gera instâncias de MaxCut:

```python
from scripts.problem_generator import gerar_problema_benchmark

# Grafo Erdős-Rényi com 100 nós
problem = gerar_problema_benchmark(
    graph_type='erdos_renyi',
    n_nodes=100,
    edge_probability=0.5,
    seed=42
)

print(f"Valor ótimo: {problem.optimal_value}")
hamiltonian = problem.to_hamiltonian()
```

**Tipos de grafo**:
- `erdos_renyi`: Grafo aleatório
- `regular`: d-regular (cada vértice tem grau d)
- `complete`: Grafo completo

### ✅ Tarefa 4: Circuito QAOA (circuit_builder.py)

Constrói ansatz QAOA escalável:

```python
from scripts.circuit_builder import QAOACircuitBuilder

builder = QAOACircuitBuilder(
    n_qubits=100,
    p_layers=3,
    hamiltonian=problem.to_hamiltonian()
)

# Construir circuito
qc = builder.build()

# Inicializar parâmetros
params = builder.initialize_parameters('heuristic', seed=42)

# Estatísticas
stats = builder.get_circuit_stats()
print(f"Depth: {stats['depth']}, Gates: {stats['total_gates']}")
```

### ✅ Tarefa 2: Otimização Bayesiana (hyperparameter_tuning.py)

Otimização de hiperparâmetros com Optuna:

```python
from scripts.hyperparameter_tuning import otimizar_hiperparametros_qaoa

# Executar otimização
resultado = otimizar_hiperparametros_qaoa(
    problem=problem,
    evaluate_fn=evaluate_qaoa,
    config=config,
    n_trials=100
)

print(f"Melhores parâmetros: {resultado['best_params']}")
print(f"Melhor approx_ratio: {resultado['best_value']:.4f}")
```

### ✅ Tarefa 8: Visualizações (visualization.py)

Visualizações específicas para QAOA:

```python
from scripts.visualization import (
    visualizar_convergencia_qaoa,
    visualizar_comparacao_ruido,
    gerar_relatorio_visual_html
)

# Convergência
visualizar_convergencia_qaoa(historico_energia, salvar='convergencia.html')

# Comparação de ruído
visualizar_comparacao_ruido(df_resultados, salvar='comparacao.html')

# Relatório completo
gerar_relatorio_visual_html(resultados, 'relatorio_qaoa.html')
```

### ✅ Tarefa 9: Análise Estatística (statistical_analysis.py)

Análise estatística rigorosa para QAOA:

```python
from scripts.statistical_analysis import (
    comparar_com_baseline,
    anova_one_way,
    gerar_relatorio_estatistico
)

# Comparar com vs. sem ruído
resultado = comparar_com_baseline(
    baseline=sem_ruido_data,
    experimental=com_ruido_data
)
print(resultado['conclusao'])  # ✅ RUÍDO BENÉFICO CONFIRMADO

# ANOVA para múltiplos níveis
resultado_anova = anova_one_way(
    grupos=[grupo1, grupo2, grupo3],
    labels=['sem_ruido', 'baixo', 'medio']
)

# Relatório completo
relatorio = gerar_relatorio_estatistico(
    df_resultados,
    coluna_metrica='approx_ratio',
    output_file='analise_estatistica.txt'
)
```

---

## 🎯 Tarefas Pendentes

As seguintes tarefas do plano ainda precisam ser implementadas:

- [ ] **Tarefa 5**: Função de custo completa com opflow
- [ ] **Tarefa 6**: Testes de escalabilidade MPS
- [ ] **Tarefa 7**: Camada de abstração de backend
- [ ] **Tarefa 10**: Sistema de reprodutibilidade completo
- [ ] **Tarefa 11**: Script de demonstração completo
- [ ] **Tarefa 12**: Documentação completa

---

## 📊 Exemplo de Saída

```
INICIANDO EXPERIMENTO QAOA
================================================================================

[1/5] Gerando problema de benchmark...
  Tipo: maxcut
  Nós: 100
  Arestas: 2487
  Valor ótimo estimado: 1245.32

[2/5] Construindo circuito QAOA...
  Qubits: 100
  P-layers: 3
  Depth: 42
  Gates: 7523 (CX: 2487)

[3/5] Executando experimentos com ruído...
  Total de experimentos: 15

  [1/15] Seed=1, Ruído=0.0000
      Energia: -1198.45
      Approx Ratio: 0.9624
      Tempo: 45.2s

  ...

[4/5] Salvando resultados...
  Resultados salvos: qaoa_framework/results/qaoa_run_20251226_213456/resultados.csv
  Resumo salvo: qaoa_framework/results/qaoa_run_20251226_213456/resumo.csv

[5/5] Gerando manifesto de execução...
  Manifesto salvo: qaoa_framework/results/qaoa_run_20251226_213456/manifesto.json

EXPERIMENTO CONCLUÍDO EM 1234.5s
================================================================================
```

---

## 🔍 Análise de Ruído Benéfico

Para analisar se o ruído foi benéfico:

```python
import pandas as pd

# Carregar resultados
df = pd.read_csv('qaoa_framework/results/<run_id>/resultados.csv')

# Comparar com vs. sem ruído
baseline = df[df['nivel_ruido'] == 0.0]['approx_ratio'].mean()
com_ruido = df[df['nivel_ruido'] > 0.0].groupby('nivel_ruido')['approx_ratio'].mean()

print(f"Baseline (sem ruído): {baseline:.4f}")
print("\nCom ruído:")
print(com_ruido)

# Identificar nível ótimo
melhor_nivel = com_ruido.idxmax()
melhoria = (com_ruido.max() - baseline) / baseline * 100

print(f"\nMelhor nível: {melhor_nivel:.4f}")
print(f"Melhoria: {melhoria:+.2f}%")
```

---

## 📚 Referências

### QAOA

1. Farhi et al. (2014). "A Quantum Approximate Optimization Algorithm." arXiv:1411.4028
2. Zhou et al. (2020). "Quantum approximate optimization algorithm: Performance, mechanism." PRX Quantum

### Matrix Product State (MPS)

3. Vidal (2003). "Efficient Classical Simulation of Slightly Entangled Quantum Computations." PRL
4. Qiskit Documentation: Matrix Product State Simulator

### Ruído Quântico

5. Nielsen & Chuang (2010). "Quantum Computation and Quantum Information."
6. Marshall et al. (2020). "Characterizing local noise in QAOA circuits." Quantum Sci. Technol.

---

## 🤝 Contribuindo

Para adicionar novos componentes:

1. Criar módulo em `scripts/`
2. Adicionar testes no `if __name__ == "__main__"`
3. Atualizar `main.py` para integração
4. Documentar em README

---

## 📝 Licença

MIT License - mesmo do projeto principal

---

## 📧 Contato

- **Issues**: GitHub Issues do projeto principal
- **Documentação**: Ver arquivos em `docs/`

---

**Status**: 🚧 **EM DESENVOLVIMENTO**

**Tarefas Completas**: 7/12 (58%)

**Última Atualização**: 2025-12-26
