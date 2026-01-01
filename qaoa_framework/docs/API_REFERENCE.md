# Referência da API - QAOA Framework

## Guia de Início Rápido

```python
from qaoa_framework.scripts import problem_generator, circuit_builder

# Gerar problema
problem = problem_generator.gerar_problema_benchmark(
    graph_type='erdos_renyi',
    n_nodes=20,
    seed=42
)

# Construir circuito
builder = circuit_builder.QAOACircuitBuilder(
    n_qubits=20,
    p_layers=3,
    hamiltonian=problem.to_hamiltonian()
)

qc = builder.build()
print(f"Circuit depth: {qc.depth()}")

```text

## Módulos Principais

1. **problem_generator**: Geração de problemas MaxCut
2. **circuit_builder**: Construção de circuitos QAOA
3. **noise_models**: Modelos de ruído quântico
4. **hyperparameter_tuning**: Otimização Bayesiana com Optuna
5. **visualization**: Visualizações interativas Plotly
6. **statistical_analysis**: Análise estatística rigorosa


## Consulte os Docstrings

Todos os módulos possuem docstrings detalhadas com exemplos de uso.

```python
from qaoa_framework.scripts import problem_generator
help(problem_generator.gerar_problema_benchmark)

```

**Versão**: 1.0  
**Autores**: Projeto Beneficial Quantum Noise

