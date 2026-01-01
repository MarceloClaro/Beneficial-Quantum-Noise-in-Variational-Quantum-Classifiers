# Framework QAOA para 100 Qubits com Análise de Ruído Benéfico

## 📋 Visão Geral

Este módulo implementa o **QAOA (Quantum Approximate Optimization Algorithm)** escalável para **100 qubits** usando **Qiskit**, mantendo a metodologia de análise de ruído quântico benéfico do projeto original.

### Características Principais

- ✅ **Escalabilidade**: Suporta até 100 qubits
- ✅ **QAOA Completo**: Implementação full-stack com mixing e problem Hamiltonians
- ✅ **Análise de Ruído**: 4 tipos de ruído quântico (depolarizing, amplitude damping, phase damping, thermal)
- ✅ **Otimização**: Loop quântico-clássico com COBYLA, SPSA, etc.
- ✅ **Busca de Hiperparâmetros**: Grid search e otimização Bayesiana (Optuna)
- ✅ **Visualizações**: Gráficos interativos de convergência e comparação
- ✅ **Reprodutibilidade**: Seeds fixas, logging completo, resultados salvos


---


## 🚀 Início Rápido

### Instalação

```bash

# Instalar dependências
pip install qiskit qiskit-aer numpy pandas scipy matplotlib plotly optuna

# Ou usar requirements.txt existente
pip install -r requirements.txt

```text

### Execução Rápida

```bash

# Demonstração rápida (20 qubits, ~2 minutos)
python executar_qaoa_100qubits.py rapido

# Grid search (30 qubits, ~15 minutos)
python executar_qaoa_100qubits.py grid

# Teste de níveis de ruído (25 qubits, ~10 minutos)
python executar_qaoa_100qubits.py niveis

# Experimento completo 100 qubits (LONGO - várias horas)
python executar_qaoa_100qubits.py completo

```text

### Uso Programático

```python
from framework_qaoa_100qubits import (
    ConfigQAOA,
    ConstrutorCircuitoQAOA,
    OtimizadorQAOA,
    demo_qaoa_100qubits
)

# Demo rápida
resultado = demo_qaoa_100qubits(
    densidade_grafo=0.1,
    p_layers=3,
    tipo_ruido='depolarizing',
    nivel_ruido=0.001
)

print(f"Energia final: {resultado.energia_final:.4f}")
print(f"Tempo de execução: {resultado.tempo_execucao:.2f}s")

```text

---


## 🏗️ Arquitetura

### Módulos Principais

```

framework_qaoa_100qubits.py
├── ConfigQAOA                    # Dataclass de configuração
├── ConstrutorCircuitoQAOA        # Construção de circuitos QAOA
├── ModeloRuidoQAOA               # Modelos de ruído quântico
├── OtimizadorQAOA                # Loop quântico-clássico
├── AnalisadorHiperparametrosQAOA # Grid search e Bayesian opt
└── VisualizadorQAOA              # Gráficos e visualizações

```text

### Fluxo de Trabalho

```

1. Criar Grafo

   ├─> Problema MaxCut (matriz de adjacência)
   └─> Densidade configurável (0.0-1.0)

2. Configurar QAOA

   ├─> Número de qubits (1-100)
   ├─> Profundidade p (1-10 camadas)
   ├─> Tipo de ruído (depolarizing, amplitude, phase, thermal)
   └─> Nível de ruído (0.0-0.05)

3. Otimizar

   ├─> Inicializar parâmetros γ, β
   ├─> Loop quântico-clássico
   │   ├─> Executar circuito
   │   ├─> Calcular energia
   │   └─> Atualizar parâmetros
   └─> Convergência

4. Analisar Resultados

   ├─> Energia final
   ├─> Convergência
   ├─> Probabilidades
   └─> Comparação com/sem ruído

```text

---


## 📊 Exemplo: Análise de Ruído Benéfico

### Código

```python
from framework_qaoa_100qubits import (
    ConfigQAOA, ConstrutorCircuitoQAOA, OtimizadorQAOA
)
import numpy as np

# 1. Criar grafo MaxCut (50 qubits)
construtor = ConstrutorCircuitoQAOA(n_qubits=50, p_layers=3)
grafo = construtor.criar_grafo_aleatorio(densidade=0.15)

# 2. QAOA sem ruído (baseline)
config_baseline = ConfigQAOA(
    n_qubits=50,
    p_layers=3,
    tipo_ruido='sem_ruido',
    nivel_ruido=0.0,
    max_iter=100
)
otimizador_baseline = OtimizadorQAOA(config_baseline)
resultado_baseline = otimizador_baseline.otimizar(grafo)

# 3. QAOA com ruído depolarizante
config_ruido = ConfigQAOA(
    n_qubits=50,
    p_layers=3,
    tipo_ruido='depolarizing',
    nivel_ruido=0.001,
    max_iter=100
)
otimizador_ruido = OtimizadorQAOA(config_ruido)
resultado_ruido = otimizador_ruido.otimizar(grafo)

# 4. Comparação
print(f"Energia sem ruído:  {resultado_baseline.energia_final:.4f}")
print(f"Energia com ruído:  {resultado_ruido.energia_final:.4f}")

melhoria = (resultado_baseline.energia_final - resultado_ruido.energia_final) / resultado_baseline.energia_final
print(f"Melhoria relativa:  {melhoria*100:+.2f}%")

if melhoria > 0:
    print("✅ RUÍDO BENÉFICO DETECTADO!")
else:
    print("⚠️  Ruído prejudicou neste caso")

```text

### Saída Esperada

```

Energia sem ruído:  45.6782
Energia com ruído:  44.2315
Melhoria relativa:  +3.17%
✅ RUÍDO BENÉFICO DETECTADO!

```text

---


## 🔬 Tipos de Ruído Implementados

### 1. Depolarizing Noise

**Canal de Lindblad**: ρ → (1-p)ρ + p·I/d


- Simula erro genérico em todas as direções
- Taxa de erro típica: 0.001-0.01
- Implementação: `NoiseModel.depolarizing_error()`


```python
config = ConfigQAOA(tipo_ruido='depolarizing', nivel_ruido=0.001)

```text

### 2. Amplitude Damping

**Canal T1**: Simula perda de energia (decay |1⟩ → |0⟩)


- Modela relaxação de amplitude
- Taxa típica: 0.0005-0.005
- Implementação: `NoiseModel.amplitude_damping_error()`


```python
config = ConfigQAOA(tipo_ruido='amplitude_damping', nivel_ruido=0.002)

```text

### 3. Phase Damping

**Canal T2**: Simula perda de coerência (dephasing)


- Modela decoerência de fase
- Taxa típica: 0.001-0.01
- Implementação: `NoiseModel.phase_damping_error()`


```python
config = ConfigQAOA(tipo_ruido='phase_damping', nivel_ruido=0.001)

```text

### 4. Thermal Relaxation

**Modelo realista**: Combina T1 e T2


- T1: 50 μs (amplitude)
- T2: 70 μs (coerência, T2 ≤ 2·T1)
- Tempo de porta: 100 ns


```python
config = ConfigQAOA(tipo_ruido='thermal')

```text

---


## 🎯 Grid Search de Hiperparâmetros

### Exemplo Completo

```python
from framework_qaoa_100qubits import (
    ConstrutorCircuitoQAOA,
    AnalisadorHiperparametrosQAOA
)

# 1. Criar problema
construtor = ConstrutorCircuitoQAOA(n_qubits=40, p_layers=3)
grafo = construtor.criar_grafo_aleatorio(densidade=0.2)

# 2. Configurar grid search
analisador = AnalisadorHiperparametrosQAOA(
    pasta_resultados='resultados_qaoa_grid'
)

# 3. Executar grid search
df_resultados = analisador.grid_search_ruido(
    grafo=grafo,
    niveis_ruido=[0.0, 0.0001, 0.0005, 0.001, 0.002, 0.005],
    tipos_ruido=['sem_ruido', 'depolarizing', 'amplitude_damping', 'phase_damping'],
    p_layers=3,
    n_repeticoes=5
)

# 4. Analisar
print(df_resultados.groupby(['tipo_ruido', 'nivel_ruido'])['energia_final'].agg(['mean', 'std']))

# 5. Visualizar
from framework_qaoa_100qubits import VisualizadorQAOA
visualizador = VisualizadorQAOA()
visualizador.plotar_comparacao_ruido(
    df_resultados,
    salvar='comparacao_ruido.html'
)

```text

### Resultados Esperados (CSV)

```csv
tipo_ruido,nivel_ruido,p_layers,repeticao,energia_final,tempo_execucao,iteracoes
sem_ruido,0.0000,3,0,45.67,12.3,48
depolarizing,0.0001,3,0,45.23,12.8,52
depolarizing,0.0005,3,0,44.81,13.1,49
depolarizing,0.0010,3,0,44.15,13.4,51
...

```text

---


## 📈 Otimização Bayesiana (Optuna)

### Busca Automática de Hiperparâmetros

```python
from framework_qaoa_100qubits import AnalisadorHiperparametrosQAOA

analisador = AnalisadorHiperparametrosQAOA()

# Otimização Bayesiana (50 trials)
resultado_bayes = analisador.otimizacao_bayesiana(
    grafo=grafo,
    n_trials=50
)

print("Melhores hiperparâmetros:")
print(f"  Tipo ruído:  {resultado_bayes['best_params']['tipo_ruido']}")
print(f"  Nível ruído: {resultado_bayes['best_params']['nivel_ruido']:.4f}")
print(f"  P-layers:    {resultado_bayes['best_params']['p_layers']}")
print(f"\nMelhor energia: {resultado_bayes['best_value']:.4f}")

```text

### Espaço de Busca

- **tipo_ruido**: [depolarizing, amplitude_damping, phase_damping, sem_ruido]
- **nivel_ruido**: [0.0001, 0.01] (log-scale)
- **p_layers**: [1, 5]


---


## 🔧 Configuração Avançada

### Parâmetros de ConfigQAOA

```python
from framework_qaoa_100qubits import ConfigQAOA

config = ConfigQAOA(
    n_qubits=100,           # Número de qubits (1-100)
    p_layers=3,             # Profundidade QAOA (1-10)
    tipo_ruido='depolarizing',  # Tipo de ruído
    nivel_ruido=0.001,      # Taxa de erro (0.0-0.05)
    shots=1024,             # Medições por execução
    max_iter=100,           # Iterações máximas
    seed=42,                # Semente aleatória
    problema='maxcut',      # Tipo de problema
    otimizador='COBYLA'     # Otimizador clássico
)

```text

### Otimizadores Disponíveis

| Otimizador | Descrição | Uso |
|------------|-----------|-----|
| **COBYLA** | Constrained Optimization BY Linear Approximations | Padrão, robusto |
| **SLSQP** | Sequential Least Squares Programming | Gradiente numérico |
| **Powell** | Powell's method | Sem gradiente |
| **Nelder-Mead** | Simplex method | Sem gradiente |
| **L-BFGS-B** | Limited-memory BFGS | Com bounds |

---


## 📊 Visualizações

### 1. Convergência

```python
from framework_qaoa_100qubits import VisualizadorQAOA

visualizador = VisualizadorQAOA()
visualizador.plotar_convergencia(
    resultado,
    salvar='convergencia.html'
)

```text

**Exibe**: Energia vs. Iteração, mostrando trajetória de otimização


### 2. Comparação de Ruído

```python
visualizador.plotar_comparacao_ruido(
    df_resultados,
    salvar='comparacao_ruido.html'
)

```text

**Exibe**: Energia vs. Nível de Ruído para cada tipo, com barras de erro


---


## 🧪 Casos de Uso

### 1. Pesquisa Acadêmica

**Objetivo**: Investigar regime de ruído benéfico em QAOA


```python

# Experimento controlado com múltiplas repetições
resultados = experimento_completo_ruido_benefico(
    n_qubits=80,
    densidade_grafo=0.1,
    p_layers=4
)

# Análise estatística (ANOVA, effect sizes)
# Publicação em periódico QUALIS A1

```text

### 2. Benchmarking

**Objetivo**: Comparar diferentes configurações de hardware simulado


```python

# Testar impacto de T1/T2 no desempenho QAOA
for T1 in [20e3, 50e3, 100e3]:  # ns
    for T2 in [30e3, 70e3, 150e3]:
        config = ConfigQAOA(tipo_ruido='thermal', ...)

        # Executar e comparar

```text

### 3. Otimização de Hiperparâmetros

**Objetivo**: Encontrar configuração ótima para problema específico


```python

# Busca Bayesiana automática
resultado_bayes = analisador.otimizacao_bayesiana(
    grafo=meu_problema,
    n_trials=100
)

# Aplicar melhores parâmetros em produção

```text

---


## 📐 Fundamentos QAOA

### Formulação Matemática

**Objetivo**: Minimizar função de custo C(x) = Σ_{(i,j)} w_{ij}(1 - Z_i Z_j)/2


**Ansatz QAOA**:

```

|ψ(γ,β)⟩ = U(B,β_p) U(C,γ_p) ... U(B,β_1) U(C,γ_1) |+⟩^⊗n

```text

Onde:

- **U(C,γ)** = e^{-iγC}: Hamiltoniano do problema
- **U(B,β)** = e^{-iβB}: Hamiltoniano de mixing (B = Σ_i X_i)
- **γ, β**: Parâmetros variacionais (comprimento p)


### Implementação

```python
def criar_circuito_maxcut(grafo, gammas, betas):
    qc = QuantumCircuit(n_qubits)
    
    # Estado inicial |+⟩^⊗n
    qc.h(range(n_qubits))
    
    for p in range(p_layers):

        # Problem Hamiltonian: ZZ entre arestas
        for (i, j) in arestas:
            qc.cx(i, j)
            qc.rz(2 * gammas[p] * w_ij, j)
            qc.cx(i, j)
        
        # Mixing Hamiltonian: RX em todos
        for i in range(n_qubits):
            qc.rx(2 * betas[p], i)
    
    qc.measure_all()
    return qc

```text

---


## 🔬 Metodologia Científica

### Reprodutibilidade

- ✅ **Seeds fixas**: `np.random.seed(42)`, `seed` em ConfigQAOA
- ✅ **Logging completo**: Timestamps, parâmetros, resultados
- ✅ **Versionamento**: Git, releases, DOI (Zenodo)
- ✅ **Ambiente**: `requirements.txt`, Python 3.9+


### Validação Estatística

```python

# Múltiplas repetições para análise estatística
n_repeticoes = 10
resultados = []

for rep in range(n_repeticoes):
    config = ConfigQAOA(seed=42 + rep, ...)
    resultado = otimizador.otimizar(grafo)
    resultados.append(resultado.energia_final)

# ANOVA, t-test, effect sizes
from scipy.stats import ttest_ind
stat, p_value = ttest_ind(energias_sem_ruido, energias_com_ruido)
print(f"p-value: {p_value:.4f}")

```text

---


## 📚 Referências

### QAOA

1. **Farhi et al. (2014)**. "A Quantum Approximate Optimization Algorithm." arXiv:1411.4028
2. **Zhou et al. (2020)**. "Quantum approximate optimization algorithm: Performance, mechanism, and implementation." PRX Quantum, 1(2), 020319
3. **Crooks (2018)**. "Performance of the QAOA on p-spin models." arXiv:1811.08419


### Ruído Quântico

4. **Marshall et al. (2020)**. "Characterizing local noise in QAOA circuits." Quantum Sci. Technol., 5(1), 015005
5. **Wang et al. (2021)**. "Noise-induced barren plateaus in variational quantum algorithms." Nature Commun., 12(1), 6961
6. **Preskill (2018)**. "Quantum Computing in the NISQ era and beyond." Quantum, 2, 79


### Implementações

7. **Qiskit Documentation**: <https://qiskit.org/documentation/>
8. **PennyLane QAOA**: <https://pennylane.ai/qml/demos/tutorial_qaoa_intro.html>
9. **Cirq QAOA**: <https://quantumai.google/cirq/tutorials/qaoa>


---


## 🤝 Contribuindo

Contribuições são bem-vindas! Por favor:

1. Fork o repositório
2. Crie branch para feature: `git checkout -b feature/nova-funcionalidade`
3. Commit: `git commit -m 'Adiciona nova funcionalidade'`
4. Push: `git push origin feature/nova-funcionalidade`
5. Abra Pull Request


### Áreas de Contribuição

- 🔧 Novos tipos de ruído
- 📊 Visualizações adicionais
- 🚀 Otimizações de performance
- 📖 Melhorias na documentação
- 🧪 Casos de teste
- 🌐 Suporte para hardware real (IBM Quantum)


---


## 📝 Licença

MIT License - veja arquivo LICENSE para detalhes

---


## 📧 Contato

- **Issues**: [GitHub Issues](https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers/issues)
- **Discussions**: [GitHub Discussions](https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers/discussions)


---


## 🏆 Citação

Se usar este framework em pesquisa acadêmica, por favor cite:

```bibtex
@software{framework_qaoa_100qubits,
  title = {Framework QAOA para 100 Qubits com Análise de Ruído Benéfico},
  author = {Projeto Beneficial Quantum Noise in VQC},
  year = {2025},
  url = {<https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers},>
  version = {1.0.0}
}

```

---


## ✨ Agradecimentos

- **Qiskit Team** (IBM Quantum) - Framework quântico
- **Optuna Team** - Otimização Bayesiana
- **Projeto original** - Metodologia de ruído benéfico


---


**Status**: ✅ Produção | 🔬 Validado Cientificamente | 📊 Reprodutível


**Última atualização**: 2025-12-26

