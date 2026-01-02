# Framework Quantum Advanced V8 - Documentação Completa

## 📋 Visão Geral

O **Framework Quantum Advanced V8** é um framework robusto e modular para otimização variacional quântica moderna, desenvolvido para resolver limitações dos frameworks anteriores e atender aos padrões de rigor científico QUALIS A1.

### Principais Características

✅ **Otimização Variacional Quântica Moderna**
- VQE (Variational Quantum Eigensolver) híbrido
- QAOA (Quantum Approximate Optimization Algorithm)
- Otimização com múltiplos métodos (ADAM, SPSA, COBYLA, etc.)

✅ **Suporte Multi-Framework**
- PennyLane (principal)
- Qiskit (em desenvolvimento)
- Cirq (em desenvolvimento)

✅ **Quantum Error Mitigation (QEM)**
- Zero-Noise Extrapolation (ZNE) com múltiplas estratégias
- TREX (Twirling-based Error Extraction)
- AUEC (Adaptive Unified Error Correction)
- Readout Mitigation

✅ **Análise Científica Rigorosa**
- Análise de complexidade computacional quântica
- Validação de fórmulas de predição de ruído
- Benchmarks contra algoritmos estado-da-arte
- Métricas de escalamento

✅ **Datasets Moleculares**
- HIV Activity Dataset (DeepChem)
- Malária Dataset (DeepChem)
- Tuberculose Dataset (DeepChem)
- Datasets clássicos (Iris, Wine, Breast Cancer)

✅ **Tracing e Telemetria**
- Logging científico QUALIS A1
- Histórico de treinamento detalhado
- Visualizações profissionais

---

## 🚀 Instalação

### 1. Dependências Básicas

```bash
pip install numpy scipy scikit-learn pandas matplotlib seaborn plotly
```

### 2. Quantum Frameworks

```bash
# PennyLane (recomendado)
pip install pennylane

# Qiskit (opcional)
pip install qiskit qiskit-aer qiskit-ibm-runtime

# Cirq (opcional)
pip install cirq cirq-google
```

### 3. DeepChem (para datasets moleculares)

```bash
# Clone e instale DeepChem
git clone https://github.com/deepchem/deepchem.git
cd deepchem
source scripts/install_deepchem_conda.sh 3.10 cpu
```

Ou via pip:
```bash
pip install deepchem
```

### 4. Otimização Bayesiana (opcional)

```bash
pip install optuna
```

---

## 📖 Como Usar

### Uso Básico

```python
from framework_quantum_advanced_v8 import (
    ExperimentConfig,
    QuantumCircuitConfig,
    NoiseConfig,
    OptimizationConfig,
    ErrorMitigationConfig,
    QuantumExperimentRunner,
    FrameworkType,
    NoiseModel,
    OptimizationMethod,
    ErrorMitigationTechnique,
)

# 1. Configurar circuito quântico
circuit_config = QuantumCircuitConfig(
    n_qubits=4,
    n_layers=2,
    n_parameters=24,  # n_qubits * n_layers * 3
    architecture="strongly_entangling",
    entanglement="full"
)

# 2. Configurar ruído
noise_config = NoiseConfig(
    noise_model=NoiseModel.DEPOLARIZING,
    noise_level=0.01,  # 1% de probabilidade de erro
    gate_error=0.001,
    readout_error=0.01
)

# 3. Configurar otimização
opt_config = OptimizationConfig(
    method=OptimizationMethod.ADAM,
    learning_rate=0.1,
    max_iterations=100,
    early_stopping=True
)

# 4. Configurar mitigação de erro
zne_config = ErrorMitigationConfig(
    technique=ErrorMitigationTechnique.ZNE,
    zne_scale_factors=[1.0, 1.5, 2.0, 2.5],
    zne_extrapolation_type="exponential"
)

# 5. Configuração completa
config = ExperimentConfig(
    framework=FrameworkType.PENNYLANE,
    circuit_config=circuit_config,
    noise_config=noise_config,
    optimization_config=opt_config,
    error_mitigation_config=zne_config,
    dataset_name="iris",
    n_shots=1024,
    seed=42
)

# 6. Executar experimento
runner = QuantumExperimentRunner(config, results_dir="./results_quantum")
results = runner.run_full_experiment()

# 7. Salvar resultados
runner.save_results("results.json")
runner.save_plots()
```

### Via Linha de Comando

```bash
# Uso padrão (Iris, 4 qubits)
python run_framework_quantum_advanced_v8.py

# Com dataset HIV
python run_framework_quantum_advanced_v8.py --dataset hiv --n_qubits 6 --n_layers 3

# Com ruído customizado
python run_framework_quantum_advanced_v8.py --noise_level 0.05 --noise_model depolarizing

# Com mitigação ZNE
python run_framework_quantum_advanced_v8.py --error_mitigation zne

# Com QAOA
python run_framework_quantum_advanced_v8.py --n_layers 5 --max_iterations 200

# Combinação completa
python run_framework_quantum_advanced_v8.py \
    --dataset malaria \
    --n_qubits 8 \
    --n_layers 4 \
    --noise_level 0.02 \
    --learning_rate 0.05 \
    --max_iterations 150 \
    --error_mitigation zne \
    --results_dir ./results_advanced
```

---

## 🔧 Configurações Detalhadas

### QuantumCircuitConfig

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `n_qubits` | int | 4 | Número de qubits |
| `n_layers` | int | 2 | Número de camadas do circuito |
| `n_parameters` | int | 24 | Total de parâmetros |
| `architecture` | str | "strongly_entangling" | Arquitetura do circuito |
| `entanglement` | str | "full" | Tipo de emaranhamento ("full", "linear") |
| `skip_final_rotation` | bool | False | Pular rotação final |

### NoiseConfig

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `noise_model` | NoiseModel | DEPOLARIZING | Modelo de ruído |
| `noise_level` | float | 0.01 | Probabilidade de erro |
| `gate_error` | float | 0.001 | Taxa de erro de gate |
| `readout_error` | float | 0.01 | Taxa de erro de readout |
| `t1_time` | float | 100.0 | Tempo T1 em μs |
| `t2_time` | float | 50.0 | Tempo T2 em μs |

### OptimizationConfig

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `method` | OptimizationMethod | ADAM | Método de otimização |
| `learning_rate` | float | 0.1 | Taxa de aprendizado |
| `max_iterations` | int | 500 | Número máximo de iterações |
| `tolerance` | float | 1e-6 | Tolerância de convergência |
| `early_stopping` | bool | True | Usar early stopping |
| `early_stopping_patience` | int | 50 | Paciência para early stopping |

### ErrorMitigationConfig

| Parâmetro | Tipo | Padrão | Descrição |
|-----------|------|--------|-----------|
| `technique` | ErrorMitigationTechnique | ZNE | Técnica de mitigação |
| `zne_scale_factors` | List[float] | [1.0, 1.5, 2.0, 2.5] | Fatores de escala ZNE |
| `zne_extrapolation_type` | str | "exponential" | Tipo de extrapolação |
| `trex_num_twirls` | int | 10 | Número de twirls para TREX |
| `auec_num_layers` | int | 5 | Número de camadas AUEC |

---

## 📊 Análise de Complexidade Quântica

O framework inclui análise rigorosa de complexidade:

```python
from framework_quantum_advanced_v8 import QuantumComplexityAnalyzer

analyzer = QuantumComplexityAnalyzer()

# Calcular profundidade do circuito
depth = analyzer.calculate_circuit_depth(
    n_qubits=4,
    n_layers=2,
    entanglement="full"
)
print(f"Circuit depth: {depth}")

# Contar gates
gates = analyzer.calculate_gate_count(
    n_qubits=4,
    n_layers=2,
    entanglement="full"
)
print(f"Single qubit gates: {gates['single_qubit']}")
print(f"Two qubit gates: {gates['two_qubit']}")

# Estimar probabilidade de barren plateau
barren_prob = analyzer.estimate_barren_plateau_probability(
    n_qubits=4,
    circuit_depth=depth,
    noise_level=0.01
)
print(f"Barren plateau probability: {barren_prob:.4f}")

# Análise completa de recursos
resources = analyzer.analyze_resource_requirements(circuit_config, n_shots=1024)
print(f"Estimated simulation time: {resources['total_estimated_time_s']:.2f}s")
```

---

## 🔍 Validação de Ruído

Framework rigoroso para validar fórmulas de predição de ruído:

```python
from framework_quantum_advanced_v8 import NoiseValidationFramework

validator = NoiseValidationFramework()

# Predizer impacto de ruído
predicted_fidelity = validator.predict_noise_impact(
    noise_level=0.01,
    circuit_depth=10,
    n_qubits=4
)
print(f"Predicted fidelity: {predicted_fidelity:.4f}")

# Validar contra dados reais
validation = validator.validate_noise_prediction(
    actual_fidelity=0.85,
    predicted_fidelity=predicted_fidelity
)
print(f"Validation passed: {validation['validation_passed']}")
print(f"Relative error: {validation['relative_error']:.4f}")

# Analisar escalamento de ruído
scaling_analysis = validator.analyze_noise_scaling(
    noise_levels=[0.001, 0.005, 0.01, 0.02],
    circuit_depths=[5, 10, 15, 20],
    measured_errors={
        (0.001, 5): 0.02,
        (0.001, 10): 0.05,
        (0.01, 5): 0.08,
        (0.01, 10): 0.15,
    }
)
print(f"Noise scaling model: {scaling_analysis['model']}")
```

---

## 🎯 Benchmarking

Comparar contra algoritmos estado-da-arte:

```python
from framework_quantum_advanced_v8 import QuantumAlgorithmBenchmark

benchmark = QuantumAlgorithmBenchmark()

# Comparar VQC vs algoritmo clássico
vqc_predictions = np.array([0.8, 0.2, 0.9, 0.1])
classical_predictions = np.array([0.7, 0.3, 0.85, 0.15])
y_true = np.array([1, 0, 1, 0])

comparison = benchmark.benchmark_against_classical(
    vqc_predictions,
    classical_predictions,
    y_true
)
print(f"VQC Accuracy: {comparison['vqc_metrics']['accuracy']:.4f}")
print(f"Classical Accuracy: {comparison['classical_metrics']['accuracy']:.4f}")
print(f"VQC Wins: {comparison['vqc_wins']}")

# Analisar escalamento
scaling = benchmark.benchmark_scaling(
    system_sizes=[2, 3, 4, 5],
    execution_times=[0.1, 0.5, 2.0, 8.0]
)
print(f"Scaling exponent: {scaling['scaling_exponent']:.2f}")
print(f"Scaling type: {scaling['scaling_type']}")
```

---

## 💾 Datasets DeepChem

Carregar datasets moleculares:

```python
from framework_quantum_advanced_v8 import DeepChemDatasetLoader

loader = DeepChemDatasetLoader()

# Carregar HIV
X_hiv, y_hiv = loader.load_dataset("hiv")
print(f"HIV dataset: X shape = {X_hiv.shape}")

# Carregar Malária
X_malaria, y_malaria = loader.load_dataset("malaria")
print(f"Malaria dataset: X shape = {X_malaria.shape}")

# Carregar Tuberculose
X_tb, y_tb = loader.load_dataset("tb")
print(f"TB dataset: X shape = {X_tb.shape}")
```

### Datasets Disponíveis

| Dataset | Amostras | Features | Descrição |
|---------|----------|----------|-----------|
| Iris | 150 | 4 | Classificação de flores |
| Wine | 178 | 13 | Análise de vinho |
| Breast Cancer | 569 | 30 | Diagnóstico de câncer |
| **HIV** | ~8000 | 1024 | Atividade contra HIV |
| **Malaria** | ~10000 | 1024 | Atividade contra Malária |
| **TB** | ~8000 | 1024 | Atividade contra TB |

---

## 📈 Estrutura de Resultados

Os resultados são salvos em JSON com a seguinte estrutura:

```json
{
  "complexity_analysis": {
    "circuit_depth": 12,
    "gate_count": {
      "single_qubit": 24,
      "two_qubit": 12,
      "total": 36,
      "parametrized": 24
    },
    "barren_plateau_probability": 0.15,
    "entanglement_entropy": 1.8,
    "classical_simulation_overhead": 16,
    "estimated_simulation_time_per_shot_ms": 0.4,
    "total_estimated_time_s": 0.4
  },
  "data_shapes": {
    "train": [100, 4],
    "val": [27, 4],
    "test": [50, 4]
  },
  "training_metrics": {
    "final_train_loss": 0.12,
    "best_val_loss": 0.15,
    "final_train_accuracy": 0.92,
    "best_val_accuracy": 0.88,
    "convergence_iterations": 100
  },
  "inference_metrics": {
    "accuracy": 0.86,
    "precision": 0.85,
    "recall": 0.87,
    "f1": 0.86,
    "auc": 0.92
  },
  "noise_analysis": {
    "predicted_fidelity": 0.956
  },
  "execution_time_seconds": 45.2
}
```

---

## 🔐 Validação Científica

O framework implementa rigorosos critérios de validação:

### 1. Validação de Convergência
- Monitoramento contínuo de loss
- Early stopping com paciência configurável
- Histórico completo de treinamento

### 2. Validação de Complexidade
- Cálculo de profundidade de circuito
- Contagem de gates
- Estimativa de probabilidade de barren plateau

### 3. Validação de Ruído
- Predição de fidelidade
- Comparação com dados empíricos
- Análise de escalamento

### 4. Validação de Performance
- Métricas de classificação (accuracy, precision, recall, F1, AUC)
- Comparação com baseline clássico
- Análise de escalamento de recursos

---

## 🧪 Exemplos Práticos

### Exemplo 1: Experimento Simples com Iris

```python
from framework_quantum_advanced_v8 import (
    ExperimentConfig,
    QuantumCircuitConfig,
    NoiseConfig,
    OptimizationConfig,
    ErrorMitigationConfig,
    QuantumExperimentRunner,
    FrameworkType,
    NoiseModel,
    OptimizationMethod,
    ErrorMitigationTechnique,
)

config = ExperimentConfig(
    framework=FrameworkType.PENNYLANE,
    circuit_config=QuantumCircuitConfig(
        n_qubits=4,
        n_layers=2,
        n_parameters=24
    ),
    noise_config=NoiseConfig(
        noise_model=NoiseModel.DEPOLARIZING,
        noise_level=0.01
    ),
    optimization_config=OptimizationConfig(
        method=OptimizationMethod.ADAM,
        learning_rate=0.1,
        max_iterations=100
    ),
    error_mitigation_config=ErrorMitigationConfig(
        technique=ErrorMitigationTechnique.ZNE
    ),
    dataset_name="iris"
)

runner = QuantumExperimentRunner(config)
results = runner.run_full_experiment()
runner.save_results()
runner.save_plots()
```

### Exemplo 2: Experimento com Dataset HIV

```bash
python run_framework_quantum_advanced_v8.py \
    --dataset hiv \
    --n_qubits 6 \
    --n_layers 3 \
    --noise_level 0.02 \
    --learning_rate 0.05 \
    --max_iterations 150 \
    --error_mitigation zne \
    --results_dir ./results_hiv_advanced
```

### Exemplo 3: Experimento com Mitigação Avançada

```python
config = ExperimentConfig(
    framework=FrameworkType.PENNYLANE,
    circuit_config=QuantumCircuitConfig(
        n_qubits=5,
        n_layers=4,
        n_parameters=60
    ),
    noise_config=NoiseConfig(
        noise_model=NoiseModel.AMPLITUDE_DAMPING,
        noise_level=0.05,
        t1_time=100.0,
        t2_time=50.0
    ),
    optimization_config=OptimizationConfig(
        method=OptimizationMethod.COBYLA,
        max_iterations=200,
        early_stopping=True,
        early_stopping_patience=30
    ),
    error_mitigation_config=ErrorMitigationConfig(
        technique=ErrorMitigationTechnique.ZNE,
        zne_scale_factors=[1.0, 1.5, 2.0, 2.5, 3.0],
        zne_extrapolation_type="polynomial"
    ),
    dataset_name="malaria",
    n_shots=2048
)
```

---

## 📝 Referências Bibliográficas

1. **Cerezo et al.** "Barren plateaus in quantum neural landscape design" Nature Reviews Physics (2021)

2. **Colless et al.** "Computation of Molecular Spectra on a Quantum Processor" Physical Review Letters (2018)

3. **Wang et al.** "Noise-Induced Barren Plateaus" Nature Communications (2021)

4. **Kandala et al.** "Hardware-efficient Variational Quantum Eigensolver for Small Molecules and Quantum Magnets" Nature (2017)

5. **Farhi et al.** "A Quantum Approximate Optimization Algorithm" arXiv:1411.4028 (2014)

---

## 🤝 Contribuições

Este framework foi desenvolvido como parte de pesquisa em computação quântica e classificadores variacionais quânticos para aplicações moleculares.

## 📄 Licença

[Especificar licença do projeto]

## 👨‍💻 Autores

Framework Científico QUALIS A1 - 2026
