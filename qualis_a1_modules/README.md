# Qualis A1 Modules - Scientific Rigor Enhancement

**Version:** 8.0-QAI  
**Purpose:** Enhance the "Beneficial Quantum Noise in VQC" framework to meet Qualis A1 publication standards (Nature, Quantum, Physical Review)


---


## 📋 Overview

This module package provides essential tools for:

1. **Mathematical Validation**: Verify quantum operator correctness
2. **Reproducibility**: Ensure experiment reproducibility
3. **Statistical Rigor**: Advanced statistical analysis with proper corrections
4. **Auditability**: Full code-to-method traceability
5. **Visualization**: Publication-quality circuit diagrams


---


## 📦 Module Structure

```text
qualis_a1_modules/
├── __init__.py                    # Package initialization
├── validation.py                   # Kraus operator validation
├── reproducibility.py              # Seed configuration and manifests
├── statistical_extensions.py       # Bonferroni, power analysis
├── auditing.py                     # Code→Method mapping
├── visualization.py                # Circuit diagram generation
└── README.md                       # This file

```

---


## 🔬 Module 1: Mathematical Validation (`validation.py`)

### Purpose
Validate mathematical correctness of quantum noise operators.

### Key Functions

#### `validar_operadores_kraus(operadores, tol=1e-10, verbose=True)`
Validates Kraus operator completeness: Σ_k K_k† K_k = I

**Example:**

```python
from qualis_a1_modules.validation import validar_operadores_kraus
import numpy as np

# Depolarizing channel operators (p=0.1)
p = 0.1
I = np.eye(2)
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]])

K0 = np.sqrt(1 - p) * I
K1 = np.sqrt(p/3) * X
K2 = np.sqrt(p/3) * Y
K3 = np.sqrt(p/3) * Z

# Validate
validar_operadores_kraus([K0, K1, K2, K3])  # Should pass

```text

#### Helper Functions:
- `obter_operadores_kraus_depolarizante(p)`: Get depolarizing channel operators
- `obter_operadores_kraus_amplitude_damping(gamma)`: Get amplitude damping operators
- `obter_operadores_kraus_phase_damping(lambda)`: Get phase damping operators


### Mathematical Foundation
All quantum channels must satisfy the trace-preserving condition:
$$\sum_k K_k^\dagger K_k = \mathbb{I}$$

This ensures that probability is conserved: Tr(ℰ(ρ)) = Tr(ρ) = 1.

### References
- Nielsen & Chuang (2010), Chapter 8: Quantum Operations
- Preskill (2018), "Quantum Computing in the NISQ era"


---


## 🔄 Module 2: Reproducibility (`reproducibility.py`)

### Purpose
Ensure complete reproducibility of experiments through seed management and execution manifests.

### Key Functions

#### `configurar_seeds_reprodutiveis(seed=42, verbose=True)`
Configures seeds for all randomness sources.

**Example:**

```python
from qualis_a1_modules.reproducibility import configurar_seeds_reprodutiveis

# Set seed for reproducibility
seed_info = configurar_seeds_reprodutiveis(seed=42, verbose=True)
print(f"Sources configured: {seed_info['configured_sources']}")

```text

#### Configured Sources:
- NumPy (`np.random.seed`)
- Python random (`random.seed`)
- PennyLane (via NumPy)
- Optuna (for hyperparameter optimization)
- scikit-learn (via NumPy)
- Qiskit (`algorithm_globals.random_seed`)


#### `gerar_manifesto_execucao(config, pasta_resultados, seed_info=None, parametros_experimento=None)`
Generates comprehensive execution manifest.

**Example:**

```python
from qualis_a1_modules.reproducibility import gerar_manifesto_execucao

config = {'version': '8.0-QAI', 'default_seed': 42}
manifesto_path = gerar_manifesto_execucao(
    config=config,
    pasta_resultados='./results',
    seed_info=seed_info,
    parametros_experimento={'n_trials': 100, 'noise_level': 0.01}
)

```text

#### Manifest Contains:
- Python version and platform
- All library versions (numpy, pennylane, qiskit, etc.)
- Git commit hash and branch
- Command line used
- Execution timestamp
- Random seeds used


### Compliance
Meets requirements of:

- Nature Scientific Data: "Code availability"
- PLoS Computational Biology: "Software and data availability"
- Quantum Journal: "Reproducibility statement"


### References
- Peng (2011), "Reproducible research in computational science"
- Sandve et al. (2013), "Ten simple rules for reproducible computational research"


---


## 📊 Module 3: Statistical Extensions (`statistical_extensions.py`)

### Purpose
Advanced statistical analysis with corrections for multiple comparisons and power analysis.

### Key Functions

#### `testes_post_hoc_com_correcao(df, grupo_col, metrica_col, metodo_correcao='bonferroni', alpha=0.05)`
Performs post-hoc tests with Bonferroni (or other) correction.

**Example:**

```python
from qualis_a1_modules.statistical_extensions import testes_post_hoc_com_correcao
import pandas as pd

df = pd.DataFrame({
    'ruido': ['depol']*10 + ['amplitude']*10 + ['phase']*10,
    'acuracia': [0.85, 0.87, 0.82, ...] # accuracies
})

resultados = testes_post_hoc_com_correcao(
    df, 'ruido', 'acuracia', metodo_correcao='bonferroni'
)
print(resultados[['grupo1', 'grupo2', 'p_valor', 'p_ajustado', 'significativo']])

```text

#### Correction Methods:
- `'bonferroni'`: α_adj = α/m (most conservative)
- `'holm'`: Holm-Bonferroni (sequential)
- `'fdr_bh'`: Benjamini-Hochberg FDR
- `'fdr_by'`: Benjamini-Yekutieli FDR


#### `analise_poder_estatistico(effect_size, n_per_group, alpha=0.05, alternative='two-sided')`
Calculates statistical power (1-β).

**Example:**

```python
from qualis_a1_modules.statistical_extensions import analise_poder_estatistico

poder_info = analise_poder_estatistico(
    effect_size=0.5,  # Medium effect (Cohen's d)
    n_per_group=30,
    alpha=0.05
)
print(f"Power: {poder_info['poder']:.2%}")
print(f"Interpretation: {poder_info['interpretacao']}")

```text

#### Power Interpretation:
- < 0.50: Inadequate
- 0.50-0.70: Low (consider increasing n)
- 0.70-0.80: Acceptable
- ≥ 0.80: Adequate (scientific standard)
- ≥ 0.90: Excellent


#### `calcular_tamanho_amostral_necessario(effect_size, poder_desejado=0.80, alpha=0.05)`
Calculates required sample size for desired power.

**Example:**

```python
from qualis_a1_modules.statistical_extensions import calcular_tamanho_amostral_necessario

n_necessario = calcular_tamanho_amostral_necessario(
    effect_size=0.5,
    poder_desejado=0.80
)
print(f"Need {n_necessario} per group ({n_necessario*2} total)")

```text

### References
- Dunn (1961), "Multiple comparisons among means"
- Cohen (1988), "Statistical Power Analysis for the Behavioral Sciences"
- Benjamini & Hochberg (1995), "Controlling the false discovery rate"


---


## 📝 Module 4: Auditing (`auditing.py`)

### Purpose
Generate explicit mappings between scientific methods and code implementation.

### Key Function

#### `gerar_tabela_codigo_metodo(pasta_resultados, formato='csv')`
Creates traceability table for auditing.

**Example:**

```python
from qualis_a1_modules.auditing import gerar_tabela_codigo_metodo

# Generate in multiple formats
for fmt in ['csv', 'json', 'markdown']:
    path = gerar_tabela_codigo_metodo('./results', formato=fmt)
    print(f"Generated: {path}")

```text

**Table Format:**

| Component | File | Class/Function | Lines | Parameters | Equation | Reference |
|-----------|------|----------------|-------|------------|----------|-----------|
| Depolarizing Channel | framework_investigativo_completo.py | RuidoDepolarizante | 1909-1934 | nivel: float | ℰ(ρ) = (1-p)ρ + p/3(XρX + YρY + ZρZ) | Nielsen & Chuang (2010) |

#### Mappings Include:
- All 11 noise models
- QNG optimizer
- Validation functions
- Statistical methods
- Line numbers and parameters


### Compliance
Enables peer reviewers to verify:

- Correct implementation of methods
- Correspondence with paper descriptions
- Parameter choices and equations


### References
- Stodden (2014), "The scientific method in practice"


---


## 📊 Module 5: Visualization (`visualization.py`)

### Purpose
Generate publication-quality quantum circuit diagrams.

### Key Functions

#### `gerar_diagrama_circuito(ansatz_func, n_qubits, n_camadas, pasta_resultados, nome_circuito='quantum_circuit', formato='png', framework='pennylane')`
Generates circuit diagram for single ansatz.

**Example:**

```python
from qualis_a1_modules.visualization import gerar_diagrama_circuito
from framework_investigativo_completo import circuito_hardware_efficient

path = gerar_diagrama_circuito(
    ansatz_func=circuito_hardware_efficient,
    n_qubits=4,
    n_camadas=2,
    pasta_resultados='./diagrams',
    nome_circuito='hardware_efficient',
    formato='png'
)

```text

#### `gerar_todos_diagramas_ansatze(pasta_resultados, n_qubits=4, n_camadas=2, formato='png')`
Generates diagrams for all 9 ansätze.

**Example:**

```python
from qualis_a1_modules.visualization import gerar_todos_diagramas_ansatze

diagramas = gerar_todos_diagramas_ansatze(
    pasta_resultados='./diagrams',
    n_qubits=4,
    n_camadas=2
)
print(f"Generated {len(diagramas)} diagrams")

```text

**Supported Ansätze:**
1. Hardware-Efficient
2. Tree Tensor Network
3. QAOA-inspired
4. Alternating Layers
5. Star Entanglement
6. Brickwork
7. Random Entangling
8. Basic (SimplifiedTwoDesign)
9. Strongly Entangling


#### Supported Frameworks:
- PennyLane (default)
- Cirq
- Qiskit


#### Supported Formats:
- PNG (default, 300 DPI)
- PDF (publication quality)
- SVG (scalable)


### References
- Bergholm et al. (2018), "PennyLane"


---


## 🎯 Usage Example: Complete Workflow

```python
import numpy as np
import pandas as pd
from qualis_a1_modules import (
    configurar_seeds_reprodutiveis,
    gerar_manifesto_execucao,
    validar_operadores_kraus,
    testes_post_hoc_com_correcao,
    analise_poder_estatistico,
    gerar_tabela_codigo_metodo,
    gerar_todos_diagramas_ansatze
)

# 1. Configure reproducibility
seed_info = configurar_seeds_reprodutiveis(seed=42)

# 2. Validate noise operators
from qualis_a1_modules.validation import obter_operadores_kraus_depolarizante
ops = obter_operadores_kraus_depolarizante(0.01)
validar_operadores_kraus(ops)

# 3. Run experiments
# ... your experiment code here ...

# 4. Statistical analysis with corrections
df_resultados = pd.DataFrame({...})  # Your results
resultados_posthoc = testes_post_hoc_com_correcao(
    df_resultados, 'noise_type', 'accuracy'
)

# 5. Power analysis
poder = analise_poder_estatistico(effect_size=0.5, n_per_group=30)
print(f"Statistical power: {poder['poder']:.2%}")

# 6. Generate audit trail
gerar_tabela_codigo_metodo('./results', formato='csv')

# 7. Generate circuit diagrams
gerar_todos_diagramas_ansatze('./results/diagrams')

# 8. Generate execution manifest
config = {'version': '8.0-QAI', 'default_seed': 42}
gerar_manifesto_execucao(config, './results', seed_info=seed_info)

print("✓ All Qualis A1 compliance artifacts generated!")

```text

---


## 📋 Qualis A1 Compliance Checklist

### Rigor Matemático (30 pts)
- [x] Docstrings with LaTeX equations (10 pts)
- [x] Kraus operator validation (10 pts)
- [x] QNG mathematical derivation (10 pts)


### Reprodutibilidade (30 pts)
- [x] Centralized seed configuration (15 pts)
- [x] Execution manifests (15 pts)


### Rigor Estatístico (20 pts)
- [x] Bonferroni correction (10 pts)
- [x] Statistical power analysis (10 pts)


### Auditoria e Transparência (20 pts)
- [x] Code→Method traceability table (10 pts)
- [x] Circuit diagram generation (5 pts)
- [ ] Multi-framework support (5 pts) - Partial (PennyLane, Qiskit ready; Cirq WIP)


**Current Score: 95/100**


---


## 🔧 Installation

No additional dependencies beyond the main framework requirements:

```bash
pip install -r requirements.txt

```

Required packages:

- numpy >= 1.23.0
- pandas >= 2.0.0
- scipy >= 1.10.0
- statsmodels >= 0.14.0
- matplotlib >= 3.5.0
- pennylane >= 0.30.0


---


## 📚 References

### Key Publications
1. Nielsen & Chuang (2010). "Quantum Computation and Quantum Information." Cambridge University Press.
2. Preskill (2018). "Quantum Computing in the NISQ era and beyond." Quantum, 2, 79.
3. Stokes et al. (2020). "Quantum Natural Gradient." Quantum, 4, 269.
4. Cohen (1988). "Statistical Power Analysis for the Behavioral Sciences." 2nd ed.
5. Peng (2011). "Reproducible research in computational science." Science, 334(6060), 1226-1227.


### Statistical Methods
- Dunn (1961). "Multiple comparisons among means." JASA, 56(293), 52-64.
- Benjamini & Hochberg (1995). "Controlling the false discovery rate." JRSS B, 57(1), 289-300.


### Quantum Noise
- Clerk et al. (2010). "Introduction to quantum noise." Rev. Mod. Phys., 82(2), 1155.
- Kandala et al. (2019). "Error mitigation extends..." Nature, 567(7749), 491-495.


---


## 📝 License

Same as parent project (see main LICENSE file).

---


## 👥 Authors

Enhancement Module Development: Qualis A1 Improvement Initiative  
Original Framework: MarceloClaro

---


## 🙏 Acknowledgments

This module enhancement was designed to meet the rigorous standards of:

- Nature family journals
- Quantum (Verein zur Förderung des Open Access Publizierens)
- Physical Review series (American Physical Society)


All improvements follow best practices from the computational science and quantum computing communities.
