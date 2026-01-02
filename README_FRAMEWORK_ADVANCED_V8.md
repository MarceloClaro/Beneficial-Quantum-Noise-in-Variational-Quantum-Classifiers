# Framework Quantum Advanced V8.0

## 📋 Visão Geral

Framework avançado de Quantum Machine Learning que integra múltiplos frameworks quânticos, técnicas de mitigação de erros de última geração e integração com DeepChem para benchmarking molecular.

### 🚀 Principais Funcionalidades

1. **Multi-Framework Support**
   - PennyLane (default)
   - Qiskit (IBM Quantum)
   - Cirq (Google Quantum)

2. **Algoritmos Híbridos VQE/QAOA**
   - Variational Quantum Eigensolver (VQE)
   - Quantum Approximate Optimization Algorithm (QAOA)
   - Implementações híbridas otimizadas

3. **Mitigação de Erros Quânticos**
   - **ZNE (Zero-Noise Extrapolation)**: Extrapolação para ruído zero
   - **TREX (Twirled Readout Error eXtinction)**: Correção de erros de leitura
   - **AUEC (Adaptive Unified Error Correction)**: Correção adaptativa unificada

4. **Integração DeepChem**
   - 3 datasets moleculares: BACE, HIV, Tox21
   - Featurização molecular automática
   - Redução de dimensionalidade via PCA

5. **Análise de Complexidade Quântica**
   - Profundidade e contagem de portas
   - Análise de conectividade
   - Avaliação de risco de "barren plateaus"
   - Estimativa de expressividade

6. **Validação de Fórmulas de Predição de Ruído**
   - Comparação entre predições teóricas e resultados experimentais
   - Benchmarking contra algoritmos clássicos

## 📦 Instalação

### Pré-requisitos

```bash
Python >= 3.10
pip >= 21.0
```

### Dependências Principais

```bash
pip install pennylane>=0.33.0
pip install qiskit>=1.0.0 qiskit-aer>=0.13.0
pip install cirq>=1.0.0
pip install numpy pandas scipy scikit-learn
pip install matplotlib plotly
```

### Instalação DeepChem

O DeepChem é essencial para os datasets moleculares. Use o script de instalação fornecido:

```bash
# Instalação via Conda (Recomendado)
bash install_deepchem.sh 3.10 cpu

# Ou instalação manual via pip
pip install deepchem
conda install -c conda-forge rdkit  # RDKit é necessário
```

Para mais detalhes, consulte:
- [DeepChem Installation Guide](https://deepchem.readthedocs.io/en/latest/get_started/installation.html)
- [DeepChem Examples](https://deepchem.readthedocs.io/en/latest/get_started/examples.html)

### Instalação do Framework

```bash
# Clone o repositório
git clone https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers.git
cd Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers

# Instale dependências
pip install -r requirements.txt

# Instale DeepChem
bash install_deepchem.sh

# Verifique a instalação
python -c "import framework_quantum_advanced_v8; print('✓ Framework carregado')"
```

## 🎯 Uso Rápido

### Execução Básica

```python
from framework_quantum_advanced_v8 import main

# Executar com configuração padrão
main()
```

### Configuração Personalizada

```python
from framework_quantum_advanced_v8 import AdvancedConfig, AdvancedVQC, DeepChemDatasetLoader

# Configuração
config = AdvancedConfig(
    framework="pennylane",
    n_qubits=4,
    n_layers=2,
    n_epochs=50,
    learning_rate=0.01,
    noise_level=0.01,
    error_mitigation="zne",  # ou "trex", "auec", "combined"
    results_dir="meus_resultados",
    verbose=True
)

# Carregar dataset
loader = DeepChemDatasetLoader(verbose=True)
dataset = loader.load_dataset("BACE", max_samples=500)

# Treinar VQC
vqc = AdvancedVQC(config)
vqc.fit(dataset['X_train'], dataset['y_train'])

# Avaliar
test_acc = vqc.score(dataset['X_test'], dataset['y_test'])
print(f"Test Accuracy: {test_acc:.4f}")
```

### Comparação Multi-Framework

```python
from framework_quantum_advanced_v8 import AdvancedConfig, AdvancedVQC

frameworks = ["pennylane", "qiskit", "cirq"]
results = []

for fw in frameworks:
    config = AdvancedConfig(framework=fw, n_epochs=30)
    vqc = AdvancedVQC(config)
    vqc.fit(X_train, y_train)
    acc = vqc.score(X_test, y_test)
    results.append({"framework": fw, "accuracy": acc})

print(pd.DataFrame(results))
```

## 📊 Estrutura de Saída

Após execução, o framework gera:

```
resultados_advanced_v8/
├── results_summary.csv          # Resumo dos resultados
├── SUMMARY.md                    # Relatório em markdown
├── complexity_BACE.md            # Análise de complexidade (BACE)
├── complexity_HIV.md             # Análise de complexidade (HIV)
└── complexity_Tox21.md           # Análise de complexidade (Tox21)
```

### Exemplo de results_summary.csv

| dataset | train_accuracy | test_accuracy | training_time | framework | error_mitigation |
|---------|----------------|---------------|---------------|-----------|------------------|
| BACE    | 0.8542         | 0.7891        | 145.32        | pennylane | zne              |
| HIV     | 0.8123         | 0.7654        | 152.18        | pennylane | zne              |
| Tox21   | 0.7892         | 0.7234        | 138.45        | pennylane | zne              |

## 🔬 Fundamentos Científicos

### Zero-Noise Extrapolation (ZNE)

ZNE é uma técnica de mitigação de erros que:

1. Executa o circuito quântico em múltiplos níveis de ruído artificialmente amplificados
2. Ajusta uma curva de extrapolação (linear, exponencial ou polinomial)
3. Extrapola o resultado para o limite de ruído zero

**Referências:**
- Temme et al. (2017). "Error mitigation for short-depth quantum circuits"
- LaRose & Mari (2021). "Mitiq: A software package for error mitigation"

### TREX (Twirled Readout Error eXtinction)

TREX corrige erros sistemáticos de medição através de:

1. Calibração da matriz de confusão de readout
2. Inversão matricial para recuperar distribuição ideal
3. Aplicação de suavização Bayesiana (opcional)

**Referências:**
- Nation et al. (2021). "Scalable mitigation of measurement errors"
- Bravyi et al. (2021). "Mitigating measurement errors in multiqubit experiments"

### AUEC (Adaptive Unified Error Correction)

AUEC é uma contribuição inovadora que unifica:

1. Correção de erros de porta (gate errors)
2. Correção de decoerência
3. Correção de erros não-estacionários (drift)

Usando:
- Filtro de Kalman Estendido Quântico
- Controle adaptativo em tempo real
- Meta-aprendizado Bayesiano

### Análise de Complexidade

O framework analisa:

- **Circuit Depth**: Número de camadas sequenciais
- **Gate Count**: Total de portas quânticas
- **Connectivity**: Grau de emaranhamento
- **Expressibility**: Capacidade de cobrir espaço de Hilbert
- **Barren Plateau Risk**: Risco de platôs de gradiente

## 📚 Datasets DeepChem

### BACE (β-secretase)

- **Tarefa**: Classificação binária de atividade inibitória
- **Amostras**: ~1,500 compostos
- **Aplicação**: Drug discovery para Alzheimer

### HIV

- **Tarefa**: Classificação binária de atividade anti-HIV
- **Amostras**: ~41,000 compostos
- **Aplicação**: Descoberta de antivirais

### Tox21

- **Tarefa**: Multi-task de toxicidade (12 ensaios)
- **Amostras**: ~8,000 compostos
- **Aplicação**: Avaliação de segurança de drogas

## 🔧 Configurações Avançadas

### Parâmetros de ZNE

```python
config = AdvancedConfig(
    error_mitigation="zne",
    zne_scale_factors=[1.0, 1.5, 2.0, 2.5, 3.0],  # Fatores de escala
)
```

Métodos de extrapolação disponíveis:
- `"linear"`: E(λ) = a + b·λ
- `"exponential"`: E(λ) = a + b·exp(-c·λ) (padrão)
- `"polynomial"`: E(λ) = a + b·λ + c·λ²

### Configuração de Ruído

```python
config = AdvancedConfig(
    noise_level=0.01,          # Nível de ruído (0 a 0.2)
    noise_type="depolarizing", # Tipo de ruído
)
```

Tipos de ruído suportados:
- `depolarizing`: Despolarização uniforme
- `amplitude_damping`: Perda de energia
- `phase_damping`: Perda de coerência de fase

### Otimização de Hiperparâmetros

```python
import optuna

def objective(trial):
    config = AdvancedConfig(
        n_qubits=trial.suggest_int("n_qubits", 2, 6),
        n_layers=trial.suggest_int("n_layers", 1, 4),
        learning_rate=trial.suggest_float("lr", 0.001, 0.1, log=True),
        noise_level=trial.suggest_float("noise", 0.0, 0.05)
    )
    
    vqc = AdvancedVQC(config)
    vqc.fit(X_train, y_train)
    return vqc.score(X_test, y_test)

study = optuna.create_study(direction="maximize")
study.optimize(objective, n_trials=50)
```

## 🧪 Testes

```bash
# Executar testes
pytest tests/test_framework_advanced_v8.py -v

# Teste rápido
python -c "from framework_quantum_advanced_v8 import main; main()"
```

## 📖 Documentação Adicional

- [TREX Error Mitigation](trex_error_mitigation.py) - Implementação completa do TREX
- [AUEC](adaptive_unified_error_correction.py) - Correção adaptativa unificada
- [Framework Investigativo](framework_investigativo_completo.py) - Framework base
- [Artigo Científico](artigo_cientifico/) - Resultados e análises

## 🤝 Integração com Artigo Científico

O framework foi desenvolvido para ser **obrigatoriamente aplicado conforme mostrado na pasta artigo_cientifico**:

```bash
# Ver estrutura do artigo
ls -la artigo_cientifico/

# Documentos principais:
# - RESUMO_EXECUTIVO_FRAMEWORK.md: Resumo executivo
# - README.md: Instruções de uso
# - artigo_completo_qualis_a1.tex: Artigo LaTeX
```

## 🔗 Referências

### Quantum Computing

1. McClean et al. (2018). "Barren plateaus in quantum neural network training landscapes." Nature Communications.
2. Cerezo et al. (2021). "Variational quantum algorithms." Nature Reviews Physics.
3. Schuld & Petruccione (2018). "Supervised Learning with Quantum Computers." Springer.

### Error Mitigation

4. Temme et al. (2017). "Error mitigation for short-depth quantum circuits." Physical Review X.
5. LaRose & Mari (2021). "Mitiq: A software package for error mitigation." Quantum.
6. Nation et al. (2021). "Scalable mitigation of measurement errors." PRX Quantum.

### DeepChem & Molecular ML

7. Ramsundar et al. (2019). "Deep Learning for the Life Sciences." O'Reilly.
8. Wu et al. (2018). "MoleculeNet: A benchmark for molecular machine learning." Chemical Science.

## 📄 Licença

Este projeto está licenciado sob a licença especificada no repositório.

## 👥 Contribuidores

- Framework Beneficial Quantum Noise Team
- Baseado em `framework_investigativo_completo.py`

## 📧 Contato

Para questões e suporte, abra uma issue no repositório GitHub.

---

**Versão:** 8.0  
**Data:** 2026-01-02  
**Status:** Production Ready ✅
