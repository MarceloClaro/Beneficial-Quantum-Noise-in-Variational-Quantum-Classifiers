# QUICKSTART: Advanced Quantum Framework V8.0

## 🚀 Quick Start (5 minutos)

### Instalação Rápida

```bash
# 1. Instalar dependências básicas
pip install numpy pandas scipy scikit-learn matplotlib plotly

# 2. Testar framework
python example_framework_v8_quick.py
```

### Instalação Completa com DeepChem

```bash
# Opção 1: Conda (Recomendado)
bash install_deepchem.sh 3.10 cpu

# Opção 2: Pip
pip install deepchem
conda install -c conda-forge rdkit
```

## 📊 Uso Básico

### Exemplo 1: Dataset Sintético

```python
from framework_quantum_advanced_v8 import AdvancedConfig, AdvancedVQC, DeepChemDatasetLoader

# Configuração
config = AdvancedConfig(
    framework="pennylane",
    n_qubits=4,
    n_layers=2,
    n_epochs=30,
    error_mitigation="zne"
)

# Carregar dados
loader = DeepChemDatasetLoader()
dataset = loader.load_dataset("BACE", max_samples=200)

# Treinar
vqc = AdvancedVQC(config)
vqc.fit(dataset['X_train'], dataset['y_train'])

# Avaliar
accuracy = vqc.score(dataset['X_test'], dataset['y_test'])
print(f"Test Accuracy: {accuracy:.4f}")
```

### Exemplo 2: Comparação Multi-Framework

```python
from framework_quantum_advanced_v8 import AdvancedConfig, AdvancedVQC

frameworks = ["pennylane", "qiskit", "cirq"]
for fw in frameworks:
    config = AdvancedConfig(framework=fw, n_epochs=20, verbose=False)
    vqc = AdvancedVQC(config)
    vqc.fit(X_train, y_train)
    acc = vqc.score(X_test, y_test)
    print(f"{fw}: {acc:.4f}")
```

### Exemplo 3: Mitigação de Erros

```python
from framework_quantum_advanced_v8 import ZeroNoiseExtrapolation

# ZNE
zne = ZeroNoiseExtrapolation(scale_factors=[1.0, 1.5, 2.0])
mitigated = zne.extrapolate([0.85, 0.80, 0.75])
print(f"Mitigated value: {mitigated:.4f}")
```

## 🔧 Configurações Comuns

### Reduzir Ruído

```python
config = AdvancedConfig(
    noise_level=0.001,  # Baixo ruído
    error_mitigation="combined"  # ZNE + TREX + AUEC
)
```

### Aumentar Capacidade

```python
config = AdvancedConfig(
    n_qubits=6,  # Mais qubits
    n_layers=3,  # Mais camadas
    n_epochs=100  # Mais treinamento
)
```

### Acelerar Treinamento

```python
config = AdvancedConfig(
    n_epochs=10,  # Menos épocas
    batch_size=64,  # Batches maiores
    verbose=False  # Sem logs detalhados
)
```

## 📈 Resultados Esperados

### Performance Típica (100 amostras, 30 épocas)

| Dataset | Train Acc | Test Acc | Tempo |
|---------|-----------|----------|-------|
| BACE    | 0.85-0.90 | 0.75-0.85| ~2min |
| HIV     | 0.80-0.85 | 0.70-0.80| ~2min |
| Tox21   | 0.78-0.82 | 0.68-0.75| ~2min |

### Complexidade do Circuito

| Qubits | Layers | Depth | Gates | Plateau Risk |
|--------|--------|-------|-------|--------------|
| 2      | 1      | 1     | 9     | LOW          |
| 4      | 2      | 2     | 30    | MEDIUM       |
| 6      | 3      | 3     | 63    | MEDIUM       |
| 8      | 4      | 4     | 108   | HIGH         |

## 🧪 Testes

```bash
# Testar tudo
python -m pytest tests/test_framework_advanced_v8.py -v

# Testar componente específico
python -m pytest tests/test_framework_advanced_v8.py::TestZeroNoiseExtrapolation -v

# Teste rápido
python example_framework_v8_quick.py
```

## 📚 Documentação Completa

- [README Completo](README_FRAMEWORK_ADVANCED_V8.md) - Documentação detalhada
- [Testes](tests/test_framework_advanced_v8.py) - Suite de testes
- [Exemplo](example_framework_v8_quick.py) - Exemplo funcional

## 🐛 Troubleshooting

### DeepChem não instala

```bash
# Use conda
conda install -c conda-forge deepchem rdkit
```

### ImportError: No module named 'pennylane'

```bash
pip install pennylane>=0.33.0
```

### Testes falhando

```bash
# Instalar dependências de teste
pip install pytest numpy pandas scipy scikit-learn
```

## 💡 Dicas

1. **Comece pequeno**: Use poucos qubits/layers primeiro
2. **Use synthetic data**: Teste sem DeepChem primeiro
3. **Monitore o tempo**: Reduz épocas para testes rápidos
4. **Veja os logs**: Use `verbose=True` para debug

## 🔗 Links Úteis

- [PennyLane Docs](https://pennylane.ai)
- [Qiskit Docs](https://qiskit.org)
- [DeepChem Docs](https://deepchem.readthedocs.io)
- [Repository Issues](https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers/issues)

## 📞 Suporte

Problemas? Abra uma issue no GitHub com:
- Versão do Python
- Comando executado
- Mensagem de erro completa
- Sistema operacional

---

**Pronto para começar?** Execute: `python example_framework_v8_quick.py` 🚀
