# Framework Qiskit - Guia de Uso Completo

## 🎯 Visão Geral

Este documento descreve o **Framework Qiskit v7.2** para análise de ruído quântico benéfico em classificadores variacionais quânticos (VQCs). Esta implementação utiliza IBM Qiskit em vez de PennyLane, mantendo toda a funcionalidade do framework original.

## 🚀 Instalação

### Pré-requisitos

```bash
# Python 3.9 ou superior
python --version
```

### Instalação Completa

```bash
# 1. Clone o repositório
git clone https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers.git
cd Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers

# 2. Instale todas as dependências (PennyLane + Qiskit)
pip install -r requirements.txt

# 3. OU instale apenas Qiskit (se já tiver outras dependências)
pip install qiskit qiskit-aer qiskit-ibm-runtime
```

### Verificar Instalação

```python
import qiskit
print(f"Qiskit version: {qiskit.__version__}")

from qiskit_aer import AerSimulator
print("✓ Qiskit-Aer disponível")

from framework_qiskit import QISKIT_AVAILABLE
if QISKIT_AVAILABLE:
    print("✓ Framework Qiskit pronto para uso!")
```

## 📚 Estrutura do Framework

### Arquivos Principais

```
framework_qiskit.py              # Framework principal Qiskit
examples/exemplo_qiskit_completo.py  # Exemplos de uso
docs/GUIA_QISKIT.md              # Este guia
```

### Componentes do Framework

1. **ClassificadorVQCQiskit**: Classificador quântico variacional
2. **Modelos de Ruído**: Depolarizante, Amplitude/Phase Damping, etc.
3. **Arquiteturas de Circuitos**: 7 arquiteturas diferentes
4. **Visualizações**: Bloch sphere, State City 3D, Q-Sphere
5. **Datasets**: Iris, Moons, Circles (compatível com PennyLane)

## 🔧 Uso Básico

### Exemplo 1: Classificador Simples

```python
from framework_qiskit import ClassificadorVQCQiskit, carregar_datasets

# Carregar dados
datasets = carregar_datasets()
dataset = datasets['moons']

# Criar classificador
vqc = ClassificadorVQCQiskit(
    n_qubits=4,
    n_camadas=2,
    arquitetura='basico',
    tipo_ruido='sem_ruido',
    n_epocas=20,
    seed=42
)

# Treinar
vqc.fit(dataset['X_train'], dataset['y_train'])

# Avaliar
acuracia = vqc.score(dataset['X_test'], dataset['y_test'])
print(f"Acurácia: {acuracia:.4f}")
```

### Exemplo 2: Com Ruído Quântico

```python
# Classificador com ruído depolarizante
vqc_ruido = ClassificadorVQCQiskit(
    n_qubits=4,
    n_camadas=2,
    arquitetura='strongly_entangling',
    tipo_ruido='depolarizante',
    nivel_ruido=0.01,  # 1% de erro
    n_epocas=20,
    seed=42
)

vqc_ruido.fit(dataset['X_train'], dataset['y_train'])
acuracia = vqc_ruido.score(dataset['X_test'], dataset['y_test'])
```

### Exemplo 3: Visualizações Quânticas

```python
from framework_qiskit import (
    visualizar_bloch_sphere,
    visualizar_state_city_3d,
    visualizar_qsphere
)

# Treinar modelo
vqc.fit(dataset['X_train'], dataset['y_train'])

# Exemplo de entrada
x = dataset['X_test'][0]

# Gerar visualizações
visualizar_bloch_sphere(vqc, x, 'bloch_sphere.png')
visualizar_state_city_3d(vqc, x, 'state_city_3d.png')
visualizar_qsphere(vqc, x, 'qsphere.png')

# Diagrama do circuito
vqc.get_circuit_diagram('circuito.png')
```

## 🏗️ Arquiteturas de Circuitos

O framework suporta 7 arquiteturas diferentes:

| Arquitetura | Descrição | N° Parâmetros | Uso Recomendado |
|-------------|-----------|---------------|-----------------|
| `basico` | Rotações RY+RZ + CNOT cadeia | 3×n_qubits×n_camadas | Baseline, datasets simples |
| `strongly_entangling` | All-to-all entanglement | 3×n_qubits×n_camadas | Datasets complexos |
| `hardware_efficient` | Apenas RY + CNOT linear | n_qubits×n_camadas | Dispositivos reais |
| `alternating_layers` | Camadas alternadas | 2×n_qubits×n_camadas | Problemas com simetria |
| `brickwork` | Padrão brickwork | 2×n_qubits×n_camadas | Datasets médios |
| `random_entangling` | Entanglement aleatório | 2×n_qubits×n_camadas | Exploração |

### Escolhendo uma Arquitetura

```python
# Para datasets simples (Iris, Moons)
vqc = ClassificadorVQCQiskit(arquitetura='basico')

# Para datasets complexos
vqc = ClassificadorVQCQiskit(arquitetura='strongly_entangling')

# Para hardware real (menos gates)
vqc = ClassificadorVQCQiskit(arquitetura='hardware_efficient')
```

## 🔬 Modelos de Ruído Quântico

### Tipos de Ruído Disponíveis

1. **Sem Ruído** (`sem_ruido`)
   - Simulação ideal
   - Baseline para comparação

2. **Depolarizante** (`depolarizante`)
   - Ruído isotrópico
   - Mais comum em dispositivos reais
   ```python
   vqc = ClassificadorVQCQiskit(
       tipo_ruido='depolarizante',
       nivel_ruido=0.01
   )
   ```

3. **Amplitude Damping** (`amplitude_damping`)
   - Perda de energia (T1 relaxation)
   - Simula decaimento para |0⟩
   ```python
   vqc = ClassificadorVQCQiskit(
       tipo_ruido='amplitude_damping',
       nivel_ruido=0.005
   )
   ```

4. **Phase Damping** (`phase_damping`)
   - Perda de coerência (T2 dephasing)
   - Preserva população, perde fase
   ```python
   vqc = ClassificadorVQCQiskit(
       tipo_ruido='phase_damping',
       nivel_ruido=0.007
   )
   ```

5. **Combinado** (`combinado`)
   - Mistura de depolarizante + amplitude damping
   - Mais realista para dispositivos reais

### Níveis de Ruído Recomendados

| Tipo de Ruído | Nível Baixo | Nível Moderado | Nível Alto |
|---------------|-------------|----------------|------------|
| Depolarizante | 0.001-0.005 | 0.01-0.02 | 0.03-0.05 |
| Amplitude Damping | 0.001-0.003 | 0.005-0.01 | 0.015-0.03 |
| Phase Damping | 0.001-0.007 | 0.01-0.02 | 0.03-0.05 |

**Região de Ruído Benéfico**: γ ∈ [0.001, 0.007] (descoberto empiricamente)

## 📊 Visualizações Avançadas

### 1. Esfera de Bloch

Visualiza o estado de qubits individuais:

```python
from framework_qiskit import visualizar_bloch_sphere

# Após treinar o modelo
x = dataset['X_test'][0]
visualizar_bloch_sphere(vqc, x, 'bloch_sphere.png')
```

**Interpretação**:
- Eixo Z: |0⟩ (polo norte) e |1⟩ (polo sul)
- Plano XY: Superposições
- Vetor Bloch: Estado puro do qubit

### 2. State City 3D

Visualização 3D da densidade de probabilidade:

```python
from framework_qiskit import visualizar_state_city_3d

visualizar_state_city_3d(vqc, x, 'state_city_3d.png')
```

**Interpretação**:
- Eixos X,Y: Estados computacionais
- Eixo Z (altura): Amplitude de probabilidade
- "Arranha-céus": Estados com alta probabilidade

### 3. Q-Sphere

Representação esférica do estado quântico completo:

```python
from framework_qiskit import visualizar_qsphere

visualizar_qsphere(vqc, x, 'qsphere.png')
```

**Interpretação**:
- Círculos: Estados da base computacional
- Tamanho: Magnitude da amplitude
- Cor: Fase quântica

### 4. Diagrama de Circuito

```python
# Formato matplotlib (PNG)
vqc.get_circuit_diagram('circuito.png')

# Formato texto (console)
print(vqc.get_circuit_diagram())
```

## 🧪 Experimentos Completos

### Script de Execução Rápida

```bash
# Executar exemplo interativo
python examples/exemplo_qiskit_completo.py
```

### Grid Search de Ruído Benéfico

```python
import pandas as pd
from framework_qiskit import ClassificadorVQCQiskit, carregar_datasets

niveis_ruido = [0.0, 0.001, 0.005, 0.01, 0.02, 0.05]
tipos_ruido = ['sem_ruido', 'depolarizante', 'phase_damping']

resultados = []
datasets = carregar_datasets()
dataset = datasets['moons']

for tipo in tipos_ruido:
    for nivel in niveis_ruido:
        if tipo == 'sem_ruido' and nivel > 0:
            continue
        
        vqc = ClassificadorVQCQiskit(
            n_qubits=4,
            n_camadas=2,
            arquitetura='basico',
            tipo_ruido=tipo,
            nivel_ruido=nivel,
            n_epocas=15,
            seed=42
        )
        
        vqc.fit(dataset['X_train'], dataset['y_train'])
        acc = vqc.score(dataset['X_test'], dataset['y_test'])
        
        resultados.append({
            'tipo_ruido': tipo,
            'nivel_ruido': nivel,
            'acuracia': acc
        })

df = pd.DataFrame(resultados)
df.to_csv('resultados_grid_search.csv', index=False)
```

## 🔄 Comparação com PennyLane

### Equivalências

| PennyLane | Qiskit | Observação |
|-----------|--------|------------|
| `qml.device('default.mixed')` | `AerSimulator(noise_model=...)` | Simulador com ruído |
| `@qml.qnode` | `QuantumCircuit` + execução manual | Circuitos parametrizados |
| `qml.RY`, `qml.RZ` | `qc.ry()`, `qc.rz()` | Rotações |
| `qml.CNOT` | `qc.cx()` | Porta controlada-NOT |
| `qml.expval(qml.PauliZ)` | Medição + cálculo de expectação | Valor esperado |
| `DepolarizingChannel` | `depolarizing_error()` | Ruído depolarizante |

### Vantagens do Qiskit

✅ **Acesso a hardware IBM Quantum** (dispositivos reais)  
✅ **Visualizações nativas ricas** (Bloch, State City, Q-Sphere)  
✅ **Comunidade IBM Quantum** e suporte empresarial  
✅ **Integração com IBM Cloud**  
✅ **Transpilação para hardware real**  

### Vantagens do PennyLane

✅ **Diferenciação automática nativa**  
✅ **Interface mais simples para ML**  
✅ **Integração com PyTorch/TensorFlow**  
✅ **Menos boilerplate code**  

## 🎯 Casos de Uso

### 1. Pesquisa Acadêmica

```python
# Experimento para artigo científico
from framework_qiskit import executar_experimento_qiskit

resultado = executar_experimento_qiskit(
    dataset_nome='iris',
    arquitetura='strongly_entangling',
    tipo_ruido='phase_damping',
    nivel_ruido=0.005,
    n_epocas=50,
    pasta_resultados='resultados_artigo',
    verbose=True
)
```

### 2. Prototipagem de Algoritmos

```python
# Testar nova arquitetura de circuito
def minha_arquitetura_custom(n_qubits, n_camadas):
    # Implementar arquitetura personalizada
    pass

# Registrar no framework
from framework_qiskit import ARQUITETURAS_QISKIT
ARQUITETURAS_QISKIT['minha_arquitetura'] = minha_arquitetura_custom
```

### 3. Educação e Demonstrações

```python
# Demo para aula/apresentação
vqc = ClassificadorVQCQiskit(
    n_qubits=2,  # Menos qubits = mais rápido
    n_camadas=1,
    arquitetura='basico',
    n_epocas=5
)

vqc.fit(X_train, y_train)
print(vqc.get_circuit_diagram())  # Mostrar circuito
```

## 🐛 Troubleshooting

### Erro: "Qiskit não está instalado"

```bash
pip install qiskit qiskit-aer qiskit-ibm-runtime
```

### Erro: "Noise model not compatible"

Certifique-se de usar `AerSimulator` com noise model:

```python
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel

# Correto
noise_model = NoiseModel()
simulator = AerSimulator(noise_model=noise_model)
```

### Aviso: "Shots muito baixo"

Aumente o número de medições para resultados mais estáveis:

```python
vqc = ClassificadorVQCQiskit(
    shots=2048  # Padrão: 1024, máximo: 8192
)
```

### Treinamento muito lento

**Soluções**:
1. Reduzir `n_epocas`
2. Usar `arquitetura='hardware_efficient'` (menos parâmetros)
3. Reduzir `shots` para prototipar
4. Usar GPU (se disponível via Qiskit-Aer-GPU)

```python
# Modo rápido (prototipagem)
vqc = ClassificadorVQCQiskit(
    n_qubits=3,
    n_camadas=1,
    n_epocas=5,
    shots=512,
    arquitetura='hardware_efficient'
)
```

## 📈 Performance e Otimização

### Dicas de Performance

1. **Use menos qubits para prototipagem**
   ```python
   vqc = ClassificadorVQCQiskit(n_qubits=3)  # Em vez de 4-5
   ```

2. **Cache de resultados**
   ```python
   # Salvar modelo treinado
   import pickle
   with open('modelo_vqc.pkl', 'wb') as f:
       pickle.dump(vqc, f)
   ```

3. **Paralelização (múltiplos experimentos)**
   ```python
   from joblib import Parallel, delayed
   
   def treinar_modelo(config):
       vqc = ClassificadorVQCQiskit(**config)
       vqc.fit(X_train, y_train)
       return vqc.score(X_test, y_test)
   
   configs = [...]  # Lista de configurações
   resultados = Parallel(n_jobs=4)(
       delayed(treinar_modelo)(cfg) for cfg in configs
   )
   ```

### Benchmarks Típicos

| Configuração | Tempo/Época | Total (20 épocas) |
|--------------|-------------|-------------------|
| 3 qubits, 1 camada | 2-3s | 40-60s |
| 4 qubits, 2 camadas | 5-8s | 100-160s |
| 5 qubits, 3 camadas | 15-25s | 300-500s |

*Testado em CPU Intel i7 2.6 GHz*

## 🔗 Recursos Adicionais

### Documentação Oficial

- [Qiskit Documentation](https://qiskit.org/documentation/)
- [Qiskit Textbook](https://qiskit.org/textbook/)
- [IBM Quantum Lab](https://quantum-computing.ibm.com/)

### Tutoriais

- [Qiskit Machine Learning](https://qiskit.org/documentation/machine-learning/)
- [Variational Algorithms](https://qiskit.org/textbook/ch-applications/vqc.html)

### Papers Relacionados

1. Schuld et al. (2020). "Circuit-centric quantum classifiers." Physical Review A.
2. Mitarai et al. (2018). "Quantum circuit learning." Physical Review A.
3. Havlíček et al. (2019). "Supervised learning with quantum-enhanced feature spaces." Nature.

## 🤝 Contribuindo

Encontrou um bug ou tem uma sugestão? Abra uma issue no GitHub!

## 📄 Licença

MIT License - Veja LICENSE para detalhes.

## ✨ Citação

Se usar este framework em sua pesquisa, por favor cite:

```bibtex
@software{framework_qiskit_vqc,
  title={Framework Qiskit para Análise de Ruído Quântico Benéfico em VQCs},
  author={Claro, Marcelo},
  year={2025},
  url={https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers}
}
```

---

**Framework Qiskit v7.2** | IBM Quantum Computing | 2025
