# Visualizações Qiskit - Resultados da Execução

## 📊 Demonstração Completa do Framework Qiskit

Este documento apresenta as visualizações geradas pelo framework Qiskit, demonstrando as capacidades exclusivas de visualização quântica em 3D.

**Data de Execução**: 24/12/2025  
**Tempo de Execução**: 1.2 segundos  
**Timeout Configurado**: 600 segundos (10 minutos)  
**Total de Visualizações**: 5 imagens


---


## 🔬 Visualizações Geradas

### 1. Circuito Quântico Básico

**Arquivo**: `01_circuito_basico.png`


Circuito VQC (Variational Quantum Classifier) básico com 2 qubits:

- **Encoding**: Rotações RY baseadas nos dados de entrada
- **Variational Layer**: Rotações parametrizadas e porta CNOT para entrelaçamento
- **Arquitetura**: Simples mas efetiva para classificação binária


![Circuito Básico](visualizacoes_qiskit_md/01_circuito_basico.png)

**Componentes**:
- `RY`: Rotações em torno do eixo Y (encoding de dados)
- `CX`: Porta CNOT (cria entrelaçamento quântico)


---


### 2. Esfera de Bloch - Sem Ruído

**Arquivo**: `02_bloch_sem_ruido.png`


Visualização 3D do estado quântico em uma **Esfera de Bloch**:

- Representa o estado de um qubit no espaço de Hilbert
- **Eixo Z**: Estados |0⟩ (polo norte) e |1⟩ (polo sul)
- **Equador**: Superposições com fase relativa


![Bloch Sphere](visualizacoes_qiskit_md/02_bloch_sem_ruido.png)

**Interpretação**:
- Posição do vetor indica o estado quântico puro
- Latitude: Probabilidade de medir |0⟩ vs |1⟩
- Longitude: Fase relativa


---


### 3. State City 3D - Sem Ruído 🏙️

**Arquivo**: `03_city3d_sem_ruido.png`


Visualização **"Plano Árido"** (State City) em 3D:

- Cada "arranha-céu" representa amplitude de um estado da base computacional
- Altura = Magnitude da amplitude
- Cor = Fase quântica


![State City 3D](visualizacoes_qiskit_md/03_city3d_sem_ruido.png)

**Características**:
- Visualização intuitiva de estados de múltiplos qubits
- Mostra distribuição de probabilidade completa
- Exclusivo do Qiskit (não disponível no PennyLane)


---


### 4. State City 3D - COM Ruído Phase Damping

**Arquivo**: `05_city3d_com_ruido.png`


Mesmo estado, mas após aplicação de **ruído quântico**:

- **Tipo de ruído**: Phase Damping (γ = 0.01)
- **Efeito físico**: Perda de coerência (dephasing)
- **Modelo**: Relaxação T₂


![State City com Ruído](visualizacoes_qiskit_md/05_city3d_com_ruido.png)

**Comparação com imagem anterior**:
- Amplitudes podem mudar ligeiramente
- Fase quântica é afetada (cores diferentes)
- Demonstra efeito do ruído realista em hardware quântico


---


### 5. Circuito com Entrelaçamento (4 qubits)

**Arquivo**: `06_circuito_entangled.png`


Circuito mais complexo com **4 qubits entrelaçados**:

- **Camada de Encoding**: 4 rotações RY
- **Camada de Entrelaçamento**: CNOTs em cadeia + cíclico
- **Camada Variacional**: Rotações RY e RZ


![Circuito Entrelaçado](visualizacoes_qiskit_md/06_circuito_entangled.png)

**Arquitetura**:
- Padrão de entrelaçamento em anel (circular)
- Maior expressividade do circuito
- Adequado para problemas de classificação mais complexos


---


## 📈 Características Exclusivas do Qiskit

### Comparação com PennyLane

| Recurso | Qiskit | PennyLane |
|---------|--------|-----------|
| **Bloch Sphere** | ✅ Nativo | ✅ Via bibliotecas externas |
| **State City 3D** | ✅ **Exclusivo** | ❌ Não disponível |
| **Q-Sphere** | ✅ **Exclusivo** | ❌ Não disponível |
| **Circuit Diagrams** | ✅ Múltiplos estilos | ✅ Limitado |
| **Hardware Real** | ✅ IBM Quantum | ⚠️ Limitado |
| **Visualizações 3D** | ✅ 4 tipos | ⚠️ 1-2 tipos |

---


## 🎯 Aplicações

### 1. Educação e Ensino
- Visualizações intuitivas de conceitos quânticos
- Demonstração de entrelaçamento e superposição
- Efeito do ruído quântico


### 2. Pesquisa
- Análise visual de estados quânticos complexos
- Comparação de circuitos com/sem ruído
- Validação de implementações


### 3. Documentação
- Figuras de alta qualidade para artigos
- Apresentações e relatórios
- Material educacional


---


## ⚙️ Configuração Técnica

### Parâmetros de Execução

```python
TIMEOUT = 600  # segundos (10 minutos)
DPI = 150      # Resolução das imagens
Backend = 'Agg'  # Matplotlib sem display

```text

### Modelo de Ruído

```python
noise_model = NoiseModel()
error_1q = phase_damping_error(0.01)  # γ = 0.01
error_2q = error_1q.tensor(error_1q)  # Para 2-qubit gates

```text

### Datasets

- **Dataset**: Make Moons (sklearn)
- **Amostras**: 50 (35 treino, 15 teste)
- **Features**: 2 (visualização 2D)
- **Classes**: 2 (classificação binária)


---


## 📊 Estatísticas

| Métrica | Valor |
|---------|-------|
| **Visualizações geradas** | 5 |
| **Tempo total** | 1.2 segundos |
| **Taxa de sucesso** | 83% (5/6 tentativas) |
| **Timeout usado** | 0.2% (1.2s / 600s) |
| **Tamanho médio por imagem** | ~150-300 KB |

---


## 🚀 Como Reproduzir

### Comando

```bash
python demo_qiskit_ultra_rapido.py

```text

### Requisitos

```bash
pip install numpy scikit-learn qiskit qiskit-aer matplotlib pylatexenc

```

### Timeout Configurável

O script aceita timeout de até **600 segundos** para garantir completude mesmo em ambientes mais lentos.

---


## 📝 Notas Técnicas

### Vantagens do Qiskit

1. **Visualizações Nativas**: Integradas diretamente no framework
2. **Hardware Real**: Preparado para execução em IBM Quantum
3. **Modelos de Ruído**: Realistas e configuráveis
4. **Qualidade de Publicação**: DPI configurável (150-300)


### Limitações

- Requer bibliotecas adicionais (pylatexenc para circuitos)
- Simulações com ruído são mais lentas (density matrix)
- Maior consumo de memória para estados grandes


---


## 🔗 Referências

- **Framework Qiskit**: `framework_qiskit.py`
- **Script de Demo**: `demo_qiskit_ultra_rapido.py`
- **Documentação Completa**: `docs/GUIA_QISKIT.md`
- **Qiskit Documentation**: <https://qiskit.org/documentation/>


---


**Status**: ✅ Execução completa em 1.2s  
**Próximo**: Integrar visualizações em documentação principal  
**Data**: 24/12/2025

