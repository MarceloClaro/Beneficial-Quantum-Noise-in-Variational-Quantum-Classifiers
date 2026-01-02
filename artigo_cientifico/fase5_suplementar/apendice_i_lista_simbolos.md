# APÊNDICE I: Lista Completa de Símbolos e Notação

**Data:** 02 de janeiro de 2026  
**Seção:** Apêndice I - Símbolos e Notação (~500 palavras)  
**Status:** Novo conteúdo para expansão Qualis A1

---

## I.1 SÍMBOLOS MATEMÁTICOS PRINCIPAIS

### I.1.1 Espaços e Conjuntos

| Símbolo | Descrição | Primeira Aparição |
|---------|-----------|-------------------|
| $\mathcal{H}$ | Espaço de Hilbert de $n$ qubits, $\mathcal{H} = \mathbb{C}^{2^n}$ | Seção 3.1 |
| $\mathcal{B}(\mathcal{H})$ | Espaço de operadores lineares em $\mathcal{H}$ | Seção 3.1 |
| $\mathcal{D}(\mathcal{H})$ | Espaço de operadores densidade (matrizes densidade) | Seção 3.1 |
| $\mathcal{X}$ | Espaço de entrada (features), $\mathcal{X} \subseteq \mathbb{R}^d$ | Seção 3.1 |
| $\mathcal{Y}$ | Espaço de saída (labels), $\mathcal{Y} = \{0, 1\}$ | Seção 3.1 |
| $\Theta$ | Espaço de parâmetros, $\Theta \subseteq \mathbb{R}^p$ | Seção 3.1 |
| $\mathcal{D}_{train}$ | Conjunto de treinamento, $\{(x_i, y_i)\}_{i=1}^N$ | Seção 3.1 |
| $\mathcal{D}_{test}$ | Conjunto de teste | Seção 6 |

### I.1.2 Estados Quânticos e Operadores

| Símbolo | Descrição |
|---------|-----------|
| $\|\psi\rangle$ | Vetor de estado puro em $\mathcal{H}$ |
| $\rho$ | Operador densidade (estado misto), $\rho \in \mathcal{D}(\mathcal{H})$ |
| $\rho_{diag}$ | Parte diagonal de $\rho$ (populações) |
| $\rho_{off}$ | Parte off-diagonal de $\rho$ (coerências) |
| $U(\theta)$ | Operador unitário parametrizado |
| $\hat{O}$ | Observável Hermitiano, $\hat{O} = \hat{O}^\dagger$ |
| $\sigma_x, \sigma_y, \sigma_z$ | Matrizes de Pauli |
| $I$ ou $\mathbb{I}$ | Operador identidade |

### I.1.3 Canais Quânticos

| Símbolo | Descrição |
|---------|-----------|
| $\Phi_\gamma$ | Canal quântico CPTP com intensidade de ruído $\gamma$ |
| $\Phi_{pd}$ | Canal de Phase Damping |
| $\Phi_{dep}$ | Canal de Depolarizing |
| $\Phi_{ad}$ | Canal de Amplitude Damping |
| $\Phi_{bf}$ | Canal de Bit Flip |
| $\Phi_{pf}$ | Canal de Phase Flip |
| $K_k$ | Operadores de Kraus, $\Phi(\rho) = \sum_k K_k \rho K_k^\dagger$ |
| $L_j$ | Operadores de Lindblad (jump operators) |

### I.1.4 Parâmetros e Hiperparâmetros

| Símbolo | Descrição | Valor Típico |
|---------|-----------|--------------|
| $n$ | Número de qubits | 4-6 |
| $p$ | Número de parâmetros treináveis | 20-80 |
| $N$ | Tamanho do conjunto de treinamento | 280 |
| $\gamma$ | Intensidade de ruído quântico | $[10^{-5}, 10^{-1}]$ |
| $\gamma^*$ | Intensidade ótima de ruído | $\sim 0.001$ |
| $\eta$ | Taxa de aprendizado | 0.01-0.1 |
| $T$ | Número de épocas de treinamento | 100-500 |
| $\delta$ | Nível de confiança estatística | 0.05 |

---

## I.2 FUNÇÕES E MÉTRICAS

### I.2.1 Funções de Perda

| Símbolo | Descrição |
|---------|-----------|
| $\mathcal{L}_{train}(\theta)$ | Perda empírica no conjunto de treinamento |
| $\mathcal{L}_{gen}(\theta)$ | Perda de generalização (erro verdadeiro) |
| $\Delta_{gen}$ | Gap de generalização, $\mathcal{L}_{gen} - \mathcal{L}_{train}$ |
| $\ell(y, \hat{y})$ | Função de perda pontual (cross-entropy, hinge) |

### I.2.2 Métricas de Complexidade

| Símbolo | Descrição |
|---------|-----------|
| $\hat{\mathcal{R}}_N(\mathcal{F})$ | Complexidade de Rademacher empírica |
| $\mathcal{C}(\rho)$ | Magnitude de coerências, $\\|\rho_{off}\\|_F$ |
| $\text{rank}_{eff}(\mathcal{F})$ | Rank efetivo da QFIM |

### I.2.3 Geometria Quântica

| Símbolo | Descrição |
|---------|-----------|
| $\mathcal{F}$ | Matriz de Informação de Fisher Quântica (QFIM) |
| $g_{ij}^{FS}$ | Métrica de Fubini-Study |
| $d_{FS}$ | Distância geodésica de Fubini-Study |
| $F(\rho, \sigma)$ | Fidelidade quântica |

---

## I.3 OPERADORES E NOTAÇÃO

### I.3.1 Operações Lineares

| Notação | Significado |
|---------|-------------|
| $\langle \cdot \rangle$ | Valor esperado |
| $\text{Tr}[\cdot]$ | Traço de operador |
| $\\|\cdot\\|$ | Norma de operador |
| $\\|\cdot\\|_F$ | Norma de Frobenius |
| $\\|\cdot\\|_1$ | Norma traço |
| $\odot$ | Produto de Hadamard (elemento-wise) |
| $\otimes$ | Produto tensorial |
| $\circ$ | Composição de funções/canais |

### I.3.2 Derivadas e Gradientes

| Notação | Significado |
|---------|-------------|
| $\partial_i$ ou $\frac{\partial}{\partial \theta_i}$ | Derivada parcial com respeito a $\theta_i$ |
| $\nabla_\theta$ | Gradiente com respeito a $\theta$ |
| $\partial_\theta \mathcal{L}$ | Vetor gradiente da perda |

### I.3.3 Expectativas e Probabilidades

| Notação | Significado |
|---------|-------------|
| $\mathbb{E}[\cdot]$ | Expectativa |
| $\mathbb{E}_{\mathcal{D}}[\cdot]$ | Expectativa sobre datasets |
| $\text{Var}[\cdot]$ | Variância |
| $P(y\|x, \theta)$ | Probabilidade condicional |

---

## I.4 HIPÓTESES E CONDIÇÕES

### I.4.1 Hipóteses Principais do Teorema 1

| Label | Descrição Curta | Condição Matemática |
|-------|-----------------|---------------------|
| **H1** | Superparametrização | $\text{rank}_{eff}(\mathcal{F}) > N$ |
| **H2** | Amostra Finita | $N < C \cdot \sqrt{p}$ |
| **H3** | Coerências Espúrias | $\\|\rho_{off}\\|_F > \epsilon = O(1/\sqrt{N})$ |

### I.4.2 Condições Auxiliares

| Label | Descrição |
|-------|-----------|
| **C1** | CPTP do Canal | $\sum_k K_k^\dagger K_k = I$ |
| **C2** | Normalização | $\text{Tr}[\rho] = 1$ |
| **C3** | Positividade | $\rho \geq 0$ (semidefinido positivo) |

---

## I.5 SCHEDULES DE RUÍDO

| Nome | Fórmula |
|------|---------|
| **Static** | $\gamma(t) = \gamma_0$ (constante) |
| **Linear** | $\gamma(t) = \gamma_{max}(1 - t/T)$ |
| **Exponential** | $\gamma(t) = \gamma_{max} e^{-\lambda t/T}$ |
| **Cosine** | $\gamma(t) = \gamma_{max} \cos^2(\pi t / 2T)$ |

---

## I.6 ABREVIAÇÕES

| Abreviação | Termo Completo |
|------------|----------------|
| **VQC** | Variational Quantum Classifier |
| **VQA** | Variational Quantum Algorithm |
| **PQC** | Parametrized Quantum Circuit |
| **NISQ** | Noisy Intermediate-Scale Quantum |
| **QFIM** | Quantum Fisher Information Matrix |
| **FS** | Fubini-Study |
| **AUEC** | Adaptive Unified Error Correction |
| **ZNE** | Zero-Noise Extrapolation |
| **QEC** | Quantum Error Correction |
| **CPTP** | Completely Positive Trace-Preserving |
| **POVM** | Positive Operator-Valued Measure |

---

## I.7 CONVENÇÕES DE NOTAÇÃO

1. **Vetores:** Minúsculas com seta ou ket: $\vec{v}$, $|v\rangle$
2. **Matrizes:** Maiúsculas sem adorno: $M$, $\rho$
3. **Operadores:** Maiúsculas com chapéu: $\hat{O}$, $\hat{H}$
4. **Espaços:** Caligráficos: $\mathcal{H}$, $\mathcal{X}$
5. **Parâmetros:** Letras gregas: $\theta$, $\gamma$, $\eta$
6. **Índices:** Subscritos: $\theta_i$, $\rho_{ij}$

---

**Contagem de Palavras:** ~550 ✅

**Status:** Apêndice I completo ✅
