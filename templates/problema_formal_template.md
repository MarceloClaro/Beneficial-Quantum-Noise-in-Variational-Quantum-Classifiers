# Formal Problem Statement Template

## Definição do Problema

**Contexto**: [Descrever o contexto e motivação]

**Seja**:
- $\mathcal{D} = \{ (x_i, y_i) \}_{i=1}^N$ um dataset com $N$ amostras, onde $x_i \in \mathbb{R}^d$ são features e $y_i \in \{0, 1, \ldots, C-1\}$ são labels
- $U(\theta)$ um circuito quântico variacional parametrizado por $\theta \in \mathbb{R}^P$
- $\mathcal{N}_p(\cdot)$ um canal quântico de ruído com parâmetro $p \in [0, 1]$
- $L(\theta, p; \mathcal{D})$ a função de custo (loss) calculada sobre o dataset

**O problema de otimização é:**

$$
(\theta^*, p^*) = \arg\min_{\theta \in \Theta, p \in [0, p_{max}]} \mathbb{E}_{(x,y) \sim \mathcal{D}_{test}} \left[ L(y, f_{\theta,p}(x)) \right]
$$

onde:
- $f_{\theta,p}(x)$ é a predição do modelo quântico com parâmetros $\theta$ sob ruído $p$
- $\Theta$ é o espaço de parâmetros permitidos
- $p_{max}$ é o limite superior de ruído aceitável

**Sujeito às restrições**:
1. $p \in [0, p_{max}]$ onde $p_{max} < 1$ (ruído não pode ser máximo)
2. Arquitetura de $U(\theta)$ é fixa (número de qubits, depth, gate set)
3. Budget computacional: $T_{max}$ iterações de otimização
4. Constraint de recursos: execuções em simulador ou hardware real

## Hipóteses

### Hipótese Principal (H₀)

**Enunciado**: Existe um regime de ruído intermediário $p^* \in (0, p_{max})$ tal que a acurácia de generalização no conjunto de teste é estritamente maior do que no regime sem ruído ($p=0$).

**Formulação matemática**:
$$
\exists p^* > 0 : \mathbb{E}[\text{Acc}(f_{\theta^*,p^*}, \mathcal{D}_{test})] > \mathbb{E}[\text{Acc}(f_{\theta^*_0, 0}, \mathcal{D}_{test})]
$$

onde $\theta^*_0$ é o ótimo sem ruído e $\theta^*$ é o ótimo com ruído $p^*$.

**Critério de aceitação**: Diferença estatisticamente significativa (teste t pareado, $\alpha = 0.05$) e tamanho de efeito $d > 0.5$ (médio a grande).

### Hipóteses Secundárias

**H₁**: Ruído quântico atua como regularizador natural, reduzindo overfitting.

**Operacionalização**:
$$
\text{Gap}_{train-test}(p^*) < \text{Gap}_{train-test}(0)
$$
onde $\text{Gap}_{train-test}(p) = \text{Acc}_{train}(p) - \text{Acc}_{test}(p)$

---

**H₂**: Schedules dinâmicos de ruído (annealing) superam configurações estáticas.

**Operacionalização**:
$$
\mathbb{E}[\text{Acc}(\text{schedule\_dynamic})] > \mathbb{E}[\text{Acc}(\text{schedule\_constant})]
$$

**Schedules a testar**:
- Constant: $p(t) = p_0$
- Linear decay: $p(t) = p_0 \cdot (1 - t/T)$
- Cosine annealing: $p(t) = p_0 \cdot \frac{1 + \cos(\pi t/T)}{2}$

---

**H₃**: O efeito benéfico de ruído varia conforme características do dataset.

**Operacionalização**: ANOVA de dois fatores (tipo de ruído × dataset) com interação significativa ($p < 0.05$).

## Variáveis e Fatores

### Variáveis Independentes (Fatores)

| Fator | Tipo | Níveis | Justificativa |
|-------|------|--------|---------------|
| Dataset | Categórico | {moons, circles, blobs, iris} | Representam diferentes distribuições e separabilidades |
| Ansatz | Categórico | {BasicEntangler, StronglyEntangling, ...} | Arquiteturas com diferentes capacidades expressivas |
| Tipo de Ruído | Categórico | {Depolarizing, AmplitudeDamping, PhaseDamping, ...} | Modelos físicos de decoerência |
| Parâmetro de Ruído (p) | Contínuo | [0, 0.5] | Intensidade do canal |
| Schedule de Ruído | Categórico | {Constant, Linear, Cosine} | Estratégias de annealing |
| Inicialização | Categórico | {random, zeros, He, Xavier} | Pontos iniciais de otimização |
| Seed Aleatória | Discreto | {42, 123, 456, 789, 1024} | Controle de estocasticidade |

### Variáveis Dependentes (Métricas)

| Métrica | Definição | Fórmula | Range |
|---------|-----------|---------|-------|
| Acurácia | Proporção de predições corretas | $\frac{TP+TN}{TP+TN+FP+FN}$ | [0, 1] |
| Perda (Loss) | Entropia cruzada binária | $-\frac{1}{N}\sum [y \log \hat{y} + (1-y)\log(1-\hat{y})]$ | $[0, \infty)$ |
| F1-Score | Média harmônica de precisão e recall | $2 \cdot \frac{P \cdot R}{P + R}$ | [0, 1] |
| Gap Train-Test | Diferença de acurácia | $\text{Acc}_{train} - \text{Acc}_{test}$ | [-1, 1] |
| Épocas até Convergência | Iterações até critério de parada | Contagem | $\mathbb{N}$ |

## Espaço Experimental

### Total de Configurações

$$
N_{configs} = N_{datasets} \times N_{ansätze} \times N_{noise\_types} \times N_{schedules} \times N_{inits} \times N_{seeds}
$$

**Para este estudo**:
$$
N_{configs} = 4 \times 7 \times 6 \times 3 \times 8 \times 5 = 20.160 \text{ configurações}
$$

### Design Experimental

- **Tipo**: Full factorial design (todos os fatores cruzados)
- **Replicações**: 5 seeds aleatórias por configuração
- **Ordem de execução**: Randomizada para controlar efeitos de ordem
- **Controles**:
  - Baseline sem ruído ($p=0$) para cada configuração
  - Validação cruzada k-fold ($k=5$) para estimar generalização

## Critérios de Sucesso

### Critérios Estatísticos

1. **Significância**: $p$-value < 0.05 (com correção de Bonferroni para múltiplas comparações)
2. **Tamanho de Efeito**: Cohen's $d > 0.5$ (efeito médio ou grande)
3. **Intervalos de Confiança**: IC de 95% não incluem zero para diferenças

### Critérios Práticos

1. **Melhoria Mínima**: Acurácia com ruído ótimo > baseline + 2% (absoluto)
2. **Robustez**: Efeito benéfico presente em ≥75% dos datasets
3. **Eficiência**: Convergência mais rápida (redução de ≥10% em épocas)

## Ameaças à Validade (Antecipadas)

### Validade Interna
- **Confounding de hardware**: Simuladores vs hardware real podem ter comportamentos diferentes
- **Ordem de execução**: Não randomizar pode introduzir viés temporal

**Mitigação**: Seeds fixas, ordem randomizada, múltiplas replicações

### Validade Externa
- **Datasets sintéticos**: Resultados podem não generalizar para problemas reais
- **Tamanho de datasets**: Pequenos datasets (N < 1000) podem não refletir cenários práticos

**Mitigação**: Incluir pelo menos um dataset real (iris), discutir scope conditions

### Validade de Construto
- **Operacionalização de "benéfico"**: Definido apenas como acurácia, mas pode haver outros benefícios (velocidade, robustez)
- **Métricas**: Acurácia pode ser enganosa em datasets desbalanceados

**Mitigação**: Usar múltiplas métricas (F1, precision, recall), reportar distribuições de classes

### Validade Estatística
- **Múltiplas comparações**: 20.160 testes aumentam chance de falsos positivos
- **Pressupostos de testes**: Normalidade, homoscedasticidade podem ser violados

**Mitigação**: Correção de Bonferroni, testes não-paramétricos quando apropriado

## Scope Conditions

**Resultados deste estudo são válidos sob as seguintes condições**:

1. **Sistemas quânticos**: $n_{qubits} \leq 8$ (limitação de simuladores clássicos)
2. **Datasets**: Problemas de classificação binária ou multiclasse com $N \leq 1000$ amostras
3. **Ruído**: Modelos idealizados de Lindblad, sem crosstalk ou ruído correlacionado
4. **Otimização**: Gradiente estocástico com learning rate fixo ou schedule simples
5. **Hardware**: Simuladores perfeitos (AerSimulator, lightning.qubit)

**Não generaliza para**:
- Hardware quântico ruidoso real (NISQ) com ruído não-Markoviano
- Problemas de otimização combinatória ou sampling
- Circuitos muito profundos (depth > 50) com muitos qubits

## Referências do Problema

[Adicionar referências que fundamentam a formulação do problema]

1. Nielsen, M. A., & Chuang, I. L. (2010). *Quantum Computation and Quantum Information*. Cambridge University Press. (Definições de canais quânticos)

2. Farhi, E., & Neven, H. (2018). Classification with Quantum Neural Networks on Near Term Processors. arXiv:1802.06002. (VQCs)

3. McClean, J. R., et al. (2018). Barren plateaus in quantum neural network training landscapes. *Nature Communications*, 9, 4812. (Desafios de otimização)

---

**Última atualização**: [Data]  
**Revisor**: [Nome]  
**Status**: [Draft | Em Revisão | Aprovado]
