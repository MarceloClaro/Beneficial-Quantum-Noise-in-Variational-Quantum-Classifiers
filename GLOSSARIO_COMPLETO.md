# Glossário de Termos - Geração de Artigos Científicos QUALIS A1

**Versão:** 1.0  
**Data:** 26/12/2025  
**Conformidade:** QUALIS A1 - PR/Quantum


---


## A

### Auditoria Técnica
Processo sistemático de inventariar e verificar todos os componentes de um projeto de pesquisa (código, dados, configurações) para garantir que o que foi implementado corresponde exatamente ao que será reportado no artigo científico.

**Exemplo:** Verificar que cada valor numérico reportado no artigo pode ser rastreado até um arquivo específico, linha de código ou log de execução.


### Ansatz
Estrutura parametrizada de um circuito quântico variacional, definindo como portas quânticas são organizadas e conectadas. Exemplos incluem: Hardware Efficient, Strongly Entangling, Random Entangling.

**Referência no Código:** `framework_investigativo_completo.py:L450-550`


---


## B

### Baseline
Configuração de controle sem intervenção experimental. No contexto deste projeto, refere-se ao VQC operando sem ruído quântico (γ=0).

**Métrica:** Acurácia baseline ≈ 50% (chance level para classificação binária)


### Beneficial Noise (Ruído Benéfico)
Fenômeno onde a introdução controlada de ruído quântico melhora o desempenho de um classificador variacional, contrariando a intuição de que ruído sempre prejudica.

**Evidência Empírica:** Acurácia máxima de 65.83% com γ=0.001431 vs. 50% sem ruído.


---


## C

### Cohen's d
Tamanho de efeito padronizado que mede a diferença entre duas médias em unidades de desvio padrão. Valores: pequeno (0.2), médio (0.5), grande (0.8), muito grande (>1.2).

**Resultado Obtido:** d = 4.03 (muito grande, praticamente definitivo)


### Conivência Código-Texto
Correspondência perfeita (100%) entre o que está implementado no código e o que é reportado no texto do artigo. Inclui: valores numéricos, configurações, métodos, e conclusões.

**Ferramenta:** `relatorio_conivencia.md` documenta esta verificação.


### Crosstalk
Tipo de ruído quântico onde operações em um qubit afetam inadvertidamente qubits vizinhos, comum em arquiteturas supercondutoras.

**Implementação:** `framework_investigativo_completo.py:L280-295`


---


## D

### Depolarizing Noise
Canal de ruído quântico que substitui o estado com probabilidade p por estado completamente misturado. Modelado por operadores de Kraus.

**Equação:** $\rho' = (1-p)\rho + p\frac{I}{2}$


### Dynamic Schedule
Estratégia onde a intensidade do ruído varia ao longo do treinamento segundo uma função (Linear, Cosine, Exponential), análogo ao annealing térmico.

**Evidência:** Schedule Cosine reduz épocas de convergência em 12.6%.


---


## E

### Effect Sizes
Métricas que quantificam a magnitude de um efeito (ex: Cohen's d, Glass's Δ, Hedges' g), independentemente do tamanho da amostra. Essenciais para avaliação da significância prática (não apenas estatística).

**Uso:** Complementa p-valores para avaliar relevância científica.


### Entanglement
Correlação quântica não-clássica entre qubits. Ansätze com maior entanglement tendem a maior expressividade mas também maior susceptibilidade ao ruído.

**Medição:** Entropia de emaranhamento, negatividade.


---


## F

### Fidelidade (Fidelity)
Medida de proximidade entre dois estados quânticos: $F(\rho, \sigma) = \text{Tr}(\sqrt{\sqrt{\rho}\sigma\sqrt{\rho}})^2$. Valores entre 0 (ortogonais) e 1 (idênticos).

**Implicação:** Fidelidade < 0.99 indica degradação significativa por ruído.


---


## G

### Generalização
Capacidade do modelo de desempenhar bem em dados não vistos durante o treinamento. Medida pela diferença entre acurácia de treino e teste.

**Gap de Generalização:** $\Delta = \text{Acc}_{\text{train}} - \text{Acc}_{\text{test}}$


---


## H

### Hardware Efficient Ansatz
Arquitetura de circuito quântico otimizada para execução em hardware específico, minimizando operações de dois qubits (CNOT) que são mais ruidosas.

**Característica:** Menor profundidade, maior taxa de sucesso em NISQ devices.


---


## I

### [INFORMAÇÃO AUSENTE]
Marcador usado quando uma informação deveria existir mas não foi encontrada no código, dados ou logs. Diferente de [NÃO DISPONÍVEL].

**Exemplo:** "Versão exata da biblioteca scipy: [INFORMAÇÃO AUSENTE]"


### Intervalo de Confiança (IC)
Faixa de valores plausíveis para um parâmetro populacional, calculado a partir da amostra. IC 95% significa que 95% dos intervalos calculados conterão o valor verdadeiro.

**Formato:** Média ± Margem de Erro, ex: 65.83% ± 2.14%


---


## K

### Kraus Operators
Operadores matemáticos que descrevem a evolução de um sistema quântico aberto (com interação ambiental). Todo canal quântico pode ser representado por operadores de Kraus.

**Propriedade:** $\sum_i K_i^\dagger K_i = I$


---


## L

### [LACUNA DE CITAÇÃO]
Marcador usado (apenas em modo R0) quando falta uma referência para suportar uma afirmação, mas não é permitido adicionar novas referências.

**Uso:** "[LACUNA DE CITAÇÃO] - necessário referência sobre quantum advantage"


---


## M

### Manifesto de Execução
Arquivo JSON contendo metadados completos de reprodutibilidade: seeds, versões de bibliotecas, comandos executados, timestamps, configurações de hardware.

**Exemplo:** `resultados_multiframework_20251226_172214/execution_manifest.json`


---


## N

### NISQ (Noisy Intermediate-Scale Quantum)
Era atual da computação quântica (2018-2030) caracterizada por 50-1000 qubits com ruído significativo, sem correção de erro quântico completa.

**Definição Preskill (2018):** Quantum systems with 50-100 qubits, noise-limited.


### [NÃO DISPONÍVEL]
Marcador usado quando a informação não pode ser gerada (ex: resultados de um pipeline que não executa por limitações computacionais).

**Exemplo:** "Resultados em hardware IBM real: [NÃO DISPONÍVEL] - acesso não disponível"


---


## O

### Overfitting
Fenômeno onde modelo ajusta ruído dos dados de treino, prejudicando generalização. Detectado por gap entre acurácia de treino e teste.

**Mitigação:** Regularização, early stopping, validação cruzada.


---


## P

### Phase Damping
Canal de ruído que causa perda de coerência de fase sem alterar populações dos estados |0⟩ e |1⟩. Modelado por operadores de Kraus específicos.

**Resultado Chave:** Melhor tipo de ruído, +3.75% vs Depolarizing (p<0.05)


---


## Q

### Quality Gate
Ponto de verificação ao final de cada fase do processo de geração de artigo para garantir que os critérios de qualidade foram atendidos antes de prosseguir para a próxima fase.

**Exemplo:** Fase 1 - verificar que cada item tem origem rastreável (arquivo/função/linha).


### Qubit
Unidade básica de informação quântica. Pode estar em superposição de estados |0⟩ e |1⟩: $|\psi\rangle = \alpha|0\rangle + \beta|1\rangle$.

---


## R

### R0 (Reference Policy - Locked)
Política de referências onde o conjunto de citações é travado (não pode ser alterado). Usado quando lista de referências já foi aprovada.

**Ação:** Marcar lacunas com [LACUNA DE CITAÇÃO].


### R1 (Reference Policy - Expanded)
Política de referências onde novas citações podem ser buscadas e adicionadas durante o processo de escrita, seguindo 7 categorias predefinidas.

**Categorias:** Fundamentos teóricos, Estado da arte, Metodologia, Benchmarks, Frameworks, Aplicações, Surveys.


### Rastreabilidade
Capacidade de traçar cada afirmação, número ou resultado em um artigo de volta à sua origem exata no código, dados ou logs de execução.

**Formato Tabela:** Seção | Afirmação | Evidência (Arquivo:Linha) | Referência


### Reprodutibilidade
Capacidade de terceiros replicarem os resultados reportados seguindo os procedimentos descritos. Requer: seeds fixas, versões documentadas, ambiente especificado.

**Padrão Ouro:** Código + dados + ambiente + seeds → mesmos resultados (±ε estatístico)


---


## S

### Schedule
Função que controla como parâmetros variam durante o treinamento. No contexto deste projeto, refere-se à variação da intensidade de ruído γ(epoch).

**Tipos:** Static (constante), Linear (γ∝epoch), Cosine (γ∝cos(πepoch/T)), Exponential (γ∝e^(-epoch))


### Seeds
Valores iniciais para geradores de números pseudo-aleatórios. Essenciais para reprodutibilidade de experimentos estocásticos.

**Seeds Usadas:** [42, 43]


---


## T

### Threats to Validity (Ameaças à Validade)
Fatores que podem comprometer a validade das conclusões de um estudo. Dividem-se em:

- **Interna:** Causalidade (confounders, viés de seleção)
- **Externa:** Generalização (limitação de datasets, escala)
- **Construto:** Medição (operacionalização de conceitos)
- **Estatística:** Inferência (poder estatístico, correções)


**Exemplo:** Uso de simuladores (ameaça à validade externa - generalização para hardware real)


---


## V

### VQC (Variational Quantum Classifier)
Circuito quântico parametrizado usado para classificação de dados, otimizado via algoritmo híbrido clássico-quântico.

**Componentes:** Encoding de dados + Ansatz variacional + Medição


### Variational Algorithm
Algoritmo híbrido onde circuito quântico (quantum processing) é otimizado por otimizador clássico (classical processing), adequado para NISQ devices.

**Exemplos:** VQE, QAOA, VQC


---


## Siglas e Acrônimos

| Sigla | Significado | Contexto |
|-------|-------------|----------|
| ABNT | Associação Brasileira de Normas Técnicas | Formatação (MODE_B) |
| ANOVA | Analysis of Variance | Teste estatístico |
| CARS | Create A Research Space | Estrutura introdução |
| CNOT | Controlled-NOT | Porta quântica 2-qubit |
| DOI | Digital Object Identifier | Referências |
| IC | Intervalo de Confiança | Estatística |
| IMRAD | Intro-Methods-Results-And-Discussion | Estrutura artigo |
| LaTeX | Lamport's TeX | Formatação (MODE_A) |
| NISQ | Noisy Intermediate-Scale Quantum | Era quântica atual |
| npj QI | Nature Partner Journal Quantum Information | Periódico alvo |
| PR | Physical Review | Família periódicos |
| QML | Quantum Machine Learning | Área de pesquisa |
| SMART | Specific-Measurable-Achievable-Relevant-Time-bound | Objetivos |
| VQA | Variational Quantum Algorithm | Classe algoritmos |
| VQC | Variational Quantum Classifier | Modelo específico |
| VQE | Variational Quantum Eigensolver | Algoritmo variacional |

---


## Notação Matemática

| Símbolo | Significado | Uso |
|---------|-------------|-----|
| γ (gamma) | Intensidade do ruído | Parâmetro principal |
| ρ (rho) | Matriz densidade | Estado quântico |
| θ (theta) | Parâmetros variacionais | Vetor de pesos |
| ⟨·⟩ | Valor esperado | Medições quânticas |
| ⊗ | Produto tensorial | Sistemas compostos |
| † | Hermitiano conjugado | Operadores quânticos |
| ℋ | Espaço de Hilbert | Espaço de estados |
| 𝒩(·) | Canal de ruído | Mapa quântico |
| ℒ | Função de perda (loss) | Otimização |
| 𝒟 | Dataset | Conjunto de dados |

---


## Referências do Glossário

Este glossário foi compilado a partir de:

1. Nielsen & Chuang (2010) - Quantum Computation and Quantum Information
2. Preskill (2018) - Quantum Computing in the NISQ era
3. Cerezo et al. (2021) - Variational Quantum Algorithms
4. Du et al. (2021) - Beneficial quantum noise in VQE
5. Documentação PennyLane, Qiskit, Cirq
6. Framework código: `framework_investigativo_completo.py`


---


**Última atualização:** 26/12/2025  
**Versão:** 1.0  
**Status:** ✅ Completo

