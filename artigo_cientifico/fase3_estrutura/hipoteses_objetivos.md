# FASE 3.2: Estruturação de Hipóteses e Objetivos

**Data:** 25 de dezembro de 2025  
**Framework SMART:** Specific, Measurable, Achievable, Relevant, Time-bound

---

## 1. HIPÓTESE PRINCIPAL (H₀)

### Formulação

**Se** ruído quântico moderado for introduzido sistematicamente em Variational Quantum Classifiers através de schedules dinâmicos,  
**Então** a acurácia de generalização em dados de teste aumentará significativamente (Δ_acc > 5%),  
**Porque** o ruído atua como regularizador estocástico que previne overfitting e suaviza o landscape de otimização.

### Fundamentação Teórica

**Base Conceitual:**
1. **Regularização Estocástica (Bishop, 1995):** Treinar com ruído ≡ regularização de Tikhonov (L2)
2. **Ruído Benéfico em VQCs (Du et al., 2021):** Demonstração empírica em dataset Moons
3. **Ressonância Estocástica (Benzi et al., 1981):** Ruído pode amplificar sinais em sistemas não-lineares

### Variáveis

| Tipo | Variável | Operacionalização |
|------|----------|-------------------|
| **Independente** | Intensidade de Ruído (γ) | γ ∈ [10⁻⁵, 10⁻¹], escala logarítmica |
| **Independente** | Schedule Dinâmico | {Static, Linear, Exponential, Cosine} |
| **Dependente** | Acurácia de Generalização | Acc_test em dataset de teste (15% dados) |
| **Controle** | Baseline (sem ruído) | γ = 0 (ruído zero) |
| **Controle** | Ansatz, Dataset, Optimizer | Fixos para comparação fair |

### Predição Quantitativa

**Esperado:**
- Regime Ótimo: γ ∈ [10⁻³, 10⁻²]
- Melhoria vs. Baseline: Δ_acc = +5% a +15%
- Schedule Ótimo: Cosine ou Exponential > Static

**Critério de Confirmação:**
- H₀ é **confirmada** se: p < 0.05 (ANOVA) E Cohen's d > 0.5 (efeito médio ou grande)
- H₀ é **refutada** se: p ≥ 0.05 OU melhoria < 2% (insignificante)

---

## 2. HIPÓTESES DERIVADAS (H₁, H₂, H₃, H₄)

### H₁: Efeito do Tipo de Ruído Quântico

**Enunciado:**
> Diferentes modelos de ruído quântico (Depolarizing, Amplitude Damping, Phase Damping, Bit Flip, Phase Flip) produzirão efeitos significativamente distintos na acurácia de VQCs, com Phase Damping e Amplitude Damping demonstrando maior benefício (Δ_acc > 7%) comparado a Depolarizing (Δ_acc ≈ 5%).

**Dimensão da Lacuna:** Generalidade (Du et al. usaram apenas Depolarizing)

**Fundamentação:**
- **Wang et al. (2021):** Demonstram que diferentes ruídos afetam trainability diferentemente
- **Lindblad Formalism:** Cada ruído tem operadores de Kraus distintos, afetando dinâmica quântica de modo único
- **Intuição Física:** 
  - Phase Damping preserva populações → menos destrutivo
  - Amplitude Damping simula T₁ decay → mais realista
  - Depolarizing é modelo simplificado → menos refinado

**Predição Quantitativa:**
```
Phase Damping:      Δ_acc = +8% ± 2%  (esperado melhor)
Amplitude Damping:  Δ_acc = +7% ± 2%
Depolarizing:       Δ_acc = +5% ± 1%  (baseline de Du et al.)
Bit Flip:           Δ_acc = +3% ± 2%
Phase Flip:         Δ_acc = +4% ± 2%
```

**Métrica de Teste:**
- ANOVA one-way: NoiseType como fator
- Post-hoc: Tukey HSD para comparações múltiplas
- Effect Size: Cohen's d entre pares

**Critério de Confirmação:**
- p < 0.05 (efeito significativo de NoiseType)
- η² > 0.06 (efeito médio ou grande)
- Ranking confirmado: Phase/Amplitude > Depolarizing > Bit/Phase

---

### H₂: Curva de Dose-Resposta (Regime Ótimo de Ruído)

**Enunciado:**
> A relação entre intensidade de ruído (γ) e acurácia segue curva não-monotônica (inverted-U), com regime ótimo em γ_opt ∈ [10⁻³, 5×10⁻³], onde acurácia é maximizada. Fora deste regime, ruído excessivo (γ > 10⁻²) degrada performance abaixo de baseline, e ruído insuficiente (γ < 10⁻⁴) não produz benefício.

**Dimensão da Lacuna:** Generalidade (mapeamento sistemático de γ)

**Fundamentação:**
- **Du et al. (2021):** Identificaram regime benéfico, mas não mapearam curva completa
- **Teoria de Regularização:** Regularização ótima equilibra bias e variance
- **Ressonância Estocástica:** Ruído ótimo amplifica sinal, excesso mascara

**Predição Quantitativa:**
```
γ < 10⁻⁴:           Δ_acc ≈ 0%        (sem efeito)
γ ∈ [10⁻³, 5×10⁻³]: Δ_acc = +10% ± 3% (regime ótimo) ✨
γ ∈ [10⁻², 5×10⁻²]: Δ_acc = -5% ± 5%  (degradação)
γ > 10⁻¹:           Δ_acc = -20% ± 5% (colapso)
```

**Métrica de Teste:**
- Polynomial Regression (grau 2 ou 3) para modelar curva
- Identificação de máximo: γ_opt = arg max Acc(γ)
- Intervalo de Confiança para γ_opt (bootstrap)

**Critério de Confirmação:**
- Curva inverted-U confirmada (R² > 0.7)
- γ_opt dentro do intervalo previsto
- Acc(γ_opt) - Acc(0) > 5% (significativo)

---

### H₃: Interação Ansatz × Tipo de Ruído

**Enunciado:**
> Existe interação significativa (p < 0.05) entre arquitetura de ansatz e tipo de ruído quântico, tal que ansätze com maior expressividade (StronglyEntangling, RandomLayers) demonstram maior benefício de ruído (Δ_acc > 12%) comparado a ansätze simples (BasicEntangling: Δ_acc ≈ 5%), devido ao efeito de mitigação de barren plateaus.

**Dimensão da Lacuna:** Interação Multi-Fatorial (não investigada sistematicamente)

**Fundamentação:**
- **Holmes et al. (2022):** Trade-off expressividade × trainability
- **Choi et al. (2022):** Ruído mitiga barren plateaus
- **Lógica:** Ansätze expressivos sofrem mais de barren plateaus → mais benefício de ruído

**Predição Quantitativa:**
```
StronglyEntangling + Phase Damping:  Δ_acc = +12% ± 3% (maior interação)
BasicEntangling + Phase Damping:     Δ_acc = +5% ± 2%  (menor interação)
Diferença (Interação):               ΔΔ_acc = +7%      (efeito de interação)
```

**Métrica de Teste:**
- ANOVA Two-Way: Ansatz × NoiseType
- Termo de Interação: F-statistic e p-value
- Interaction Plot: Visualização de linhas não-paralelas
- Effect Size: η²_partial para interação

**Critério de Confirmação:**
- p_interaction < 0.05
- η²_interaction > 0.04 (efeito pequeno-médio)
- Interaction plot mostra linhas não-paralelas (visual)

---

### H₄: Superioridade de Schedules Dinâmicos

**Enunciado:**
> Schedules dinâmicos de ruído (Linear, Exponential, Cosine) superam estratégia estática (Static) em termos de acurácia final (Δ_acc > 3% vs. Static) e velocidade de convergência (épocas até convergência reduzida em 20-30%), com Cosine Annealing demonstrando melhor performance global devido ao equilíbrio entre exploração inicial e refinamento final.

**Dimensão da Lacuna:** Dinâmica (inovação metodológica original)

**Fundamentação:**
- **Kirkpatrick et al. (1983):** Simulated Annealing - temperatura decrescente melhora convergência
- **Loshchilov & Hutter (2016):** Cosine Annealing em deep learning - sucesso empírico
- **Analogia:** Ruído alto (início) → exploração; Ruído baixo (final) → refinamento

**Predição Quantitativa:**
```
Cosine Annealing:    Δ_acc = +11% ± 2% (melhor)
Exponential:         Δ_acc = +9% ± 2%
Linear:              Δ_acc = +8% ± 2%
Static (γ_opt):      Δ_acc = +7% ± 2%  (baseline)
```

**Métrica de Teste:**
- ANOVA one-way: Schedule como fator
- Post-hoc: Tukey HSD
- Convergence Analysis: Número de épocas até tolerance (δ_loss < 10⁻⁵)

**Critério de Confirmação:**
- p < 0.05 (efeito significativo de Schedule)
- Cosine > Static (p_pairwise < 0.05, d > 0.5)
- Redução em épocas: 20-30% (quantificável)

**Nota Importante:**
Se H₄ for **parcialmente refutada** (Schedules não superam Static significativamente), isso também é resultado científico válido e deve ser reportado honestamente.

---

## 3. OBJETIVOS ESPECÍFICOS (Framework SMART)

### Objetivo 1: Quantificar Benefício de Ruído em Múltiplos Datasets

**S - Specific:**
> Quantificar a melhoria de acurácia proporcionada por ruído quântico ótimo (γ_opt) em comparação com baseline (γ=0) para 4 datasets distintos (Moons, Circles, Iris, Wine).

**M - Measurable:**
- Métrica: Δ_acc = Acc(γ_opt) - Acc(0)
- Target: Δ_acc > 5% em ≥3 datasets
- IC 95%: Reportar intervalo de confiança para cada dataset

**A - Achievable:**
- Du et al. (2021) alcançaram ~5% em Moons
- Expectativa: Generalização para outros datasets é viável
- Recursos: Framework implementado, hardware suficiente

**R - Relevant:**
- Testa H₁ (generalidade)
- Lacuna: Du et al. usaram 1 dataset apenas
- Importância: Demonstra generalidade do fenômeno

**T - Time-bound:**
- Grid search: ~20-30 horas por dataset (total: 80-120h)
- Bayesian optimization: ~4-6 horas por dataset (total: 16-24h)

**Alinhamento com Lacuna:** Dimensão 1 (Generalidade)

---

### Objetivo 2: Mapear Curva Dose-Resposta de Ruído

**S - Specific:**
> Mapear sistematicamente a relação γ → Acc(γ) para 11 intensidades de ruído logaritmicamente espaçadas (10⁻⁵ a 10⁻¹) e identificar regime ótimo (γ_opt) com IC 95%.

**M - Measurable:**
- Métrica: Curva Acc(γ) + γ_opt + IC_95%(γ_opt)
- Target: γ_opt ∈ [10⁻³, 10⁻²] (hipótese)
- Fit: R² > 0.7 para polynomial regression

**A - Achievable:**
- 11 valores de γ × 5 repetições = 55 experimentos por configuração
- Viável com otimização Bayesiana (subset inteligente)

**R - Relevant:**
- Testa H₂ (curva dose-resposta)
- Lacuna: Mapeamento completo não foi feito
- Aplicação: Permite engenharia de ruído ótimo

**T - Time-bound:**
- Por configuração (ansatz + noise type + dataset): ~2-3 horas
- Total para 10 configurações representativas: ~20-30 horas

**Alinhamento com Lacuna:** Dimensão 1 (Generalidade) + Aplicação Prática

---

### Objetivo 3: Investigar Interações Multi-Fatoriais

**S - Specific:**
> Conduzir ANOVA multifatorial (7 fatores: Ansatz, NoiseType, NoiseStrength, Schedule, Dataset, InitStrategy, CircuitDepth) para identificar fatores significativos (p < 0.05) e interações de ordem 2 (ex: Ansatz × NoiseType).

**M - Measurable:**
- Métrica: F-statistic, p-value, η² para cada fator e interação
- Target: ≥3 interações significativas (p < 0.05, η² > 0.04)
- Output: Tabela ANOVA completa

**A - Achievable:**
- Dados de 8.280 experimentos (amostra grande)
- Statsmodels implementa ANOVA multifatorial
- Computacionalmente viável (análise post-hoc em minutos)

**R - Relevant:**
- Testa H₃ (interação Ansatz × Ruído)
- Lacuna: Análise de interações ausente na literatura
- Rigor: Padrão QUALIS A1

**T - Time-bound:**
- Coleta de dados: Já concluída (8.280 experimentos)
- Análise estatística: ~2-4 horas (ANOVA + post-hoc)

**Alinhamento com Lacuna:** Dimensão 3 (Interação Multi-Fatorial)

---

### Objetivo 4: Validar Superioridade de Schedules Dinâmicos

**S - Specific:**
> Comparar 4 estratégias de schedule (Static, Linear, Exponential, Cosine) em termos de acurácia final e convergência, usando γ_opt identificado no Objetivo 2.

**M - Measurable:**
- Métrica 1: Δ_acc = Acc(Schedule) - Acc(Static)
- Métrica 2: Épocas até convergência (δ_loss < 10⁻⁵)
- Target: Cosine > Static (p < 0.05, d > 0.5)

**A - Achievable:**
- 4 schedules × 5 repetições × 10 configurações = 200 experimentos
- Tempo: ~10-15 horas (parallelizable)

**R - Relevant:**
- Testa H₄ (superioridade de schedules dinâmicos)
- Lacuna: Inovação metodológica original (não existe na literatura)
- Contribuição: Se validado, estabelece novo paradigma

**T - Time-bound:**
- Experimentos: ~10-15 horas
- Análise: ~2 horas

**Alinhamento com Lacuna:** Dimensão 2 (Dinâmica) - **INOVAÇÃO METODOLÓGICA** ✨

---

## 4. TABELA DE ALINHAMENTO: Lacuna → Hipótese → Objetivo → Métrica

| Lacuna | Hipótese | Objetivo | Métrica | Critério de Sucesso |
|--------|----------|----------|---------|---------------------|
| **Generalidade (Dataset)** | H₁ (Tipo de Ruído) | Obj. 1 (Múltiplos Datasets) | Δ_acc por dataset | Δ_acc > 5% em ≥3 datasets, p < 0.05 |
| **Generalidade (Ruído)** | H₂ (Dose-Resposta) | Obj. 2 (Mapear γ → Acc) | γ_opt + IC 95% | γ_opt ∈ [10⁻³, 10⁻²], R² > 0.7 |
| **Interação Multi-Fatorial** | H₃ (Ansatz × Ruído) | Obj. 3 (ANOVA) | F, p, η² de interação | p_interaction < 0.05, η² > 0.04 |
| **Dinâmica (Schedules)** | H₄ (Dinâmico > Estático) | Obj. 4 (Validar Schedules) | Δ_acc, Épocas | Cosine > Static (p < 0.05, d > 0.5) |

---

## 5. CHECKLIST SMART PARA CADA OBJETIVO

### Objetivo 1
- [x] **S**pecific: Quantificar Δ_acc em 4 datasets ✅
- [x] **M**easurable: Δ_acc > 5%, IC 95% ✅
- [x] **A**chievable: Framework pronto, 16-24h ✅
- [x] **R**elevant: Testa generalidade (Lacuna 1) ✅
- [x] **T**ime-bound: 16-24 horas ✅

### Objetivo 2
- [x] **S**pecific: Mapear 11 valores de γ, identificar γ_opt ✅
- [x] **M**easurable: γ_opt + IC 95%, R² > 0.7 ✅
- [x] **A**chievable: 20-30 horas ✅
- [x] **R**elevant: Engenharia de ruído ótimo ✅
- [x] **T**ime-bound: 20-30 horas ✅

### Objetivo 3
- [x] **S**pecific: ANOVA multifatorial, 7 fatores ✅
- [x] **M**easurable: F, p, η² ✅
- [x] **A**chievable: Dados prontos, análise 2-4h ✅
- [x] **R**elevant: Rigor QUALIS A1 (Lacuna 3) ✅
- [x] **T**ime-bound: 2-4 horas ✅

### Objetivo 4
- [x] **S**pecific: Comparar 4 schedules ✅
- [x] **M**easurable: Δ_acc, épocas, p, d ✅
- [x] **A**chievable: 200 experimentos, 10-15h ✅
- [x] **R**elevant: Inovação metodológica (Lacuna 2) ✅
- [x] **T**ime-bound: 10-15 horas ✅

---

## 6. PREVISÃO DE RESULTADOS E CENÁRIOS

### Cenário Ótimo (Todas Hipóteses Confirmadas)

✅ **H₁ Confirmada:** Diferentes ruídos têm efeitos distintos (p < 0.001)  
✅ **H₂ Confirmada:** Curva inverted-U clara, γ_opt ≈ 0.003  
✅ **H₃ Confirmada:** Interação Ansatz × Ruído significativa (p < 0.01)  
✅ **H₄ Confirmada:** Cosine > Static (+3%, p < 0.05)  

**Implicação:** Fenômeno de ruído benéfico é robusto, generalizado, e otimizável via schedules dinâmicos. Contribuição científica de alto impacto.

### Cenário Realista (3 de 4 Confirmadas)

✅ **H₁ Confirmada**  
✅ **H₂ Confirmada**  
✅ **H₃ Confirmada**  
❌ **H₄ Refutada:** Schedules dinâmicos não superam Static significativamente  

**Implicação:** Fenômeno de ruído benéfico é robusto e generalizado. Schedules dinâmicos não oferecem vantagem clara (resultado negativo, mas cientificamente válido). Discussão: Por que schedules não ajudam? Possível explicação: Ruído ótimo é estável durante treinamento.

### Cenário Parcial (2 de 4 Confirmadas)

✅ **H₁ Confirmada**  
✅ **H₂ Confirmada**  
❌ **H₃ Refutada:** Interação Ansatz × Ruído não significativa  
❌ **H₄ Refutada**  

**Implicação:** Fenômeno de ruído benéfico existe, mas é mais simples do que hipotetizado (efeitos aditivos, não interações). Schedules não ajudam. Contribuição: Simplificação do modelo conceitual.

### Cenário Crítico (≤1 Confirmada)

❌ **Maioria das Hipóteses Refutadas**  

**Implicação:** Fenômeno de ruído benéfico é específico de contexto (Du et al. 2021 foi caso especial). Resultado negativo importante: Delimita fronteiras de validade. Publicável em periódicos que valorizam resultados negativos (e.g., Scientific Reports, PLOS ONE).

---

## 7. CRITÉRIOS DE DECISÃO PÓS-ANÁLISE

### Se H₀ (Principal) for Confirmada:
✅ Submeter a periódico de alto impacto (Nature Comms, npj QI, Quantum)  
✅ Enfatizar generalização do fenômeno de Du et al. (2021)  
✅ Destacar rigor metodológico e inovação (schedules)  

### Se H₀ for Refutada, mas H₁-H₃ Confirmadas:
✅ Submeter a periódico especializado (PRX Quantum, Quantum Sci. Technol.)  
✅ Focar em caracterização detalhada do fenômeno  
✅ Discussão honesta sobre limitações  

### Se Maioria das Hipóteses for Refutada:
✅ Submeter como "Negative Results" (Scientific Reports, PLOS ONE)  
✅ Enfatizar importância de delimitar fronteiras de validade  
✅ Contribuição: Evitar que outros percam tempo em direções infrutíferas  

---

**Documento gerado automaticamente pelo framework de análise QUALIS A1**  
**Última atualização:** 25/12/2025
