# Guia Completo de Geração de Artigos Científicos QUALIS A1

**Framework de Rastreabilidade Total para Periódicos de Alto Impacto**  
**Versão:** 2.0  
**Data:** 26/12/2025  
**Conformidade:** Nature, Science, Quantum, Physical Review, npj QI

---

## 📖 SUMÁRIO EXECUTIVO

Este guia implementa o mega-prompt completo para geração de artigos científicos com 100% de conivência entre código/dados e texto, garantindo reprodutibilidade, auditabilidade e máxima avaliação por bancas de revisão QUALIS A1.

**Documentos Principais:**
- `config_artigo.json` - Configuração inicial
- `GLOSSARIO_COMPLETO.md` - Termos técnicos
- `FAQ_TROUBLESHOOTING_COMPLETO.md` - Perguntas frequentes
- `CHECKLIST_AUDITORIA_COMPLETO.md` - Sistema de pontuação 0-100
- `FLUXOGRAMA_R0_R1.md` - Políticas de referências
- Este documento (GUIA_COMPLETO_GERACAO_ARTIGOS.md)

---

## 🎯 VISÃO GERAL

### Objetivo Geral

Gerar um artigo científico completo, rigoroso e auditável, pronto para submissão a periódicos de alto impacto (Nature, Science, Quantum, Physical Review, QUALIS A1), com 100% de conivência entre código/dados e texto.

### Princípios Fundamentais

1. **NÃO inventar detalhes:** Se algo não estiver em código/dados/logs, usar **[INFORMAÇÃO AUSENTE]**
2. **NÃO inventar números:** Todo valor quantitativo deve ter lastro verificável; caso contrário **[NÃO DISPONÍVEL]**
3. **Se R0:** Proibido alterar conjunto de referências; quando faltar base, usar **[LACUNA DE CITAÇÃO]**
4. **Reprodutibilidade:** Reportar HW/SW, versões, seeds, configs, scripts e comandos
5. **Auditoria:** Cada seção exige rastreabilidade: **Seção → Evidência → Origem**

---

## 📁 ESTRUTURA DO PROCESSO (6 FASES)

```
artigo_cientifico/
├── config_artigo.json              ← Configuração inicial
│
├── fase1_analise/                   ← Auditoria técnica
│   ├── analise_codigo_inicial.md
│   ├── tabela_componentes.md
│   └── mapa_execucao.md
│
├── fase2_bibliografia/              ← Literatura e referências
│   ├── referencias_compiladas.md
│   ├── sintese_literatura.md
│   └── taxonomia_estado_da_arte.md
│
├── fase3_estrutura/                 ← Projeto do artigo
│   ├── problema_formal.md
│   ├── titulos_palavras_chave.md
│   └── hipoteses_objetivos.md
│
├── fase4_secoes/                    ← Redação principal
│   ├── resumo_abstract.md
│   ├── introducao_completa.md
│   ├── revisao_literatura_completa.md
│   ├── metodologia_completa.md
│   ├── resultados_completo.md
│   ├── discussao_completa.md
│   ├── conclusao_completa.md
│   └── agradecimentos_referencias.md
│
├── fase5_suplementar/               ← Material suplementar
│   ├── tabelas_suplementares.md
│   ├── figuras_suplementares.md
│   └── notas_metodologicas_adicionais.md
│
└── fase6_consolidacao/              ← Auditoria final
    ├── relatorio_conivencia.md
    ├── rastreabilidade_completa.md
    ├── tabela_codigo_metodo.md
    ├── artigo_completo_final.md
    └── sumario_executivo_final.md
```

---

## 🚀 INÍCIO RÁPIDO

### Passo 1: Configuração Inicial

```bash
# 1. Clonar repositório (se necessário)
git clone https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers.git
cd Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers

# 2. Copiar e editar configuração
cp config_artigo.json config_artigo_custom.json
nano config_artigo_custom.json
```

### Passo 2: Preencher config_artigo.json

```json
{
  "output_mode": "MODE_A",           // MODE_A (inglês/LaTeX) ou MODE_B (português/ABNT)
  "reference_policy": "R1",          // R0 (travadas) ou R1 (expandidas)
  "editorial_profile": "PROFILE_PR_QUANTUM",  // PROFILE_PR_QUANTUM ou PROFILE_GENERAL
  "target_journals": {
    "primary": "Quantum",
    "secondary": ["Physical Review A", "npj Quantum Information"]
  },
  "inputs": {
    "code_path": "/path/to/your/code",
    "data_path": "/path/to/your/data",
    "artifacts_path": "/path/to/figures_logs"
  },
  "user_inputs": {
    "research_question": "Sua pergunta de pesquisa aqui",
    "hypotheses": ["H₀: ...", "H₁: ..."],
    "approved_references": []  // Apenas para R0
  }
}
```

### Passo 3: Executar o Framework

```bash
# Opção A: Executar todas as fases sequencialmente
python gerador_artigo_completo.py --config config_artigo_custom.json

# Opção B: Executar fase por fase (recomendado)
python gerador_artigo_completo.py --config config_artigo_custom.json --fase 1
# ... revisar saída ...
python gerador_artigo_completo.py --config config_artigo_custom.json --fase 2
# ... e assim por diante ...
```

---

## 📝 FASE 1: AUDITORIA TÉCNICA

### Objetivo

Produzir um inventário técnico completo e verificável do projeto.

### Outputs Obrigatórios

1. **analise_codigo_inicial.md** (15-25 KB)
   - Estrutura técnica (arquivos, linhas, módulos, classes)
   - Componentes experimentais (fatores, níveis, total configs)
   - Metodologia implementada
   - Inovações e contribuições

2. **tabela_componentes.md** (5-10 KB)
   - Resumo executivo técnico em tabelas

3. **mapa_execucao.md** (8-12 KB)
   - Passo a passo reprodutível do pipeline

4. **manifesto_execucao.json** (auto-gerado)
   - Lista de comandos, seeds, paths, configs

### Template: analise_codigo_inicial.md

```markdown
# Análise de Código Inicial

## 1. Estrutura Técnica
- **Contagem:** X arquivos, Y linhas, Z diretórios
- **Dependências:** [Listar com versões]
  ```
  PennyLane==0.38.0
  Qiskit==1.0.2
  NumPy==1.24.3
  ```
- **Módulos:** [Tabela: arquivo | objetivo | I/O | deps]
- **Classes:** [Tabela: nome | propósito | métodos | instâncias]
- **Modelos:** [Lista: nome canônico | alias no código]
- **Técnicas Analíticas:** [Lista: técnica | referência no código]

## 2. Componentes Experimentais
- **Fatores:** [Tabela: fator | definição operacional]
  | Fator | Níveis | Definição |
  |-------|--------|-----------|
  | Ansatz | 7 | Arquitetura do circuito (Basic, Strongly Entangling, ...) |
  | Ruído | 6 | Tipo de canal (Depolarizing, Phase Damping, ...) |
  
- **Total de Configurações:** 4×7×6×2×8 = 2.688
- **Datasets:** [Tabela: nome | tamanho | features | pré-proc | licença]
- **Métricas:** [Tabela: métrica | definição formal | função no código]

## 3. Metodologia Implementada
- **Pré-processamento:** StandardScaler, PCA (se usado)
- **Treino/Otimização:** Adam, lr=0.01, 50 épocas, early stopping patience=5
- **Validação:** 80/20 train/test split, seeds=[42,43]
- **Plano Estatístico:** ANOVA 3-way, Tukey HSD post-hoc, IC 95%, Cohen's d

## 4. Inovações e Contribuições
- **Código Novo:** Dynamic noise schedules (Linear, Cosine) - primeira vez na literatura
- **Diferenças vs Baselines:** Du et al. (2021) usaram apenas 1 dataset, 1 ansatz, 1 ruído; 
  nós generalizamos para 4×7×6=168 combinações
```

### Quality Gate F1

- [ ] Cada item tem origem (arquivo/função/linha/config/log)
- [ ] O total de configurações foi calculado e conferido
- [ ] Ausências explicitadas com [INFORMAÇÃO AUSENTE]/[NÃO DISPONÍVEL]
- [ ] Pelo menos 3 revisores verificaram independentemente

### Comandos Úteis

```bash
# Contar arquivos e linhas
find . -name "*.py" | wc -l
find . -name "*.py" -exec wc -l {} + | tail -1

# Listar dependências
pip freeze > requirements_atual.txt

# Extrair classes e funções
grep -r "^class " *.py
grep -r "^def " *.py

# Calcular total de configs
python -c "print(4*7*6*2*8)"  # datasets × ansatz × ruído × schedules × inits
```

**Tempo Estimado:** 8-12 horas

---

## 📚 FASE 2: BIBLIOGRAFIA

### Objetivo

Compilar e sintetizar a literatura relevante seguindo política R0 ou R1.

### Outputs Obrigatórios

1. **referencias_compiladas.md** (20-30 KB)
   - 35-60 referências organizadas em 7 categorias (R1) ou lista aprovada (R0)
   - Cada referência com: Autor, Ano, Título, Periódico, DOI, Justificativa

2. **sintese_literatura.md** (15-25 KB)
   - Consensos, divergências, lacunas, posicionamento
   - Análise crítica (não lista de resumos!)

3. **taxonomia_estado_da_arte.md** (10-15 KB)
   - Tabela: abordagem × pressupostos × custo × falhas × evidência

### Fluxograma R0 vs R1

**Ver:** `FLUXOGRAMA_R0_R1.md` para detalhes completos

**Resumo:**
- **R0:** Lista pré-aprovada, marcar [LACUNA DE CITAÇÃO] se faltar
- **R1:** Buscar em 7 categorias, adicionar novas referências com DOI

### Template: referencias_compiladas.md (R1)

```markdown
# Referências Compiladas (Política R1)

## Categoria 1: Referências Fundacionais (5-8 refs)

### [F1] Nielsen & Chuang (2010) - Quantum Computation
**Referência Completa (ABNT):**
NIELSEN, M. A.; CHUANG, I. L. **Quantum computation and quantum information**. 
10th ed. Cambridge: Cambridge University Press, 2010.

**DOI:** 10.1017/CBO9780511976667  
**Citações:** 45,000+  
**Justificativa:** Fundamento teórico de computação quântica, define conceitos de 
qubits, portas, entrelaçamento usados na metodologia.

## Categoria 2: Estado da Arte (8-12 refs, últimos 2-3 anos)

### [A1] Du et al. (2021) - Beneficial Noise
**Referência Completa (ABNT):**
DU, Yuxuan et al. **Quantum noise can help quantum sensing**. 
*Physical Review Letters*, v. 128, n. 8, p. 080506, 2021. 
DOI: https://doi.org/10.1103/PhysRevLett.128.080506

**Citações:** 47  
**Relevância:** ⭐⭐⭐⭐⭐  
**Justificativa:** Primeiro trabalho a demonstrar ruído benéfico em contexto quântico, 
fundamento teórico para nossa hipótese principal.

**Achados-Chave:**
- Ruído pode melhorar sensibilidade de sensores quânticos
- Regime ótimo identificado teoricamente
- Validado em simulações VQE

**Contraponto:** 
Trabalho limitado a 1 dataset, 1 tipo de ruído. Nossa contribuição: 
generalização sistemática.

[... continuar para 7 categorias ...]

## Resumo Estatístico
- **Total de referências:** 45
- **Cobertura DOI:** 38/45 (84.4%)
- **Distribuição temporal:** 
  - 2019-2024: 28 (62%)
  - 2010-2018: 12 (27%)
  - <2010: 5 (11%)
```

### Quality Gate F2

- [ ] Cada técnica central do pipeline tem referência ou [LACUNA]
- [ ] Síntese contém contraste/avaliação (não lista de resumos)
- [ ] Pelo menos 5 referências com contrapontos identificados
- [ ] DOI presente em ≥80% das referências
- [ ] Referências organizadas por categoria (R1) ou verificadas na lista (R0)

**Tempo Estimado:** 6-10 horas (R0) ou 15-25 horas (R1)

---

## 🎯 FASE 3: PROJETO DO ARTIGO

### Objetivo

Estruturar o artigo com título, hipóteses e objetivos testáveis.

### Outputs Obrigatórios

1. **problema_formal.md** (3-5 KB)
   - Definição matemática formal do problema
   - Notação LaTeX rigorosa

2. **titulos_palavras_chave.md** (2-3 KB)
   - 3 opções de título (A1-compatível)
   - 6 palavras-chave

3. **hipoteses_objetivos.md** (5-8 KB)
   - H₀ + H₁-H₄ formalmente definidas
   - 4 objetivos SMART
   - Tabela evidência/teste para cada hipótese

### Template: problema_formal.md

```markdown
# Formal Problem Statement

## Definição Matemática

**Seja:**
- $\mathcal{D} = \{(x_i, y_i)\}_{i=1}^N$ um dataset com $N$ amostras, 
  $x_i \in \mathbb{R}^d$, $y_i \in \{0,1\}$
- $U(\theta)$ um circuito quântico parametrizado por $\theta \in \mathbb{R}^P$
- $\mathcal{N}_\gamma(\cdot)$ um canal de ruído quântico com intensidade $\gamma \in [0, \gamma_{max}]$
- $f(x; \theta, \gamma) = \langle 0^{\otimes n} | U^\dagger(\theta) \mathcal{N}_\gamma(|x\rangle) U(\theta) | 0^{\otimes n} \rangle$
- $\mathcal{L}(\theta, \gamma) = -\sum_{i=1}^N [y_i \log f(x_i; \theta, \gamma) + (1-y_i)\log(1-f(x_i; \theta, \gamma))]$

**O problema é encontrar:**

$$
(\theta^*, \gamma^*) = \arg\min_{\theta, \gamma} \mathbb{E}_{(x,y) \sim \mathcal{D}_{test}} [\mathcal{L}(y, f(x; \theta, \gamma))]
$$

**Sujeito a:**
- $\gamma \in [0, \gamma_{max}]$ onde $\gamma_{max}$ é limitado por degradação de fidelidade
- Arquitetura de $U(\theta)$ fixa (ansatz pré-definido)
- $P$ parâmetros treináveis

## Hipótese Principal (H₀)

**Formal:**
$$
\exists \gamma^* > 0 : \mathbb{E}[\text{Acc}_{test}(\gamma^*)] > \mathbb{E}[\text{Acc}_{test}(0)]
$$

**Em palavras:**
Existe um nível de ruído quântico $\gamma^* > 0$ onde a acurácia de generalização 
do VQC supera significativamente a acurácia sem ruído.

## Restrições e Escopo

1. **Escala:** Limitado a $n \leq 10$ qubits (hardware NISQ disponível)
2. **Ruídos:** Modelados por canais de Lindblad (Markovianos)
3. **Datasets:** Toy problems de classificação binária (generalização para real a validar)
4. **Simuladores:** PennyLane/Qiskit/Cirq (validação em hardware real futura)
```

### Template: hipoteses_objetivos.md

```markdown
# Hipóteses e Objetivos

## Hipóteses Formais

### H₀ (Principal): Existência de Regime Benéfico
**Afirmação:** Existe $\gamma^* > 0$ tal que VQC com ruído supera VQC sem ruído.  
**Teste:** t-test pareado entre acurácias com $\gamma^*$ vs $\gamma=0$  
**Métrica:** Acurácia de teste, IC 95%  
**Evidência Esperada:** Diferença > 5 pontos percentuais, p < 0.01, Cohen's d > 0.5  
**Origem:** H₀ testável via experimentos controlados

### H₁: Ruído como Regularizador
**Afirmação:** Ruído reduz overfitting (gap treino-teste)  
**Teste:** Comparar $\Delta_{train-test}$ com vs sem ruído  
**Métrica:** $|\text{Acc}_{train} - \text{Acc}_{test}|$  
**Evidência Esperada:** Gap menor com ruído (p < 0.05)  
**Origem:** Analogia com dropout (Srivastava et al., 2014)

### H₂: Superioridade Phase Damping
**Afirmação:** Phase Damping > Depolarizing  
**Teste:** ANOVA + Tukey HSD post-hoc  
**Métrica:** Acurácia média por tipo de ruído  
**Evidência Esperada:** Diferença > 2%, p < 0.05  
**Origem:** Intuição física (preserva populações)

### H₃: Vantagem de Schedules Dinâmicos
**Afirmação:** Schedules Cosine/Linear > Static  
**Teste:** ANOVA fatorial 2×K (schedule × outro fator)  
**Métrica:** Épocas até convergência, acurácia final  
**Evidência Esperada:** Redução ≥10% épocas, p < 0.05  
**Origem:** Analogia com simulated annealing

### H₄: Independência de Plataforma
**Afirmação:** Efeito benéfico em PennyLane, Qiskit, Cirq  
**Teste:** Teste de Friedman (medidas repetidas)  
**Métrica:** Acurácia em mesma config, 3 frameworks  
**Evidência Esperada:** p < 0.001 (efeito consistente)  
**Origem:** Validação de generalidade

## Objetivos SMART

### O1: Identificar Regime Ótimo de Ruído
**S (Specific):** Determinar $\gamma_{opt}$ para cada combinação (ansatz, dataset, ruído)  
**M (Measurable):** Curva acurácia vs γ, identificar máximo  
**A (Achievable):** Grid search em [0, 0.01] com 20 pontos  
**R (Relevant):** Responde H₀  
**T (Time-bound):** Fase de experimentos (semanas 3-6)

### O2: Quantificar Efeito de Regularização
**S:** Calcular $\Delta_{gap} = gap_{sem\_ruido} - gap_{com\_ruido}$  
**M:** Redução percentual do gap  
**A:** Comparar 100 configs com vs sem ruído  
**R:** Responde H₁  
**T:** Análise (semana 7)

### O3: Comparar Tipos de Ruído
**S:** Ranking dos 6 tipos por acurácia média  
**M:** Tabela com médias ± IC 95%, Cohen's d entre pares  
**A:** ANOVA multifatorial + post-hoc  
**R:** Responde H₂  
**T:** Análise (semana 7)

### O4: Validar em Múltiplos Frameworks
**S:** Replicar top-5 configs em PennyLane, Qiskit, Cirq  
**M:** Acurácia, tempo de execução, desvio entre frameworks  
**A:** Execução controlada (mesmos seeds)  
**R:** Responde H₄  
**T:** Validação (semana 8)
```

### Quality Gate F3

- [ ] Problema formal compatível com execução real do código
- [ ] Cada hipótese tem teste/métrica correspondente no pipeline
- [ ] Objetivos seguem critério SMART
- [ ] Notação matemática consistente e definida
- [ ] Evidências esperadas são quantitativas (não vagas)

**Tempo Estimado:** 4-6 horas

---

## ✍️ FASE 4: REDAÇÃO

### Objetivo

Redigir as seções principais do artigo com rigor técnico e matemático.

### Outputs Obrigatórios

1. **resumo_abstract.md** (250-300 palavras, estrutura IMRAD)
2. **introducao_completa.md** (1.500-2.500 palavras, estrutura CARS)
3. **revisao_literatura_completa.md** (2.000-3.000 palavras)
4. **metodologia_completa.md** (2.500-4.000 palavras, inclui Algorithm 1)
5. **resultados_completo.md** (2.000-3.000 palavras, 9+ tabelas)
6. **discussao_completa.md** (2.500-4.000 palavras)
7. **conclusao_completa.md** (500-800 palavras)
8. **agradecimentos_referencias.md** (refs em ABNT ou APA)

### Template: Algorithm 1 (LaTeX para metodologia)

```latex
\begin{algorithm}[H]
\caption{Experimental Pipeline for Beneficial Noise Analysis}\label{alg:pipeline}
\begin{algorithmic}[1]
\REQUIRE $\mathcal{D}_{train}, \mathcal{D}_{val}, \mathcal{D}_{test}$, 
         Configurations $\mathcal{C}$, Seeds $\mathcal{S}$
\STATE Initialize results table $R \leftarrow \emptyset$
\FOR{each configuration $c \in \mathcal{C}$}
  \STATE Extract $(ansatz, noise\_type, \gamma, schedule, init) \leftarrow c$
  \FOR{each seed $s \in \mathcal{S}$}
    \STATE Set random seed to $s$
    \STATE Initialize quantum model $M_c$ with configuration $c$
    \STATE Apply encoding: $|\psi_{in}\rangle \leftarrow \text{Encode}(x, M_c)$
    \STATE Apply ansatz: $|\psi_{var}\rangle \leftarrow U(\theta) |\psi_{in}\rangle$
    \STATE Apply noise: $\rho_{noisy} \leftarrow \mathcal{N}_\gamma(|\psi_{var}\rangle)$
    \STATE Measure: $\langle Z \rangle \leftarrow \text{Measure}(\rho_{noisy})$
    \STATE Optimize: $\theta^* \leftarrow \text{Adam}(\mathcal{L}, \theta_0, lr)$
    \STATE Evaluate: $m \leftarrow \text{Test}(M_c(\theta^*), \mathcal{D}_{test})$
    \STATE Append $(c, s, m)$ to $R$
  \ENDFOR
\ENDFOR
\STATE Aggregate: $R_{avg} \leftarrow \text{Mean}(R, \text{by}=c)$
\RETURN $R, R_{avg}$
\end{algorithmic}
\end{algorithm}
```

### Template: Tabela Código→Método

```markdown
| Componente do Método | Arquivo/Função/Linha | Parâmetros | Artefatos Gerados |
|----------------------|----------------------|------------|-------------------|
| Definição do Ansatz | `circuits.py:L45-78`, `StronglyEntanglingCircuit` | `n_qubits=4, depth=2` | Objeto PennyLane QNode |
| Canal de Ruído Phase Damping | `noise.py:L15-30`, `phase_damping_channel` | `gamma=0.001431` | Operadores de Kraus K₀, K₁ |
| Otimizador Adam | `train.py:L120`, `qml.AdamOptimizer` | `lr=0.01, beta1=0.9, beta2=0.999` | - |
| Avaliação no Teste | `evaluate.py:L55-70`, `evaluate_model` | `X_test, y_test` | `metrics.json` com acc, precision, recall |
| Geração de Figuras | `plot_utils.py:L100-150`, `plot_accuracy_vs_gamma` | `results_df, figsize=(10,6)` | `figura2b_beneficial_noise.png` |
```

### Estrutura da Introdução (CARS - Create A Research Space)

**Move 1: Estabelecer Território**
- Contextualização: Era NISQ, problema do ruído
- Importância: Computação quântica promete vantagens
- Pesquisa anterior: Trabalhos sobre mitigação de ruído

**Move 2: Estabelecer Nicho (Gap)**
- Lacuna: Ruído sempre visto como obstáculo
- Pergunta: E se ruído pudesse ser benéfico?
- Justificativa: Analogia com dropout em ML clássico

**Move 3: Ocupar o Nicho**
- Nossa abordagem: Investigação sistemática multiparamétrica
- Hipóteses: H₀-H₄
- Contribuições: Primeira validação multiframework, dynamic schedules

### Quality Gate F4

- [ ] Sem números sem lastro (verificar rastreabilidade)
- [ ] R0 respeitado (se aplicável)
- [ ] Methods contém: notação, equações, Algorithm 1, mapa código→método
- [ ] Results: tabelas numeradas, legendas descritivas, valores com IC
- [ ] Discussion: interpretação, não repetição de resultados
- [ ] Conclusion: responde objetivos, limitações honestas, trabalhos futuros específicos

**Tempo Estimado:** 20-30 horas (maior esforço do projeto)

---

## 📊 FASE 5: MATERIAL SUPLEMENTAR

### Objetivo

Fornecer informações adicionais para reprodutibilidade e profundidade.

### Outputs Obrigatórios

1. **tabelas_suplementares.md** (5 tabelas detalhadas)
2. **figuras_suplementares.md** (8 figuras especificadas)
3. **notas_metodologicas_adicionais.md** (detalhes técnicos)

### Tabelas Suplementares Obrigatórias

#### Tabela S1: Configurações Completas (CSV)
```csv
config_id,dataset,ansatz,noise_type,gamma,schedule,init,seed,acc_train,acc_test,epochs,time_sec
1,iris,basic,depolarizing,0.001,static,xavier,42,0.95,0.87,45,12.3
2,iris,basic,depolarizing,0.001,static,xavier,43,0.96,0.88,43,12.1
...
2688,moons,random,correlated,0.01,cosine,kaiming,43,0.62,0.58,50,18.7
```

**Requisitos:**
- Total de linhas = Total configs × seeds (ex: 2.688 × 2 = 5.376)
- Todas as colunas documentadas em `notas_metodologicas_adicionais.md`
- CSV disponível como arquivo separado

#### Tabela S2: Comparação com Estado da Arte
```markdown
| Trabalho | Ano | Datasets | Ansätze | Ruídos | Melhor Acc | Nossa Melhoria |
|----------|-----|----------|---------|--------|------------|----------------|
| Du et al. | 2021 | 1 (Iris) | 1 (Basic) | 1 (Depol) | 89.0% | +3.8% (92.8%) |
| Wang et al. | 2022 | 1 (MNIST) | 1 (Conv) | 0 (Sem ruído) | 94.5% | N/A (domínio diferente) |
| Este trabalho | 2024 | 4 | 7 | 6 | 65.8% (Moons) | - (novo dataset) |
```

#### Tabela S3: Hiperparâmetros Explorados
```markdown
| Hiperparâmetro | Range | # Valores | Método de Busca |
|----------------|-------|-----------|-----------------|
| Learning rate | [0.001, 0.1] | 5 | Log-uniform |
| Batch size | {16, 32, 64} | 3 | Grid |
| Gamma (ruído) | [0, 0.01] | 20 | Linear |
| Depth (camadas) | {1, 2, 3, 4} | 4 | Grid |
| N_qubits | {4, 6, 8} | 3 | Grid |
```

#### Tabela S4: Testes Post-Hoc (Bonferroni)
```markdown
| Comparação | Mean Diff | p-value | p-adj (Bonf) | Significativo? |
|------------|-----------|---------|--------------|----------------|
| Phase Dam vs Depol | +3.75% | 0.0001 | 0.0015 | ✅ Sim (α=0.0033) |
| Phase Dam vs Ampl Dam | +2.13% | 0.0156 | 0.2340 | ❌ Não |
| Depol vs Ampl Dam | -1.62% | 0.0892 | 1.0000 | ❌ Não |
```

#### Tabela S5: Análise de Sensibilidade
```markdown
| Parâmetro | Δ Parâmetro | Δ Acurácia | Sensibilidade (Δacc/Δparam) |
|-----------|-------------|------------|------------------------------|
| Gamma | +10% | +0.5% | 0.05 (baixa) |
| Learning rate | +10% | -2.3% | -0.23 (alta) |
| Depth | +1 camada | +1.8% | 1.80 (alta) |
```

### Figuras Suplementares (8 figuras)

```markdown
# Figuras Suplementares

## Figura S1: Heatmap de Acurácia (Ansatz × Ruído)
**Descrição:** Heatmap 7×6 mostrando acurácia média para cada combinação ansatz-ruído.  
**Eixos:** X=Ruído (6 tipos), Y=Ansatz (7 tipos), Color=Acurácia [0,1]  
**Escala:** Viridis colormap, 300 DPI, formato PNG  
**Achado-chave:** Strongly Entangling + Phase Damping = melhor combo (canto superior direito)  
**Script:** `plot_utils.py:L200-250`

## Figura S2: Curvas de Aprendizado (Train vs Test)
**Descrição:** 6 subplots (um por tipo de ruído) com curvas acc_train e acc_test vs época.  
**Eixos:** X=Época [0,50], Y=Acurácia [0,1]  
**Linhas:** Azul=Treino, Laranja=Teste, área sombreada=IC 95%  
**Achado-chave:** Ruído reduz gap treino-teste (menos overfitting)  
**Script:** `plot_utils.py:L300-380`

[... 6 figuras adicionais ...]
```

### Quality Gate F5

- [ ] Tabela S1 confere com total calculado na Fase 1
- [ ] Cada tabela/figura aponta script/config/log ou [INFORMAÇÃO AUSENTE]
- [ ] Nada "core" da metodologia ficou apenas no suplemento
- [ ] Todas as figuras têm especificação completa (eixos, escalas, colormap, DPI)
- [ ] CSV da Tabela S1 disponível e validado (sem linhas faltando)

**Tempo Estimado:** 8-12 horas

---

## ✅ FASE 6: CONSOLIDAÇÃO

### Objetivo

Garantir consistência, rastreabilidade e conformidade final.

### Outputs Obrigatórios

1. **relatorio_conivencia.md** (% de conivência + discrepâncias)
2. **rastreabilidade_completa.md** (tabela preenchida)
3. **tabela_codigo_metodo.md** (mapeamento completo)
4. **artigo_completo_final.md** ou **.tex** (consolidado)
5. **sumario_executivo_final.md** (visão geral)

### Template: Tabela de Rastreabilidade

```markdown
| Seção | Afirmação/Número | Evidência (Arquivo:Linha) | Referência (Autor, Ano) |
|-------|------------------|---------------------------|-------------------------|
| 4.5 Results | Acurácia máxima 65.83% | `resultados_multiframework_20251226_172214/resultados_completos.json:row_1523:col_accuracy` | - |
| 4.4 Methods | Equação de Lindblad (Eq. 3) | `noise.py:L15-20` | (Lindblad, 1976) |
| 4.6 Discussion | Ruído como regularizador | `fase2_bibliografia/sintese_literatura.md:L145-167` | (Du et al., 2021; Srivastava et al., 2014) |
| 4.5 Results | Cohen's d = 4.03 | `analysis/statistics.py:L89` + `results_stats.json:effect_size` | - |
| 3.2 Methods | Seeds [42, 43] | `config_artigo.json:L18` + `framework_investigativo_completo.py:L12` | - |
```

**Requisitos:**
- Mínimo 20 entradas rastreáveis
- 100% das afirmações quantitativas principais rastreadas
- Arquivo:Linha específico (não vago como "código")
- Referências quando aplicável

### Checklist de Consistência

```bash
# Verificar números consistentes
python tools/check_consistency.py --code framework_investigativo_completo.py \
  --text artigo_cientifico/fase4_secoes/resultados_completo.md \
  --output relatorio_conivencia.md

# Exemplo output esperado:
# ✅ Acurácia 65.83%: Código=65.827%, Texto=65.83% (Δ=0.003%, OK)
# ✅ Seeds [42,43]: Ambos documentos consistentes
# ❌ Total configs: Código=2688, Texto=2680 (DISCREPÂNCIA!)
```

### Quality Gate Final

- [ ] Consistência ≥ 95% (meta: 100% após ajustes)
- [ ] Citação↔referência 100% (nenhuma citação órfã)
- [ ] Reprodutibilidade: ambiente + scripts + seeds documentados
- [ ] Ameaças à validade + scope conditions explicitados
- [ ] Checklist de auditoria ≥ 90/100 pontos
- [ ] Código público (GitHub) com README e LICENSE
- [ ] Dados públicos ou justificativa de privacidade

**Tempo Estimado:** 6-8 horas

---

## 🎓 CONSOLIDAÇÃO FINAL E SUBMISSÃO

### Checklist Pré-Submissão (ABNT/LaTeX)

#### Para MODE_A (Inglês/LaTeX)

- [ ] Compilar LaTeX sem erros/warnings
  ```bash
  cd artigo_cientifico/latex_template
  pdflatex npj_qi_submission.tex
  bibtex npj_qi_submission
  pdflatex npj_qi_submission.tex
  pdflatex npj_qi_submission.tex
  ```
- [ ] Verificar formatação:
  - [ ] Margens: 2.5cm (Nature/npj QI)
  - [ ] Fonte: Times New Roman 11pt ou Computer Modern
  - [ ] Espaçamento: 1.5 linhas
  - [ ] Numeração: Páginas, seções, equações, tabelas, figuras
- [ ] Figuras:
  - [ ] Resolução ≥ 300 DPI
  - [ ] Formato: PDF (vetorial) ou PNG/TIFF (raster)
  - [ ] Legendas completas e auto-explicativas
  - [ ] Citadas no texto antes de aparecerem

#### Para MODE_B (Português/ABNT)

- [ ] Verificar ABNT NBR 6023 (Referências):
  - [ ] Autor em CAIXA ALTA
  - [ ] Título em **negrito**
  - [ ] Periódico em *itálico*
  - [ ] DOI presente quando disponível
- [ ] Verificar ABNT NBR 10520 (Citações):
  - [ ] (AUTOR, ano) para citações indiretas
  - [ ] "Citação direta entre aspas" (AUTOR, ano, p. X)
- [ ] Elementos pré-textuais:
  - [ ] Capa
  - [ ] Folha de rosto
  - [ ] Resumo (português)
  - [ ] Abstract (inglês)
  - [ ] Lista de figuras
  - [ ] Lista de tabelas
  - [ ] Sumário

### Submissão ao Periódico

#### Documentos Necessários

1. **Manuscrito principal** (PDF ou LaTeX source)
2. **Material suplementar** (ZIP com CSVs, figuras suplementares, código)
3. **Cover letter** (1 página)
4. **Author contributions statement**
5. **Data availability statement**
6. **Conflict of interest statement**
7. **Funding information**

#### Template Cover Letter

```markdown
Dear Editor,

We submit for your consideration the manuscript entitled "From Obstacle to 
Opportunity: Harnessing Beneficial Quantum Noise in Variational Classifiers" 
for publication in [Journal Name].

This work presents the first systematic investigation of beneficial quantum 
noise across multiple frameworks (PennyLane, Qiskit, Cirq), demonstrating 
that controlled noise can improve generalization in variational quantum 
classifiers. Our key contributions include:

1. First multi-framework validation (3 platforms, identical configs)
2. Novel dynamic noise schedules (Cosine, Linear) - original to this work
3. Comprehensive statistical analysis (ANOVA, effect sizes, post-hoc tests)
4. Complete reproducibility (code, data, seeds publicly available)

We believe this work is of high interest to [Journal Name] readers because 
[justificativa específica do periódico].

All authors have approved the manuscript and agree with its submission. 
The manuscript has not been published elsewhere and is not under consideration 
by another journal.

Sincerely,
[Nome do Corresponding Author]
[Afiliação]
[Email]
```

---

## 🔧 FERRAMENTAS E AUTOMAÇÃO

### Scripts Principais

```bash
# 1. Gerador completo (todas as 6 fases)
python gerador_artigo_completo.py --config config_artigo.json

# 2. Verificar consistência código-texto
python tools/check_consistency.py --code framework_*.py --text artigo_cientifico/fase4_secoes/

# 3. Gerar Tabela S1 automaticamente
python tools/generate_s1.py --results resultados_multiframework_*/resultados_completos.json --output fase5_suplementar/tabela_s1_configuracoes.csv

# 4. Consolidar artigo final
python tools/build_paper.py --mode latex --output artigo_completo_final.tex

# 5. Verificar auditoria
python tools/audit_checker.py --checklist CHECKLIST_AUDITORIA_COMPLETO.md --artigo artigo_cientifico/
```

### Ambientes Recomendados

**Para escrita (MODE_A - LaTeX):**
- Overleaf (online, colaborativo)
- TeXstudio (desktop, avançado)
- VSCode + LaTeX Workshop extension

**Para escrita (MODE_B - Markdown/ABNT):**
- Typora (WYSIWYG Markdown)
- VSCode + Markdown All in One
- Obsidian (knowledge graph)

**Para gerenciamento de referências:**
- Zotero (grátis, open-source)
- Mendeley (grátis, integrado)
- EndNote (pago, robusto)

**Para versionamento:**
- Git + GitHub (essencial)
- GitLab (alternativa)
- Zenodo (DOI para código/dados)

---

## 📊 CRONOGRAMA COMPLETO

Ver: `CRONOGRAMA_ESTIMADO_COMPLETO.md` para detalhes

**Resumo:**
- **Fase 1 (Auditoria):** 8-12h → Semana 1
- **Fase 2 (Bibliografia):** 6-10h (R0) ou 15-25h (R1) → Semana 1-2
- **Fase 3 (Projeto):** 4-6h → Semana 2
- **Fase 4 (Redação):** 20-30h → Semanas 3-5
- **Fase 5 (Suplementar):** 8-12h → Semana 6
- **Fase 6 (Consolidação):** 6-8h → Semana 6
- **Total:** 52-78h (6-10 dias úteis, 1-2 meses calendário)

---

## 🆘 TROUBLESHOOTING

Ver: `FAQ_TROUBLESHOOTING_COMPLETO.md` para 30+ perguntas e respostas

**Top 5 Problemas Frequentes:**

1. **"Meu código não tem seeds fixas"**
   → Adicionar seeds, ou executar N vezes e reportar média ± DP

2. **"Não tenho logs de execução"**
   → Executar pipeline e gerar, ou marcar resultados como [NÃO DISPONÍVEL]

3. **"Faltam referências para afirmações"**
   → Se R0: marcar [LACUNA DE CITAÇÃO]; Se R1: buscar em 7 categorias

4. **"Execução muito lenta"**
   → Usar Bayesian Optimization, paralelizar, ou reduzir configurações

5. **"Números inconsistentes entre código e texto"**
   → Usar `tools/check_consistency.py` para detectar e corrigir

---

## 📚 REFERÊNCIAS PRINCIPAIS DO GUIA

1. **Metodologia Científica:**
   - Creswell, J. W. (2014). *Research Design*. 4th ed. SAGE.
   - Yin, R. K. (2017). *Case Study Research*. 6th ed. SAGE.

2. **Escrita Acadêmica:**
   - Swales, J. M. (1990). *Genre Analysis*. Cambridge University Press.
   - Belcher, W. L. (2019). *Writing Your Journal Article in Twelve Weeks*. 2nd ed.

3. **Reprodutibilidade:**
   - Wilkinson et al. (2016). The FAIR Guiding Principles. *Scientific Data*.
   - Peng, R. D. (2011). Reproducible research. *Science*, 334(6060), 1226-1227.

4. **Estatística:**
   - Cohen, J. (1988). *Statistical Power Analysis*. 2nd ed.
   - Field, A. (2013). *Discovering Statistics Using IBM SPSS*. 4th ed.

5. **Quantum ML:**
   - Cerezo et al. (2021). Variational Quantum Algorithms. *Nature Reviews Physics*.
   - Preskill, J. (2018). Quantum Computing in the NISQ era. *Quantum*, 2, 79.

---

## ✅ VALIDAÇÃO FINAL

**Antes de submeter, verificar:**

- [ ] Pontuação no Checklist de Auditoria ≥ 90/100
- [ ] Conivência código-texto ≥ 95%
- [ ] Formato ABNT/LaTeX 100% correto
- [ ] Todas as hipóteses testadas e respondidas
- [ ] Todos os objetivos atingidos
- [ ] Limitações honestamente discutidas
- [ ] Trabalhos futuros específicos (não genéricos como "testar em outros datasets")
- [ ] Código e dados disponíveis publicamente (ou justificativa)
- [ ] Material Suplementar completo (5 tabelas + 8 figuras + notas)
- [ ] Sumário executivo criado
- [ ] Todos os arquivos de auditoria entregues
- [ ] Manuscrito compilado e revisado por ≥2 co-autores

---

**PARABÉNS! Você completou o framework de geração de artigos QUALIS A1 com rastreabilidade total.**

**Próximos passos:**
1. Submeter ao periódico escolhido
2. Aguardar revisão por pares
3. Responder a reviews (usar rastreabilidade para respostas precisas)
4. Celebrar aceitação! 🎉

---

**Versão:** 2.0  
**Última Atualização:** 26/12/2025  
**Mantenedor:** Framework Geração Artigos QUALIS A1  
**Licença:** MIT  
**Status:** ✅ Completo e testado

**Para dúvidas:**  
- Consulte `FAQ_TROUBLESHOOTING_COMPLETO.md`  
- Abra issue no GitHub  
- Email: [inserir contato]
