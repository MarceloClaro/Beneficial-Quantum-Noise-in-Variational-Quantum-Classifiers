# Checklist de Auditoria Completo (0-100 pontos)

**Sistema de Pontuação QUALIS A1**  
**Versão:** 1.0  
**Data:** 26/12/2025  
**Conformidade:** Nature, Science, Physical Review, Quantum, npj QI

---

## 📊 VISÃO GERAL

Este checklist fornece um sistema objetivo de avaliação para artigos científicos, garantindo conformidade com padrões QUALIS A1. A pontuação máxima é 100 pontos, distribuídos em 4 categorias principais.

**Interpretação da Pontuação:**
- **90-100:** 🥇 Excelente - Pronto para submissão a periódicos top
- **80-89:** 🥈 Muito Bom - Necessita pequenos ajustes
- **70-79:** 🥉 Bom - Necessita ajustes moderados
- **60-69:** ⚠️ Satisfatório - Necessita melhorias significativas
- **<60:** ❌ Insatisfatório - Requer revisão profunda

---

## 🎯 CATEGORIAS E PONTUAÇÃO

### CATEGORIA 1: REPRODUTIBILIDADE (30 pontos)

#### 1.1 Ambiente Documentado (10 pontos)

**Critérios:**
- [ ] **4 pts:** Sistema operacional especificado (nome, versão, arquitetura)
- [ ] **3 pts:** Hardware documentado (CPU, RAM, GPU se usado, storage)
- [ ] **3 pts:** Todas as bibliotecas com versões exatas (ex: PennyLane==0.38.0)

**Como verificar:**
```bash
# Verificar se existe ENVIRONMENT.md ou seção de Environment em Methods
grep -r "Operating System\|CPU\|RAM\|Python.*[0-9]" artigo_cientifico/
```

**Exemplos:**

✅ **10/10 - Excelente:**
```markdown
## Computational Environment
- **OS:** Ubuntu 22.04 LTS (kernel 5.15.0-91)
- **CPU:** Intel Xeon E5-2680 v4 @ 2.40GHz (28 cores)
- **RAM:** 128 GB DDR4 ECC
- **Python:** 3.9.18
- **PennyLane:** 0.38.0
- **NumPy:** 1.24.3
- **SciPy:** 1.10.1
```

❌ **3/10 - Insatisfatório:**
```markdown
Python 3.x, libraries as needed
```

---

#### 1.2 Seeds Fixas e Reportadas (10 pontos)

**Critérios:**
- [ ] **5 pts:** Seeds explicitamente fixas no código (ex: `np.random.seed(42)`)
- [ ] **3 pts:** Valores das seeds reportados no artigo
- [ ] **2 pts:** Justificativa da escolha das seeds ou múltiplas seeds testadas

**Como verificar:**
```bash
# Verificar código
grep -n "seed\|random_state" framework_*.py

# Verificar artigo
grep -n "seed\|random state" artigo_cientifico/fase4_secoes/metodologia_completa.md
```

**Exemplos:**

✅ **10/10 - Excelente:**
```python
# Código
SEEDS = [42, 43, 44]  # Multiple seeds for robustness
for seed in SEEDS:
    np.random.seed(seed)
    torch.manual_seed(seed)
    result = run_experiment()
```
```markdown
<!-- Artigo -->
Para garantir reprodutibilidade, fixamos três seeds aleatórias (42, 43, 44) 
e reportamos média ± desvio padrão entre as execuções.
```

❌ **2/10 - Insatisfatório:**
```python
# Sem seed fixa
result = run_experiment()
```

---

#### 1.3 Pipeline Executável (10 pontos)

**Critérios:**
- [ ] **4 pts:** Script principal executa sem erros
- [ ] **3 pts:** Comandos de execução documentados (README ou Methods)
- [ ] **3 pts:** Tempo de execução estimado reportado

**Como verificar:**
```bash
# Testar execução
python framework_investigativo_completo.py --help
python framework_investigativo_completo.py --quick-test

# Verificar documentação
grep -n "python.*framework\|bash\|sh" README.md
```

**Exemplos:**

✅ **10/10 - Excelente:**
```bash
# Modo rápido (teste - 10 minutos)
python framework_investigativo_completo.py --quick-test --seed 42

# Modo completo (48-72 horas)
python framework_investigativo_completo.py --full --seeds 42 43 44
```

❌ **3/10 - Insatisfatório:**
```
(Sem instruções de execução)
```

---

### CATEGORIA 2: RASTREABILIDADE (30 pontos)

#### 2.1 Tabela de Rastreabilidade Completa (15 pontos)

**Critérios:**
- [ ] **5 pts:** Tabela existe e está preenchida
- [ ] **5 pts:** Todas as afirmações quantitativas têm origem rastreável
- [ ] **3 pts:** Formato consistente (Seção | Afirmação | Evidência | Referência)
- [ ] **2 pts:** Evidências específicas (arquivo:linha, não apenas "código")

**Como verificar:**
```bash
# Verificar existência
ls artigo_cientifico/fase6_consolidacao/rastreabilidade_completa.md

# Contar entradas
grep -c "^\|" artigo_cientifico/fase6_consolidacao/rastreabilidade_completa.md
```

**Template obrigatório:**
```markdown
| Seção | Afirmação/Número | Evidência (Arquivo:Linha) | Referência |
|-------|------------------|---------------------------|------------|
| 4.5 Results | Acurácia 65.83% | `resultados.csv:row_1523:col_acc` | - |
| 4.4 Methods | Eq. de Lindblad | `noise.py:L15-20` | (Lindblad, 1976) |
```

**Pontuação:**
- ≥20 entradas rastreáveis: 15 pts
- 10-19 entradas: 10 pts
- 5-9 entradas: 5 pts
- <5 entradas: 0 pts

---

#### 2.2 Mapa Código→Método Completo (15 pontos)

**Critérios:**
- [ ] **5 pts:** Tabela Código→Método existe
- [ ] **5 pts:** Todos os componentes metodológicos mapeados
- [ ] **3 pts:** Inclui arquivo, função/classe, linhas, parâmetros
- [ ] **2 pts:** Inclui artefatos gerados (figuras, tabelas, logs)

**Como verificar:**
```bash
ls artigo_cientifico/fase6_consolidacao/tabela_codigo_metodo.md
```

**Template obrigatório:**
```markdown
| Componente do Método | Arquivo/Função/Linha | Parâmetros | Artefatos Gerados |
|----------------------|----------------------|------------|-------------------|
| Ansatz Strongly Entangling | `circuits.py:L45-78` | `n_qubits=4, depth=2` | Objeto PennyLane |
| Canal Phase Damping | `noise.py:L15-30` | `gamma=0.001431` | Operador Kraus |
| Otimizador Adam | `train.py:L120` | `lr=0.01, beta1=0.9` | - |
| Avaliação Teste | `evaluate.py:L55-70` | `X_test, y_test` | `metrics.json` |
```

**Pontuação:**
- ≥15 componentes mapeados: 15 pts
- 10-14 componentes: 10 pts
- 5-9 componentes: 5 pts
- <5 componentes: 0 pts

---

### CATEGORIA 3: RIGOR ESTATÍSTICO (20 pontos)

#### 3.1 Testes Apropriados (5 pontos)

**Critérios:**
- [ ] **2 pts:** Teste de hipótese apropriado (t-test, ANOVA, Mann-Whitney, etc.)
- [ ] **2 pts:** Verificação de pressupostos (normalidade, homoscedasticidade)
- [ ] **1 pt:** Reportar estatística de teste e graus de liberdade

**Como verificar:**
```bash
grep -i "t-test\|ANOVA\|Kruskal\|Mann-Whitney\|Friedman" artigo_cientifico/fase4_secoes/resultados_completo.md
```

**Exemplos:**

✅ **5/5 - Excelente:**
```markdown
Aplicamos ANOVA de três fatores (ruído × ansatz × dataset) após verificar 
normalidade (Shapiro-Wilk, p=0.18) e homoscedasticidade (Levene, p=0.24). 
Encontramos efeito principal significativo do ruído (F(5,2682)=28.3, p<0.001).
```

❌ **1/5 - Insatisfatório:**
```markdown
Os resultados foram analisados estatisticamente.
```

---

#### 3.2 Correção para Múltiplas Comparações (5 pontos)

**Critérios:**
- [ ] **3 pts:** Correção aplicada (Bonferroni, Benjamini-Hochberg, Holm, etc.)
- [ ] **2 pts:** Justificativa da escolha do método de correção

**Como verificar:**
```bash
grep -i "Bonferroni\|Benjamini\|Hochberg\|FDR\|Holm\|multiple.*comparison" artigo_cientifico/fase4_secoes/resultados_completo.md
```

**Exemplos:**

✅ **5/5 - Excelente:**
```markdown
Para múltiplas comparações post-hoc (15 pares), aplicamos correção de 
Bonferroni (α=0.05/15=0.0033). Phase Damping superou Depolarizing 
(p_adj=0.0012 < 0.0033, significativo após correção).
```

❌ **0/5 - Insatisfatório:**
```markdown
Phase Damping foi melhor (p=0.03).
<!-- 15 comparações, sem correção → falso positivo! -->
```

---

#### 3.3 Intervalos de Confiança (5 pontos)

**Critérios:**
- [ ] **3 pts:** IC reportados para médias principais (95% ou 99%)
- [ ] **2 pts:** IC incluem margem de erro ou [lower, upper]

**Como verificar:**
```bash
grep -E "IC 95%|\[.*,.*\]|±.*%" artigo_cientifico/fase4_secoes/resultados_completo.md
```

**Exemplos:**

✅ **5/5 - Excelente:**
```markdown
Acurácia ótima: 65.83% [IC 95%: 63.69%, 67.97%] ou 65.83% ± 2.14%
```

❌ **0/5 - Insatisfatório:**
```markdown
Acurácia ótima: 65.83%
<!-- Sem indicação de incerteza -->
```

---

#### 3.4 Tamanhos de Efeito (5 pontos)

**Critérios:**
- [ ] **3 pts:** Effect sizes calculados (Cohen's d, η², Cliff's Δ)
- [ ] **1 pt:** Interpretação qualitativa (pequeno/médio/grande)
- [ ] **1 pt:** IC para effect sizes (se possível)

**Como verificar:**
```bash
grep -i "Cohen.*d\|eta.*squared\|Cliff.*delta\|effect.*size" artigo_cientifico/fase4_secoes/resultados_completo.md
```

**Exemplos:**

✅ **5/5 - Excelente:**
```markdown
O tamanho de efeito foi muito grande (Cohen's d=4.03 [IC 95%: 3.21, 4.85]), 
indicando diferença substancial entre configurações com e sem ruído.
```

❌ **0/5 - Insatisfatório:**
```markdown
A diferença foi estatisticamente significativa (p<0.001).
<!-- Sem magnitude do efeito -->
```

---

### CATEGORIA 4: TRANSPARÊNCIA (20 pontos)

#### 4.1 Código Disponível Publicamente (10 pontos)

**Critérios:**
- [ ] **5 pts:** Repositório GitHub/GitLab público com código completo
- [ ] **2 pts:** README com instruções de instalação e uso
- [ ] **2 pts:** Licença especificada (MIT, Apache, GPL, etc.)
- [ ] **1 pt:** DOI via Zenodo/Figshare para versionamento

**Como verificar:**
```bash
# Verificar .git
ls -la .git

# Verificar LICENSE
cat LICENSE

# Verificar README
grep -i "installation\|usage\|how to" README.md
```

**Pontuação:**
- Público + README + LICENSE + DOI: 10 pts
- Público + README + LICENSE: 8 pts
- Público + README: 5 pts
- Privado ou sem instruções: 0 pts

---

#### 4.2 Dados Disponíveis Publicamente (5 pontos)

**Critérios:**
- [ ] **3 pts:** Dados brutos ou processados disponíveis (GitHub, Zenodo, OSF)
- [ ] **1 pt:** Metadados descritivos (formato, tamanho, descrição das colunas)
- [ ] **1 pt:** Licença de dados (CC BY, CC0, Open Database License)

**Como verificar:**
```bash
ls resultados_multiframework_*/resultados_completos.json
ls resultados_multiframework_*/execution_manifest.json
```

**Pontuação:**
- Dados + metadados + licença: 5 pts
- Dados + metadados: 4 pts
- Dados apenas: 2 pts
- Dados privados/indisponíveis: 0 pts

---

#### 4.3 Limitações e Ameaças à Validade Discutidas (5 pontos)

**Critérios:**
- [ ] **2 pts:** Seção explícita "Limitations" ou "Threats to Validity"
- [ ] **1 pt:** Ameaças internas identificadas (confounders, vieses)
- [ ] **1 pt:** Ameaças externas identificadas (generalização, escala)
- [ ] **1 pt:** Ameaças de construto e estatísticas identificadas

**Como verificar:**
```bash
grep -i "limitation\|threat.*validity\|weakness\|caveat" artigo_cientifico/fase4_secoes/discussao_completa.md
```

**Exemplos:**

✅ **5/5 - Excelente:**
```markdown
## 5.7 Limitations and Threats to Validity

**Internal Validity (Causalidade):**
- Confounding: Learning rate e ruído co-variaram em algumas configurações

**External Validity (Generalização):**
- Escala: Experimentos limitados a 4 qubits (hardware NISQ atual)
- Datasets: Toy problems (Moons, Iris); aplicabilidade a problemas reais não verificada

**Construct Validity (Medição):**
- Acurácia como métrica única; outras métricas (F1, AUC) não exploradas

**Statistical Validity (Inferência):**
- Poder estatístico: Com N=2.688, poder > 0.99 para detectar d≥0.5
```

❌ **1/5 - Insatisfatório:**
```markdown
Nossos resultados têm algumas limitações que serão abordadas em trabalhos futuros.
<!-- Vago, sem especificidade -->
```

---

## 📋 CHECKLIST CONSOLIDADO

### Instruções de Uso

1. **Preencher cada item** marcando [x] quando satisfeito
2. **Calcular pontuação** ao final de cada categoria
3. **Somar pontos** para obter score total (0-100)
4. **Interpretar resultado** segundo escala (Excelente/Muito Bom/Bom/Satisfatório/Insatisfatório)
5. **Priorizar melhorias** nos itens de menor pontuação

---

### CATEGORIA 1: REPRODUTIBILIDADE (30 pontos)

#### 1.1 Ambiente Documentado (10 pontos)
- [ ] **4 pts:** Sistema operacional especificado
- [ ] **3 pts:** Hardware documentado
- [ ] **3 pts:** Bibliotecas com versões exatas

**Subtotal 1.1:** _____ / 10

#### 1.2 Seeds Fixas (10 pontos)
- [ ] **5 pts:** Seeds fixas no código
- [ ] **3 pts:** Seeds reportadas no artigo
- [ ] **2 pts:** Justificativa ou múltiplas seeds

**Subtotal 1.2:** _____ / 10

#### 1.3 Pipeline Executável (10 pontos)
- [ ] **4 pts:** Script executa sem erros
- [ ] **3 pts:** Comandos documentados
- [ ] **3 pts:** Tempo estimado reportado

**Subtotal 1.3:** _____ / 10

**TOTAL CATEGORIA 1:** _____ / 30

---

### CATEGORIA 2: RASTREABILIDADE (30 pontos)

#### 2.1 Tabela de Rastreabilidade (15 pontos)
- [ ] **5 pts:** Tabela existe e preenchida
- [ ] **5 pts:** Todas afirmações quantitativas rastreáveis
- [ ] **3 pts:** Formato consistente
- [ ] **2 pts:** Evidências específicas (arquivo:linha)

**Subtotal 2.1:** _____ / 15

#### 2.2 Mapa Código→Método (15 pontos)
- [ ] **5 pts:** Tabela existe
- [ ] **5 pts:** Todos componentes mapeados
- [ ] **3 pts:** Inclui arquivo/função/linhas/parâmetros
- [ ] **2 pts:** Inclui artefatos gerados

**Subtotal 2.2:** _____ / 15

**TOTAL CATEGORIA 2:** _____ / 30

---

### CATEGORIA 3: RIGOR ESTATÍSTICO (20 pontos)

#### 3.1 Testes Apropriados (5 pontos)
- [ ] **2 pts:** Teste de hipótese apropriado
- [ ] **2 pts:** Verificação de pressupostos
- [ ] **1 pt:** Reportar estatística e graus de liberdade

**Subtotal 3.1:** _____ / 5

#### 3.2 Correção Múltiplas Comparações (5 pontos)
- [ ] **3 pts:** Correção aplicada
- [ ] **2 pts:** Justificativa da escolha

**Subtotal 3.2:** _____ / 5

#### 3.3 Intervalos de Confiança (5 pontos)
- [ ] **3 pts:** IC reportados (95% ou 99%)
- [ ] **2 pts:** IC com margem de erro ou [lower, upper]

**Subtotal 3.3:** _____ / 5

#### 3.4 Tamanhos de Efeito (5 pontos)
- [ ] **3 pts:** Effect sizes calculados
- [ ] **1 pt:** Interpretação qualitativa
- [ ] **1 pt:** IC para effect sizes

**Subtotal 3.4:** _____ / 5

**TOTAL CATEGORIA 3:** _____ / 20

---

### CATEGORIA 4: TRANSPARÊNCIA (20 pontos)

#### 4.1 Código Público (10 pontos)
- [ ] **5 pts:** Repositório público
- [ ] **2 pts:** README com instruções
- [ ] **2 pts:** Licença especificada
- [ ] **1 pt:** DOI via Zenodo/Figshare

**Subtotal 4.1:** _____ / 10

#### 4.2 Dados Públicos (5 pontos)
- [ ] **3 pts:** Dados disponíveis
- [ ] **1 pt:** Metadados descritivos
- [ ] **1 pt:** Licença de dados

**Subtotal 4.2:** _____ / 5

#### 4.3 Limitações Discutidas (5 pontos)
- [ ] **2 pts:** Seção explícita de Limitações
- [ ] **1 pt:** Ameaças internas identificadas
- [ ] **1 pt:** Ameaças externas identificadas
- [ ] **1 pt:** Ameaças de construto/estatísticas

**Subtotal 4.3:** _____ / 5

**TOTAL CATEGORIA 4:** _____ / 20

---

## 🎯 PONTUAÇÃO FINAL

**SOMA TOTAL:** _____ / 100

**CLASSIFICAÇÃO:**
- [ ] 🥇 90-100: Excelente - Pronto para submissão
- [ ] 🥈 80-89: Muito Bom - Pequenos ajustes
- [ ] 🥉 70-79: Bom - Ajustes moderados
- [ ] ⚠️ 60-69: Satisfatório - Melhorias significativas
- [ ] ❌ <60: Insatisfatório - Revisão profunda

---

## 📊 ANÁLISE POR CATEGORIA

### Gráfico de Radar (preencher manualmente)

```
Reprodutibilidade: _____ / 30 = _____% 
Rastreabilidade:   _____ / 30 = _____%
Rigor Estatístico: _____ / 20 = _____%
Transparência:     _____ / 20 = _____%
```

### Prioridades de Melhoria

1. **Categoria com menor pontuação:** ________________
2. **Item específico com 0 pontos:** ________________
3. **Meta para próxima revisão:** ________________

---

## 🔧 PLANO DE AÇÃO (se <90 pontos)

### Se Reprodutibilidade < 25/30
- [ ] Adicionar `ENVIRONMENT.md` completo
- [ ] Fixar seeds e documentar
- [ ] Testar execução do pipeline

### Se Rastreabilidade < 25/30
- [ ] Criar `rastreabilidade_completa.md`
- [ ] Preencher `tabela_codigo_metodo.md`
- [ ] Mapear todas afirmações quantitativas

### Se Rigor Estatístico < 16/20
- [ ] Adicionar correção para múltiplas comparações
- [ ] Calcular e reportar intervalos de confiança
- [ ] Calcular e interpretar tamanhos de efeito

### Se Transparência < 16/20
- [ ] Tornar repositório público
- [ ] Adicionar seção explícita de Limitações
- [ ] Obter DOI via Zenodo

---

## 📝 NOTAS ADICIONAIS

### Caso Este Projeto (Beneficial Quantum Noise)

**Pontuação Estimada:** 91-95 / 100 (🥇 Excelente)

**Pontos Fortes:**
- ✅ Código público no GitHub
- ✅ Seeds fixas documentadas [42, 43]
- ✅ Ambiente detalhado (3 frameworks)
- ✅ Rastreabilidade completa
- ✅ ANOVA multifatorial + post-hoc
- ✅ Effect sizes (Cohen's d = 4.03)
- ✅ ICs reportados

**Oportunidades de Melhoria (+5 pontos):**
- [ ] Adicionar DOI via Zenodo (+1 pt)
- [ ] Expandir seção de Limitações com ameaças específicas (+2 pts)
- [ ] Licença de dados explícita para CSV (+1 pt)
- [ ] Verificação de pressupostos ANOVA detalhada (+1 pt)

---

## 📚 REFERÊNCIAS DO CHECKLIST

1. **Reprodutibilidade:**
   - Wilkinson et al. (2016). The FAIR Guiding Principles. *Scientific Data*.
   - Peng, R. D. (2011). Reproducible research in computational science. *Science*, 334(6060), 1226-1227.

2. **Rigor Estatístico:**
   - Cohen, J. (1988). *Statistical power analysis for the behavioral sciences*. 2nd ed.
   - Field, A. (2013). *Discovering statistics using IBM SPSS statistics*. 4th ed.

3. **Transparência:**
   - Nosek et al. (2015). Promoting an open research culture. *Science*, 348(6242), 1422-1425.
   - Miguel et al. (2014). Promoting transparency in social science research. *Science*, 343(6166), 30-31.

---

**Versão:** 1.0  
**Última Atualização:** 26/12/2025  
**Compatível com:** QUALIS A1, Nature, Science, PR, Quantum, npj QI  
**Licença:** CC BY 4.0
