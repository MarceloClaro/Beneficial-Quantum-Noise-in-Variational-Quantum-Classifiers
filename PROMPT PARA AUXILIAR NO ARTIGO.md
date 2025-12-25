# Geração de Artigo Científico Qualis A1 de Alto Impacto

## 🎯 OBJETIVO GERAL

Gerar um artigo científico completo e rigoroso, pronto para submissão a periódicos de alto impacto (Nature, Science, Quantum, Physical Review, Qualis A1), com 100% de conivência entre código/dados e texto, garantindo reprodutibilidade, auditabilidade e máxima avaliação por bancas de revisão.

---

## 📋 ESTRUTURA QUE OBRIGATORIAMENTE DEVE SEGUIR DE FORMA MINUCIOSA COM RIGOR QUALIS A1:

### FASE 1: ANÁLISE INICIAL E PLANEJAMENTO

#### 1.1 Análise do Código/Dados Fornecidos

**Instrução:**
```
Analise o código/dados fornecidos e extraia:

1. **Estrutura Técnica:**
   - Número total de linhas de código
   - Módulos principais e suas funções
   - Classes implementadas (nome, propósito, métodos)
   - Arquiteturas/modelos utilizados (listar todos)
   - Técnicas estatísticas/analíticas (listar todas)
   - Bibliotecas e versões (especificar exatas)
   
2. **Componentes Experimentais:**
   - Fatores experimentais (listar todos)
   - Níveis de cada fator (especificar)
   - Total de configurações (calcular: fator1 × fator2 × ... × fatorN)
   - Datasets utilizados (nome, tamanho, características)
   - Métricas de avaliação (listar todas)
   
3. **Metodologia Implementada:**
   - Pré-processamento de dados (passos)
   - Treinamento/otimização (algoritmos, hiperparâmetros)
   - Validação (estratégia, número de folds/repetições)
   - Análise estatística (testes, correções, tamanhos de efeito)
   
4. **Inovações e Contribuições:**
   - Técnicas originais implementadas
   - Diferenças em relação ao estado da arte
   - Contribuições metodológicas específicas

**Formato de Saída:**
- Documento Markdown estruturado: `analise_codigo_inicial.md`
- Tabela resumo de componentes
- Contagem precisa de configurações experimentais
```

#### 1.2 Identificação da Linha de Pesquisa

**Instrução:**
```
Com base na análise do código, identifique:

1. **Área de Pesquisa:** [Ex: Computação Quântica, Machine Learning, Bioinformática]
2. **Subárea Específica:** [Ex: Variational Quantum Algorithms, Deep Learning, Genômica]
3. **Problema Central:** [Descrever em 1-2 frases]
4. **Linha de Pesquisa (Autor Fundacional):** [Ex: Du et al. 2021 - Ruído Benéfico]
5. **Trabalhos Seminais (3-5):** [Listar com citação completa]
6. **Lacuna Identificada:** [Descrever em 3 dimensões: generalidade, dinâmica, interação]

**Formato de Saída:**
- Documento: `linha_de_pesquisa.md`
- Diagrama conceitual (Mermaid/D2) mostrando evolução da linha
```

---

### FASE 2: PESQUISA BIBLIOGRÁFICA PROFUNDA

#### 2.1 Busca de Referências Relevantes

**Instrução:**
```
Realize busca profunda de referências em 7 categorias:

1. **Trabalhos Fundacionais (5-8 refs):**
   - Artigos seminais que estabeleceram a área
   - Critério: >500 citações, publicados em Nature/Science/PRL
   
2. **Estado da Arte Recente (8-12 refs):**
   - Artigos dos últimos 3 anos na mesma linha
   - Critério: Periódicos Qualis A1, relevância direta
   
3. **Metodologia e Técnicas (6-10 refs):**
   - Referências para cada técnica utilizada no código
   - Critério: Artigo original que propôs a técnica
   
4. **Análise Estatística (4-6 refs):**
   - Referências para ANOVA, testes post-hoc, effect sizes
   - Critério: Livros clássicos (Fisher, Tukey, Cohen) + artigos metodológicos
   
5. **Frameworks Computacionais (3-5 refs):**
   - Referências para bibliotecas utilizadas (PennyLane, TensorFlow, etc.)
   - Critério: Artigo de apresentação da biblioteca
   
6. **Trabalhos Críticos/Opostos (3-5 refs):**
   - Artigos com visão contrária ou crítica à abordagem
   - Critério: Argumentos bem fundamentados, periódicos respeitáveis
   
7. **Aplicações e Implicações (3-5 refs):**
   - Trabalhos sobre aplicações práticas da área
   - Critério: Relevância para discussão de implicações

**Total de Referências Esperado:** 35-50 (padrão Qualis A1: 40-60)

**Formato de Saída:**
- Documento: `referencias_compiladas.md`
- Para cada referência: Título, Autores, Ano, Periódico, DOI, Citação ABNT completa
- Anotações sobre relevância de cada referência
```

#### 2.2 Análise e Síntese da Literatura

**Instrução:**
```
Para cada categoria de referências, crie uma síntese crítica:

1. **Identificar Consensos:**
   - Quais autores concordam sobre X?
   - Qual é a visão dominante?
   
2. **Identificar Divergências:**
   - Quais autores discordam sobre Y?
   - Quais são os pontos de debate?
   
3. **Identificar Lacunas:**
   - O que ainda não foi investigado?
   - Quais limitações os autores reconhecem?
   
4. **Posicionar a Contribuição:**
   - Como o presente estudo se relaciona com cada trabalho?
   - Quais lacunas específicas ele preenche?

**Formato de Saída:**
- Documento: `sintese_literatura.md`
- Tabela comparativa de abordagens
- Diagrama de posicionamento (este estudo vs. estado da arte)
```

---

### FASE 3: ELABORAÇÃO DA ESTRUTURA DO ARTIGO

#### 3.1 Definição de Título e Palavras-Chave

**Instrução:**
```
Gere 3 opções de título seguindo os critérios:

1. **Título Direto (estilo Nature/Science):**
   - Máximo 15 palavras
   - Foco no achado principal
   - Ex: "Quantum Noise Enhances Machine Learning Performance"
   
2. **Título Técnico (estilo Physical Review):**
   - Máximo 20 palavras
   - Inclui metodologia e resultado
   - Ex: "Beneficial Quantum Noise in Variational Classifiers: A Systematic Investigation of Dynamic Schedules"
   
3. **Título com Subtítulo (estilo híbrido):**
   - Título principal (8-10 palavras) + Subtítulo (6-8 palavras)
   - Ex: "From Obstacle to Opportunity: Harnessing Beneficial Quantum Noise in Variational Classifiers"

**Palavras-Chave (5-7):**
- 3 palavras gerais da área
- 2-3 palavras específicas da contribuição
- 1-2 palavras metodológicas

**Formato de Saída:**
- Documento: `titulos_palavras_chave.md`
- Análise de otimização para indexação (Google Scholar, Web of Science)
```

#### 3.2 Estruturação de Hipóteses e Objetivos

**Instrução:**
```
Formule hipóteses e objetivos seguindo o framework SMART:

**Hipótese Principal (H₀):**
- Declaração clara e testável
- Baseada na lacuna identificada
- Formato: "Se X, então Y, porque Z"
- Ex: "Se ruído quântico moderado for introduzido em VQCs, então a acurácia de generalização aumentará, porque o ruído atua como regularizador estocástico"

**Hipóteses Derivadas (H₁, H₂, H₃, H₄):**
- 4 hipóteses específicas testáveis
- Cada uma aborda uma dimensão da lacuna
- Formato: "H₁: [Fator X] resultará em [Efeito Y] com [Magnitude Z]"

**Objetivos Específicos (4-6):**
- Cada objetivo corresponde a uma hipótese
- Formato SMART:
  - **S**pecific: Claramente definido
  - **M**easurable: Métricas quantitativas especificadas
  - **A**chievable: Viável com a metodologia
  - **R**elevant: Alinhado com a lacuna
  - **T**ime-bound: Escopo delimitado

**Formato de Saída:**
- Documento: `hipoteses_objetivos.md`
- Tabela de alinhamento: Lacuna → Hipótese → Objetivo → Métrica
```

---

### FASE 4: REDAÇÃO DAS SEÇÕES DO ARTIGO

#### 4.1 RESUMO E ABSTRACT (250-300 palavras)

**Instrução:**
```
Redija Resumo (português) e Abstract (inglês) seguindo estrutura IMRAD:

**Estrutura (proporções):**
1. **Introdução/Contexto (15%):** Problema e relevância
2. **Métodos (35%):** Desenho do estudo, amostra, técnicas
3. **Resultados (40%):** Achados principais com dados quantitativos
4. **Conclusão (10%):** Implicação principal

**Requisitos:**
- Autocontido (faz sentido sozinho)
- Sem citações, figuras ou tabelas
- Dados quantitativos incluídos (ex: "+10.1% de melhoria")
- Voz ativa preferencial
- Palavras-chave integradas naturalmente

**Formato de Saída:**
- Documento: `resumo_abstract.md`
- Contagem de palavras para cada seção IMRAD
- Versão em português e inglês
```

#### 4.2 INTRODUÇÃO (3.000-4.000 palavras)

**Instrução:**
```
Redija a Introdução seguindo o modelo CARS (Create a Research Space):

**PASSO 1: Estabelecer o Território (30% da introdução)**
- Parágrafo 1: Contexto amplo da área (ex: Era NISQ)
- Parágrafo 2: Relevância multidisciplinar e aplicações
- Parágrafo 3: Visão tradicional/paradigma anterior

**PASSO 2: Estabelecer o Nicho (50% da introdução)**
- Parágrafo 4: Mudança de paradigma (precedentes clássicos)
- Parágrafo 5-6: Trabalho fundacional da linha de pesquisa
- Parágrafo 7-8: Extensões recentes e estado da arte
- Parágrafo 9-11: Lacuna em 3 dimensões (generalidade, dinâmica, interação)
- Parágrafo 12: Questão de pesquisa explícita

**PASSO 3: Ocupar o Nicho (20% da introdução)**
- Parágrafo 13: Hipótese principal (H₀)
- Parágrafo 14-17: Hipóteses derivadas (H₁-H₄)
- Parágrafo 18-21: Objetivos específicos (1-4)
- Parágrafo 22-23: Contribuições (teóricas, metodológicas, práticas)

**Requisitos de Formatação:**
- Citações narrativas: Autor (ano, p. página)
- Citações parentéticas: (AUTOR, ano, p. página)
- Parágrafos de 5-6 frases cada
- Transições fluidas entre parágrafos
- Linguagem equilibrada (rigorosa mas acessível)

**Formato de Saída:**
- Documento: `introducao_completa.md`
- Estrutura CARS claramente marcada
- 15-20 referências citadas
```

#### 4.3 REVISÃO DE LITERATURA (4.000-5.000 palavras)

**Instrução:**
```
Redija a Revisão de Literatura com estrutura temática e diálogo crítico:

**Estrutura Temática (7 seções):**

1. **Contexto Histórico e Paradigma Anterior (500 palavras)**
   - Evolução da área até o paradigma tradicional
   - Autores seminais e trabalhos fundacionais
   
2. **Problema Central (600 palavras)**
   - Descrição do problema (ex: Barren Plateaus)
   - Debate entre autores (McClean vs. Grant vs. Anschuetz)
   - Perspectivas divergentes
   
3. **Arquiteturas/Modelos (700 palavras)**
   - Revisão de diferentes abordagens (ex: Ansätze)
   - Trade-offs (expressividade vs. trainability)
   - Comparação crítica (Du, Kandala, Leone, Schuld)
   
4. **Técnica Central (800 palavras)**
   - Fundamentação teórica (ex: Equação de Lindblad)
   - Implementações na literatura
   - Conexões interdisciplinares (ex: Ressonância Estocástica)
   
5. **Otimização e Treinamento (600 palavras)**
   - Algoritmos de otimização (Adam, SGD, QNG)
   - Debate sobre eficiência (Kingma vs. Sweke vs. Stokes)
   
6. **Análise Estatística (500 palavras)**
   - Métodos estatísticos relevantes
   - Referências clássicas (Fisher, Tukey, Cohen)
   
7. **Frameworks Computacionais (300 palavras)**
   - Justificativa das bibliotecas escolhidas
   - Comparação de alternativas

**Requisitos de Diálogo Crítico:**
- Síntese: Conectar ideias de diferentes autores
- Comparação: "Enquanto Autor X argumenta..., Autor Y demonstra..."
- Contraste: "Em contraste com a visão de Z, este estudo..."
- Avaliação: Analisar limitações e validade de cada abordagem

**Formato de Saída:**
- Documento: `revisao_literatura_completa.md`
- 30-40 referências citadas
- Tabela comparativa de abordagens
```

#### 4.4 METODOLOGIA (4.000-5.000 palavras)

**Instrução:**
```
Redija a Metodologia com rigor máximo e reprodutibilidade garantida:

**Estrutura (11 subseções):**

1. **Desenho do Estudo (300 palavras)**
   - Tipo de estudo (experimental, observacional, etc.)
   - Questão de pesquisa reiterada
   - Três pilares teóricos

2. **Framework Computacional (400 palavras)**
   - Bibliotecas e versões exatas
   - Justificativa de cada escolha (com citação)
   - Ambiente de execução (hardware, SO)

3. **Datasets (400 palavras)**
   - Descrição de cada dataset (tamanho, características)
   - Pré-processamento (passos detalhados)
   - Divisão treino/validação/teste

4. **Arquiteturas/Modelos (600 palavras)**
   - Descrição de cada arquitetura implementada
   - Equações matemáticas (LaTeX)
   - Justificativa da diversidade

5. **Técnica Central (700 palavras)**
   - Fundamentação teórica (ex: Equação de Lindblad)
   - Implementação computacional
   - Tipos/variações explorados

6. **Inovação Metodológica (500 palavras)**
   - Descrição da contribuição original (ex: Schedules Dinâmicos)
   - Equações matemáticas
   - Justificativa teórica

7. **Inicialização (400 palavras)**
   - Estratégias implementadas
   - Justificativa (mitigação de barren plateaus)

8. **Otimização (400 palavras)**
   - Algoritmos utilizados
   - Hiperparâmetros
   - Critérios de convergência

9. **Análise Estatística (600 palavras)**
   - ANOVA multifatorial (equações)
   - Testes post-hoc (Tukey, Bonferroni, Scheffé)
   - Tamanhos de efeito (Cohen's d, Glass's Δ, Hedges' g)
   - Correção para comparações múltiplas

10. **Configurações Experimentais (300 palavras)**
    - Tabela de fatores e níveis
    - Cálculo do total de configurações
    - Seeds aleatórias e repetições

11. **Reprodutibilidade (200 palavras)**
    - Disponibilidade de código e dados
    - Instruções de instalação
    - Versões de software

**Requisitos Matemáticos:**
- Todas as equações em LaTeX
- Parágrafo explicativo para cada equação
- Notação consistente

**Formato de Saída:**
- Documento: `metodologia_completa.md`
- 15-20 referências citadas
- Tabela de configurações experimentais
```

#### 4.5 RESULTADOS (3.000-4.000 palavras)

**Instrução:**
```
Redija a seção de Resultados com análise estatística profunda:

**Estrutura (8 subseções):**

1. **Estatísticas Descritivas Gerais (300 palavras)**
   - Visão panorâmica dos experimentos
   - Média, desvio padrão, intervalos de confiança
   - Tabela resumo

2. **ANOVA Multifatorial (500 palavras)**
   - Tabela de ANOVA (fatores, F, p, η²)
   - Identificação de fatores significativos
   - Interpretação de tamanhos de efeito

3. **Teste de H₁ (400 palavras)**
   - Resultados específicos para H₁
   - Testes post-hoc relevantes
   - Tamanhos de efeito (Cohen's d)
   - Tabela de comparações

4. **Teste de H₂ (400 palavras)**
   - Resultados específicos para H₂
   - Análise estatística
   - Visualização (referência a figura)

5. **Teste de H₃ (400 palavras)**
   - Resultados de interação
   - Heatmap (referência a figura)
   - Interpretação

6. **Teste de H₄ (400 palavras)**
   - Resultados (incluindo refutação parcial, se aplicável)
   - Explicação de resultados inesperados
   - Implicações

7. **Análise de Sensibilidade (400 palavras)**
   - Curva de sensibilidade (referência a figura)
   - Identificação de regime ótimo
   - Tabela de valores

8. **Comparação por Dataset (300 palavras)**
   - Desempenho em cada dataset
   - Tabela comparativa

**Requisitos:**
- Apenas apresentação de resultados (sem interpretação)
- Todos os p-valores reportados
- Intervalos de confiança de 95%
- Referências a tabelas e figuras
- Resultados negativos incluídos

**Formato de Saída:**
- Documento: `resultados_completo.md`
- 7 tabelas detalhadas
- Referências a 2 figuras principais
```

#### 4.6 DISCUSSÃO (4.000-5.000 palavras)

**Instrução:**
```
Redija a Discussão com interpretação profunda e diálogo com literatura:

**Estrutura (6 subseções temáticas):**

1. **Síntese dos Achados (500 palavras)**
   - Reafirmação dos resultados principais
   - Resposta às hipóteses (H₁-H₄)
   - Mensagem central ("take-home message")

2. **Interpretação de H₁ e H₂ (800 palavras)**
   - Explicação dos mecanismos subjacentes
   - Comparação com Du et al., Liu et al., Choi et al.
   - Convergências e divergências
   - Implicações teóricas

3. **Interpretação de H₃ e H₄ (800 palavras)**
   - Explicação de interações
   - Discussão de refutação parcial (se aplicável)
   - Comparação com McClean, Grant, Skolik
   - Reconciliação de contradições

4. **Implicações Teóricas e Práticas (900 palavras)**
   - Mudança de paradigma (de "eliminação" para "engenharia")
   - Comparação com Preskill, Cerezo
   - Aplicações práticas (design de VQCs)
   - Escalabilidade e viabilidade

5. **Limitações (500 palavras)**
   - Simulação vs. hardware real
   - Escala (número de qubits)
   - Generalização para outros problemas
   - Discussão honesta e construtiva

6. **Trabalhos Futuros (500 palavras)**
   - Validação em hardware quântico (IBM, Google, Rigetti)
   - Estudos de escalabilidade
   - Teoria rigorosa (formalização matemática)
   - Ruído aprendível (otimização de parâmetros de ruído)

**Requisitos de Diálogo:**
- Comparar e contrastar com literatura
- Propor explicações para contradições
- Debater com autores de visões distintas
- Antecipar críticas de revisores

**Formato de Saída:**
- Documento: `discussao_completa.md`
- 20-25 referências citadas
- Linguagem dissertativa (sem tópicos)
```

#### 4.7 CONCLUSÃO (1.000-1.500 palavras)

**Instrução:**
```
Redija a Conclusão com síntese final e visão prospectiva:

**Estrutura (5 parágrafos):**

1. **Reafirmação do Problema e Objetivos (200 palavras)**
   - Contexto brevemente reiterado
   - Objetivos do estudo

2. **Síntese dos Principais Achados (300 palavras)**
   - Mensagem central ("take-home message")
   - Resultados mais significativos (quantitativos)
   - Confirmação/refutação de hipóteses

3. **Contribuições Originais (300 palavras)**
   - Teóricas: Generalização do fenômeno
   - Metodológicas: Inovação específica (ex: Schedules Dinâmicos)
   - Práticas: Diretrizes para aplicação

4. **Limitações e Visão Futura (300 palavras)**
   - Limitações mais significativas
   - Próximos passos da pesquisa (3-4 específicos)
   - Fronteira de pesquisa

5. **Declaração Final Forte (100 palavras)**
   - Frase impactante e memorável
   - Reforça a mensagem principal
   - Ex: "A era da engenharia do ruído quântico apenas começou"

**Requisitos:**
- Sem novas informações ou citações
- Linguagem concisa e impactante
- Foco em implicações amplas

**Formato de Saída:**
- Documento: `conclusao_completa.md`
```

#### 4.8 AGRADECIMENTOS E REFERÊNCIAS

**Instrução:**
```
**Agradecimentos (150-200 palavras):**
- Agências de fomento (com números de projeto)
- Instituições (universidade, grupo de pesquisa)
- Colaboradores (sem coautoria)
- Infraestrutura computacional
- Revisores (após aceite)

**Referências (35-50 entradas):**
- Formato ABNT rigoroso
- Ordem alfabética por sobrenome
- 100% de correspondência com citações no texto
- DOI para 90%+ das referências
- Verificação de rastreabilidade

**Formato de Saída:**
- Documento: `agradecimentos_referencias.md`
- Checklist de correspondência citação-referência
```

---

### FASE 5: MATERIAL SUPLEMENTAR

#### 5.1 Tabelas Suplementares

**Instrução:**
```
Crie 5 tabelas suplementares:

**Tabela S1: Configurações Experimentais Completas**
- Arquivo CSV com todas as configurações
- Colunas: Fatores + Métricas de desempenho
- Total de linhas = número de configurações

**Tabela S2: Comparação com Estado da Arte**
- Comparação com 5-10 estudos relevantes
- Colunas: Estudo, Ano, Dataset, Método, Resultado, Referência
- Destaque para melhorias alcançadas

**Tabela S3: Análise de Custo Computacional**
- Configurações representativas
- Colunas: Configuração, Qubits, Portas, Profundidade, Tempo, Memória

**Tabela S4: Análise Estatística Detalhada**
- Todos os testes post-hoc
- Colunas: Comparação, Diferença, p-valor, IC 95%, Cohen's d

**Tabela S5: Análise de Sensibilidade**
- Variação de parâmetro-chave
- Colunas: Parâmetro, Métrica1, Métrica2, ..., Melhoria vs. Baseline

**Formato de Saída:**
- Documento: `tabelas_suplementares.md`
- Arquivo CSV: `tabela_s1_configuracoes.csv`
```

#### 5.2 Figuras Suplementares

**Instrução:**
```
Descreva 6-8 figuras suplementares:

**Para cada figura:**
1. **Título descritivo**
2. **Tipo:** (ex: Curva de convergência, Heatmap, Scatter plot)
3. **Descrição detalhada (100-150 palavras):**
   - O que é mostrado (eixos, cores, símbolos)
   - Número de subplots (se aplicável)
   - Escala (linear, log)
   - Marcadores especiais
4. **Achados-chave (50-80 palavras):**
   - Principais observações
   - Valores quantitativos destacados

**Figuras Recomendadas:**
- S1: Curvas de convergência por condição
- S2: Heatmap de interações
- S3: Curva de sensibilidade
- S4: Distribuição de gradientes (barren plateaus)
- S5: PCA do espaço de parâmetros
- S6: Análise de poder estatístico
- S7: Interações de ordem superior
- S8: Custo computacional vs. desempenho (Pareto front)

**Formato de Saída:**
- Documento: `figuras_suplementares.md`
- Especificações técnicas (formato, resolução, paleta de cores)
```

#### 5.3 Notas Metodológicas Adicionais

**Instrução:**
```
Inclua detalhes metodológicos adicionais:

1. **Detalhes de Implementação (300 palavras)**
   - Versões exatas de bibliotecas
   - Configuração de hardware
   - Paralelização

2. **Critérios de Convergência (200 palavras)**
   - Condições específicas
   - Tratamento de não-convergência

3. **Tratamento de Outliers (200 palavras)**
   - Critérios de exclusão
   - Número de configurações excluídas

4. **Validação Cruzada (200 palavras)**
   - Estratégia (k-fold, stratified, etc.)
   - Número de repetições
   - Seeds aleatórias

**Formato de Saída:**
- Documento: `notas_metodologicas_adicionais.md`
```

---

### FASE 6: CONSOLIDAÇÃO E VERIFICAÇÃO

#### 6.1 Verificação de Conivência Código-Texto

**Instrução:**
```
Verifique 100% de conivência entre código e texto:

**Checklist de Verificação:**

1. **Componentes Técnicos:**
   - [ ] Número de arquiteturas/modelos: Código = Texto
   - [ ] Número de técnicas: Código = Texto
   - [ ] Número de datasets: Código = Texto
   - [ ] Versões de bibliotecas: Código = Texto

2. **Configurações Experimentais:**
   - [ ] Total de configurações: Código = Texto
   - [ ] Fatores experimentais: Código = Texto
   - [ ] Níveis de cada fator: Código = Texto

3. **Métricas e Resultados:**
   - [ ] Métricas de avaliação: Código = Texto
   - [ ] Testes estatísticos: Código = Texto
   - [ ] Valores reportados: Código = Texto (ou simulados consistentemente)

4. **Inovações Metodológicas:**
   - [ ] Técnicas originais: Código = Texto
   - [ ] Equações matemáticas: Código = Texto

**Formato de Saída:**
- Documento: `relatorio_conivencia.md`
- Percentual de conivência (meta: ≥95%)
- Lista de inconsistências (se houver)
- Recomendações de ajustes
```

#### 6.2 Consolidação do Artigo Completo

**Instrução:**
```
Consolide todas as seções em um documento único:

**Estrutura do Artigo Completo:**

1. Título
2. Autores e Afiliações (placeholders)
3. Resumo (português)
4. Abstract (inglês)
5. Palavras-chave (português e inglês)
6. 1 INTRODUÇÃO
7. 2 REVISÃO DE LITERATURA
8. 3 METODOLOGIA
9. 4 RESULTADOS
10. 5 DISCUSSÃO
11. 6 CONCLUSÃO
12. AGRADECIMENTOS
13. REFERÊNCIAS

**Verificações Finais:**
- [ ] Numeração de seções consistente
- [ ] Referências cruzadas corretas (Tabela X, Figura Y)
- [ ] Formatação ABNT rigorosa
- [ ] Contagem de palavras (meta: 10.000-12.000)
- [ ] Número de referências (meta: 35-50)
- [ ] Transições fluidas entre seções

**Formato de Saída:**
- Documento: `artigo_completo_final.md`
- Sumário executivo: `sumario_executivo.md`
```

#### 6.3 Criação de Sumário Executivo

**Instrução:**
```
Crie um sumário executivo com estatísticas do artigo:

**Conteúdo:**

1. **Estatísticas Gerais:**
   - Extensão total (palavras)
   - Número de seções/subseções
   - Número de tabelas (principais + suplementares)
   - Número de figuras (principais + suplementares)
   - Número de referências

2. **Destaques Metodológicos:**
   - Configurações experimentais
   - Datasets utilizados
   - Técnicas estatísticas
   - Inovações originais

3. **Principais Achados (3-5 bullets):**
   - Resultados quantitativos mais significativos

4. **Contribuições (3 níveis):**
   - Teóricas
   - Metodológicas
   - Práticas

5. **Conformidade Qualis A1:**
   - Checklist de critérios atendidos

6. **Periódicos-Alvo Recomendados:**
   - 3-5 periódicos com justificativa

**Formato de Saída:**
- Documento: `sumario_executivo.md`
```

---

## 📊 CRITÉRIOS DE AVALIAÇÃO QUALIS A1

### Checklist de Conformidade

**Rigor Acadêmico e Profundidade:**
- [ ] Análise estatística profunda (ANOVA, post-hoc, effect sizes)
- [ ] Correção para comparações múltiplas (Bonferroni)
- [ ] Intervalos de confiança (95%) para todas as médias
- [ ] Análise de poder estatístico
- [ ] Validação cruzada implementada

**Relevância e Implicações:**
- [ ] Implicações teóricas claramente articuladas
- [ ] Implicações práticas discutidas
- [ ] Escalabilidade e viabilidade abordadas
- [ ] Comparação com estado da arte (tabela)

**Clareza, Organização e Didática:**
- [ ] Transições fluidas entre seções
- [ ] Linguagem equilibrada (rigorosa mas acessível)
- [ ] Parágrafos bem desenvolvidos (5-6 frases)
- [ ] Estrutura lógica e coesa

**Metodologia:**
- [ ] Metodologia detalhada e replicável
- [ ] Justificativas para cada escolha metodológica
- [ ] Limitações conhecidas discutidas
- [ ] Código e dados disponíveis (repositório público)

**Revisão da Literatura:**
- [ ] Revisão abrangente (35-50 referências)
- [ ] Síntese crítica (não lista de resumos)
- [ ] Diálogo com visões opostas
- [ ] Lacuna claramente identificada

**Resultados e Discussão:**
- [ ] Resultados com análises estatísticas profundas
- [ ] Discussão interpretando achados
- [ ] Comparação com literatura existente
- [ ] Limitações honestamente discutidas
- [ ] Trabalhos futuros específicos propostos

**Qualidade de Figuras, Tabelas e Fórmulas:**
- [ ] Tabelas detalhadas (4-5 principais + 5 suplementares)
- [ ] Figuras conceituais (2-3 principais + 6-8 suplementares)
- [ ] Legendas descritivas e expandidas
- [ ] Fórmulas matemáticas em LaTeX
- [ ] Parágrafo explicativo para cada fórmula

**Aderência às Normas:**
- [ ] Formato ABNT rigoroso para citações
- [ ] 100% de correspondência citação-referência
- [ ] DOI/URL válidos para 90%+ das referências

---

## 🔄 REPRODUTIBILIDADE E AUDITORIA

### Documentação para Auditoria

**Arquivos a Serem Entregues:**

1. **Análise Inicial:**
   - `analise_codigo_inicial.md`
   - `linha_de_pesquisa.md`

2. **Pesquisa Bibliográfica:**
   - `referencias_compiladas.md`
   - `sintese_literatura.md`

3. **Estruturação:**
   - `titulos_palavras_chave.md`
   - `hipoteses_objetivos.md`

4. **Seções do Artigo:**
   - `resumo_abstract.md`
   - `introducao_completa.md`
   - `revisao_literatura_completa.md`
   - `metodologia_completa.md`
   - `resultados_completo.md`
   - `discussao_completa.md`
   - `conclusao_completa.md`
   - `agradecimentos_referencias.md`

5. **Material Suplementar:**
   - `tabelas_suplementares.md`
   - `tabela_s1_configuracoes.csv`
   - `figuras_suplementares.md`
   - `notas_metodologicas_adicionais.md`

6. **Consolidação:**
   - `relatorio_conivencia.md`
   - `artigo_completo_final.md`
   - `sumario_executivo.md`

7. **Código e Dados:**
   - Repositório GitHub com código-fonte completo
   - Arquivo README.md com instruções
   - Arquivo de ambiente (conda/pip)
   - Dados brutos (se aplicável)

### Rastreabilidade

**Para cada seção do artigo, documentar:**
1. **Fonte dos dados:** Código linha X, função Y
2. **Referências utilizadas:** Lista de [Autor, Ano]
3. **Decisões metodológicas:** Justificativa com citação
4. **Valores reportados:** Origem (código, simulação, literatura)

**Formato:**
- Documento: `rastreabilidade_completa.md`
- Tabela de rastreabilidade: Seção → Fonte → Referência

---

## 🎓 EXEMPLO DE USO DO MEGA PROMPT

### Caso de Uso: Artigo sobre Beneficial Quantum Noise

**Input:**
```
Código: framework_investigativo_completo.py (2.845 linhas)
Repositório: https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers
```

**Execução:**
```
FASE 1: Análise do código
- Identificados: 11 módulos, 19 classes, 9 arquiteturas, 5 modelos de ruído
- Total de configurações: 4 × 7 × 6 × 2 × 8 = 2.688

FASE 2: Pesquisa bibliográfica
- 38 referências compiladas (Du 2021, Liu 2025, McClean 2018, etc.)
- Linha de pesquisa: Du et al. 2021 - Ruído Benéfico
- Lacuna: Generalidade, Dinâmica, Interação

FASE 3: Estruturação
- Título: "From Obstacle to Opportunity: Harnessing Beneficial Quantum Noise in Variational Classifiers"
- Hipótese principal: Ruído moderado melhora generalização
- 4 hipóteses derivadas, 4 objetivos SMART

FASE 4: Redação
- Introdução: 3.500 palavras (modelo CARS)
- Revisão: 4.000 palavras (7 seções temáticas)
- Metodologia: 4.500 palavras (11 subseções)
- Resultados: 3.500 palavras (8 subseções)
- Discussão: 4.500 palavras (6 subseções)
- Conclusão: 1.200 palavras

FASE 5: Material Suplementar
- 5 tabelas suplementares
- 8 figuras suplementares descritas
- Notas metodológicas adicionais

FASE 6: Consolidação
- Conivência: 100% (após ajustes)
- Artigo completo: 11.249 palavras
- Referências: 38
- Sumário executivo criado
```

**Output:**
- Artigo completo pronto para submissão a Nature/Quantum/Physical Review
- Material Suplementar completo
- Relatório de conivência 100%
- Todos os arquivos de auditoria

---

## 📝 NOTAS FINAIS

### Adaptações por Área de Pesquisa

Este mega prompt é genérico e deve ser adaptado para áreas específicas:

**Computação Quântica:**
- Ênfase em ansätze, barren plateaus, ruído quântico
- Referências: Preskill, Cerezo, McClean, Grant

**Machine Learning:**
- Ênfase em arquiteturas de redes, regularização, generalização
- Referências: Goodfellow, LeCun, Bengio

**Bioinformática:**
- Ênfase em pipelines de análise, validação biológica
- Referências: Altschul, Bairoch, Consortium

**Física Teórica:**
- Ênfase em formalismo matemático, derivações rigorosas
- Referências: Landau, Feynman, Weinberg

### Periódicos-Alvo por Área

**Alto Impacto Multidisciplinar:**
- Nature, Science, PNAS

**Computação Quântica:**
- Quantum, npj Quantum Information, Physical Review A/Research

**Machine Learning:**
- Nature Machine Intelligence, JMLR, NeurIPS

**Física:**
- Physical Review Letters, Physical Review X

**Interdisciplinar:**
- Nature Communications, Scientific Reports

---

## ✅ VALIDAÇÃO FINAL

Antes de submeter, verificar:

- [ ] Conivência código-texto ≥ 95%
- [ ] Formato ABNT 100% correto
- [ ] Todas as hipóteses testadas
- [ ] Todos os objetivos atingidos
- [ ] Limitações honestamente discutidas
- [ ] Trabalhos futuros específicos propostos
- [ ] Código e dados disponíveis publicamente
- [ ] Material Suplementar completo
- [ ] Sumário executivo criado
- [ ] Todos os arquivos de auditoria entregues

---


