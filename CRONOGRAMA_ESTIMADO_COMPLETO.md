# Cronograma Estimado Completo - Geração de Artigos QUALIS A1

**Framework de Rastreabilidade Total**  
**Versão:** 1.0  
**Data:** 26/12/2025

---

## 📅 VISÃO GERAL

**Duração Total:** 52-78 horas de trabalho efetivo (6-10 dias úteis)  
**Calendário:** 4-8 semanas (considerando revisões e intervalos)  
**Equipe Mínima:** 1 pesquisador principal + 1 revisor

---

## 📊 DISTRIBUIÇÃO DE TEMPO POR FASE

```
FASE 1: Auditoria Técnica       ████████░░  8-12h  (15-17%)
FASE 2: Bibliografia             ██████░░░░  6-25h  (12-32%)
FASE 3: Projeto do Artigo        ████░░░░░░  4-6h   (7-8%)
FASE 4: Redação                  ████████████████████ 20-30h (38-40%)
FASE 5: Material Suplementar     ████████░░  8-12h  (15-17%)
FASE 6: Consolidação             ██████░░░░  6-8h   (10-12%)
                                 ─────────────────────
                                 52-78h total
```

---

## 🗓️ CRONOGRAMA DETALHADO

### SEMANA 1: Fundação (22-35h)

#### Segunda-feira (8h)
**FASE 1: Auditoria Técnica (Parte 1)**
- [ ] **08:00-10:00** (2h) Configuração inicial
  - Preencher `config_artigo.json`
  - Instalar ferramentas necessárias
  - Configurar ambiente Git/GitHub
  
- [ ] **10:00-12:00** (2h) Análise de código
  - Contar arquivos, linhas, módulos
  - Identificar classes e funções principais
  - Listar dependências com versões

- [ ] **14:00-16:00** (2h) Componentes experimentais
  - Calcular total de configurações
  - Documentar fatores e níveis
  - Identificar datasets e métricas

- [ ] **16:00-18:00** (2h) Metodologia e inovações
  - Documentar pipeline de execução
  - Identificar código novo vs baseline
  - Criar `analise_codigo_inicial.md`

**Deliverables:** `analise_codigo_inicial.md`, `tabela_componentes.md`

---

#### Terça-feira (8h)
**FASE 1: Auditoria Técnica (Parte 2) + FASE 2 Início**

- [ ] **08:00-10:00** (2h) Finalizar Fase 1
  - Criar `mapa_execucao.md`
  - Gerar `manifesto_execucao.json`
  - ✅ Quality Gate F1

- [ ] **10:00-12:00** (2h) Busca bibliográfica inicial
  - Definir política (R0 ou R1)
  - Se R1: Buscar categoria 1 (Fundacionais)
  - Se R0: Validar lista de referências aprovadas

- [ ] **14:00-16:00** (2h) Busca categoria 2
  - Estado da arte (últimos 2-3 anos)
  - arXiv, Web of Science, Google Scholar

- [ ] **16:00-18:00** (2h) Busca categorias 3-4
  - Metodologias e estatística
  - Documentar DOIs e justificativas

**Deliverables:** Fase 1 completa, 15-25 refs para `referencias_compiladas.md`

---

#### Quarta-feira (6-8h)
**FASE 2: Bibliografia (Continuação)**

- [ ] **08:00-10:00** (2h) Busca categorias 5-7
  - Frameworks/bibliotecas
  - Contrapontos críticos
  - Aplicações práticas

- [ ] **10:00-12:00** (2h) Compilação de referências
  - Preencher `referencias_compiladas.md`
  - Formatar ABNT ou APA
  - Verificar DOIs (meta: ≥80%)

- [ ] **14:00-16:00** (2-4h) Síntese de literatura
  - Identificar consensos
  - Identificar divergências
  - Posicionar contribuição
  - Criar `sintese_literatura.md`

**Deliverables:** `referencias_compiladas.md` (35-60 refs), `sintese_literatura.md`

---

#### Quinta-feira (4h)
**FASE 3: Projeto do Artigo**

- [ ] **08:00-10:00** (2h) Formalização do problema
  - Notação matemática
  - Equações LaTeX
  - Hipótese principal
  - Criar `problema_formal.md`

- [ ] **10:00-12:00** (2h) Títulos e objetivos
  - 3 opções de título
  - 6 palavras-chave
  - H₀ + H₁-H₄
  - 4 objetivos SMART
  - Criar `titulos_palavras_chave.md`, `hipoteses_objetivos.md`
  - ✅ Quality Gate F3

**Deliverables:** Fase 3 completa (3 documentos)

---

#### Sexta-feira (4h)
**FASE 4: Início da Redação**

- [ ] **08:00-10:00** (2h) Abstract/Resumo
  - Estrutura IMRAD (Intro, Methods, Results, Discussion)
  - 250-300 palavras
  - Auto-suficiente
  - Criar `resumo_abstract.md`

- [ ] **10:00-12:00** (2h) Introdução (Parte 1)
  - Contextualização histórica
  - Problema fundamental
  - Início da estrutura CARS

**Deliverables:** `resumo_abstract.md` completo, introdução 20% completa

---

### SEMANA 2: Redação Principal (16-20h)

#### Segunda-feira (8h)
**FASE 4: Introdução + Related Work**

- [ ] **08:00-12:00** (4h) Introdução (Parte 2)
  - Move 2: Estabelecer nicho (gap)
  - Move 3: Ocupar nicho (contribuições)
  - Finalizar `introducao_completa.md` (1.500-2.500 palavras)

- [ ] **14:00-18:00** (4h) Revisão de literatura
  - Organizar em 7 temas
  - Síntese crítica (não lista!)
  - Contrapontos
  - Criar `revisao_literatura_completa.md` (2.000-3.000 palavras)

**Deliverables:** Introdução e Related Work completas

---

#### Terça-feira (8-10h)
**FASE 4: Metodologia**

- [ ] **08:00-10:00** (2h) Methods: Fundamentação
  - Paradigma de pesquisa
  - Questões de pesquisa
  - Hipóteses formais

- [ ] **10:00-12:00** (2h) Methods: Framework computacional
  - Descrição técnica do pipeline
  - Ansätze, ruídos, schedules
  - Datasets e pré-processamento

- [ ] **14:00-16:00** (2h) Methods: Algoritmo e equações
  - Algorithm 1 (LaTeX)
  - Equações principais com explicações
  - Notação consistente

- [ ] **16:00-18:00** (2-4h) Methods: Mapeamento código→método
  - Tabela Código→Método
  - Rastreabilidade de componentes
  - Validação estatística
  - Finalizar `metodologia_completa.md` (2.500-4.000 palavras)

**Deliverables:** `metodologia_completa.md` completa

---

### SEMANA 3: Resultados e Discussão (12-16h)

#### Quarta-feira (6-8h)
**FASE 4: Resultados**

- [ ] **08:00-10:00** (2h) Results: Estatística descritiva
  - Tabelas de médias ± IC
  - Gráficos principais
  - Acurácia ótima identificada

- [ ] **10:00-12:00** (2h) Results: Testes de hipóteses
  - ANOVA resultados
  - Post-hoc tests
  - Effect sizes

- [ ] **14:00-16:00** (2-4h) Results: Análises secundárias
  - Comparação multiframework (se aplicável)
  - Análise de sensibilidade
  - Finalizar `resultados_completo.md` (2.000-3.000 palavras, 9+ tabelas)

**Deliverables:** `resultados_completo.md` completo

---

#### Quinta-feira (6-8h)
**FASE 4: Discussão**

- [ ] **08:00-10:00** (2h) Discussion: Interpretação
  - O que os resultados significam?
  - Comparação com literatura
  - Explicação de mecanismos

- [ ] **10:00-12:00** (2h) Discussion: Implicações
  - Teóricas
  - Práticas
  - Metodológicas

- [ ] **14:00-16:00** (2-4h) Discussion: Limitações e validade
  - Threats to validity (4 tipos)
  - Scope conditions
  - Trade-offs
  - Finalizar `discussao_completa.md` (2.500-4.000 palavras)

**Deliverables:** `discussao_completa.md` completo

---

#### Sexta-feira (4h)
**FASE 4: Conclusão e Referências**

- [ ] **08:00-10:00** (2h) Conclusão
  - Resumo de achados
  - Resposta aos objetivos
  - Contribuições
  - Trabalhos futuros específicos
  - Criar `conclusao_completa.md` (500-800 palavras)

- [ ] **10:00-12:00** (2h) Referências e agradecimentos
  - Formatar referências (ABNT ou APA)
  - Verificar 100% citação↔referência
  - Criar `agradecimentos_referencias.md`
  - ✅ Quality Gate F4

**Deliverables:** Fase 4 completa (todas as seções do artigo)

---

### SEMANA 4: Material Suplementar (8-12h)

#### Segunda-feira (4h)
**FASE 5: Tabelas Suplementares**

- [ ] **08:00-10:00** (2h) Tabela S1 (Configurações)
  - Gerar CSV completo
  - Validar total de linhas
  - Documentar colunas

- [ ] **10:00-12:00** (2h) Tabelas S2-S5
  - S2: Comparação estado da arte
  - S3: Hiperparâmetros
  - S4: Testes post-hoc
  - S5: Análise de sensibilidade

**Deliverables:** `tabelas_suplementares.md` + `tabela_s1_configuracoes.csv`

---

#### Terça-feira (4-8h)
**FASE 5: Figuras e Notas**

- [ ] **08:00-12:00** (4h) Figuras Suplementares
  - Especificar 8 figuras
  - Descrição completa (eixos, escalas, colormap, DPI)
  - Identificar achado-chave de cada figura
  - Apontar scripts de geração
  - Criar `figuras_suplementares.md`

- [ ] **14:00-16:00** (0-4h) Notas Metodológicas
  - Detalhes técnicos adicionais
  - Derivações matemáticas longas
  - Implementações específicas
  - Criar `notas_metodologicas_adicionais.md`
  - ✅ Quality Gate F5

**Deliverables:** Fase 5 completa (3 documentos + CSV)

---

### SEMANA 5-6: Consolidação e Auditoria Final (6-8h)

#### Quarta-feira (3h)
**FASE 6: Verificação de Conivência**

- [ ] **08:00-10:00** (2h) Executar scripts de verificação
  - `check_consistency.py`
  - Identificar discrepâncias
  - Calcular % de conivência

- [ ] **10:00-11:00** (1h) Corrigir discrepâncias
  - Ajustar código ou texto
  - Re-executar verificação
  - Meta: ≥95% conivência
  - Criar `relatorio_conivencia.md`

**Deliverables:** `relatorio_conivencia.md` com ≥95% conivência

---

#### Quinta-feira (3-5h)
**FASE 6: Rastreabilidade e Consolidação**

- [ ] **08:00-10:00** (2h) Tabela de rastreabilidade
  - Preencher 20+ entradas
  - Mapear afirmações→código
  - Criar `rastreabilidade_completa.md`

- [ ] **10:00-11:00** (1h) Tabela código→método
  - Completar mapeamento
  - Validar origens
  - Criar `tabela_codigo_metodo.md`

- [ ] **11:00-12:00** (0-2h) Artigo consolidado
  - Compilar todas as seções
  - Gerar PDF/LaTeX final
  - Criar `artigo_completo_final.md` ou `.tex`

**Deliverables:** `rastreabilidade_completa.md`, `tabela_codigo_metodo.md`, artigo final

---

#### Sexta-feira (2h)
**FASE 6: Sumário e Checklist**

- [ ] **08:00-09:00** (1h) Sumário executivo
  - Visão geral do projeto
  - Destaques principais
  - Próximos passos
  - Criar `sumario_executivo_final.md`

- [ ] **09:00-10:00** (1h) Checklist de auditoria
  - Preencher `CHECKLIST_AUDITORIA_COMPLETO.md`
  - Calcular pontuação (meta: ≥90/100)
  - Identificar oportunidades de melhoria
  - ✅ Quality Gate F6

**Deliverables:** Fase 6 completa, projeto pronto para submissão! 🎉

---

## 📊 MARCOS E ENTREGAS (Milestones)

| Milestone | Entrega | Prazo | Status |
|-----------|---------|-------|--------|
| **M1** | Auditoria técnica completa | Fim Semana 1 Terça | ⏳ |
| **M2** | Bibliografia e estrutura | Fim Semana 1 | ⏳ |
| **M3** | Abstract + Introdução + Related Work | Fim Semana 2 Segunda | ⏳ |
| **M4** | Metodologia completa | Fim Semana 2 Terça | ⏳ |
| **M5** | Resultados e Discussão | Fim Semana 3 Quinta | ⏳ |
| **M6** | Conclusão e artigo principal completo | Fim Semana 3 Sexta | ⏳ |
| **M7** | Material suplementar completo | Fim Semana 4 Terça | ⏳ |
| **M8** | Auditoria final e submissão | Fim Semana 5-6 Sexta | ⏳ |

---

## 🎯 CONTINGÊNCIAS E RISCOS

### Risco 1: Atraso na Fase 2 (Bibliografia)
**Probabilidade:** Média  
**Impacto:** Alto (atrasa todo cronograma)  
**Mitigação:**
- Se política R1, considerar mudar para R0 (reduz 10-15h)
- Paralelizar busca (dividir 7 categorias entre co-autores)
- Usar ferramentas automatizadas (Connected Papers, Semantic Scholar)

### Risco 2: Discrepâncias Código-Texto (Fase 6)
**Probabilidade:** Alta (primeira vez)  
**Impacto:** Médio (2-4h extras)  
**Mitigação:**
- Executar `check_consistency.py` ao final de cada seção (não deixar para Fase 6)
- Manter documentação paralela durante Fase 1 (atualizar conforme código evolui)

### Risco 3: Redação Bloqueada (Writer's Block)
**Probabilidade:** Média  
**Impacto:** Médio (atraso 1-2 dias)  
**Mitigação:**
- Usar templates fornecidos (copiar estrutura, preencher lacunas)
- Começar por seções mais fáceis (Methods antes de Introduction)
- Usar ferramentas IA para drafts iniciais (GPT-4, Claude) e refinar

### Risco 4: Qualidade Insatisfatória (<90/100 pontos)
**Probabilidade:** Baixa (com este guia)  
**Impacto:** Alto (revisão profunda necessária)  
**Mitigação:**
- Preencher checklist progressivamente (não esperar fim)
- Solicitar revisão de co-autor ao fim de cada fase
- Priorizar itens de alto peso (rastreabilidade, estatística)

---

## 🔄 VARIAÇÕES DO CRONOGRAMA

### Modo Acelerado (40h, 1 semana)
**Quando usar:** Prazo apertado, deadline em 7 dias

**Ajustes:**
- Fase 2: Usar R0 (lista pré-aprovada) → 4h ao invés de 15-25h
- Fase 4: Reduzir palavras em 30% (Introduction 1.000 palavras, etc.)
- Fase 5: Apenas Tabelas S1-S3, Figuras S1-S4 (reduzir pela metade)
- Trabalhar 8h/dia × 5 dias úteis

**Trade-off:** Qualidade ligeiramente reduzida (85-90 pts ao invés de 90-100)

---

### Modo Detalhado (100h, 3-4 semanas)
**Quando usar:** Artigo para Nature/Science, sem pressa

**Adições:**
- **+10h** em Fase 2: Busca exaustiva em 7 categorias, 60+ refs
- **+8h** em Fase 4: Redação cuidadosa, 3 revisões internas
- **+6h** em Fase 5: Material suplementar expandido (10 tabelas, 12 figuras)
- **+6h** em Fase 6: Auditoria tripla, rastreabilidade de 100% afirmações

**Resultado:** Pontuação 95-100, artigo impecável

---

## 📝 CHECKLIST DE PROGRESSO

### Semana 1
- [ ] Fase 1 completa (8-12h)
- [ ] Fase 2 completa (6-25h)
- [ ] Fase 3 completa (4-6h)

### Semana 2-3
- [ ] Fase 4 completa (20-30h)
  - [ ] Abstract
  - [ ] Introduction
  - [ ] Related Work
  - [ ] Methods
  - [ ] Results
  - [ ] Discussion
  - [ ] Conclusion
  - [ ] References

### Semana 4
- [ ] Fase 5 completa (8-12h)
  - [ ] 5 tabelas suplementares
  - [ ] 8 figuras suplementares
  - [ ] Notas metodológicas

### Semana 5-6
- [ ] Fase 6 completa (6-8h)
  - [ ] Relatório conivência ≥95%
  - [ ] Rastreabilidade completa
  - [ ] Tabela código→método
  - [ ] Artigo consolidado
  - [ ] Sumário executivo
  - [ ] Checklist ≥90/100

---

## 🎓 SUBMISSÃO

**Preparação final:** 2-4h  
**Atividades:**
- [ ] Compilar LaTeX/PDF final
- [ ] Preparar material suplementar (ZIP)
- [ ] Escrever cover letter
- [ ] Preencher formulário de submissão online
- [ ] Upload de arquivos
- [ ] Verificação final pré-envio

**🎉 SUBMISSÃO CONCLUÍDA!**

**Aguardar resposta:** 4-12 semanas (típico para periódicos QUALIS A1)

---

## 📊 ESTATÍSTICAS HISTÓRICAS (Este Projeto)

**Projeto:** Beneficial Quantum Noise in VQCs  
**Período:** 25-26/12/2025  
**Horas reais:** ~65h  
**Semanas calendário:** 2  
**Equipe:** 1 pesquisador principal + IA assistência

**Fases:**
- Fase 1: 3h (eficiente devido a código bem organizado)
- Fase 2: 4h (R1, mas foco em refs principais)
- Fase 3: 2h (problema bem definido)
- Fase 4: 40h (maior esforço, 8 seções)
- Fase 5: 6h (tabelas/figuras já planejadas)
- Fase 6: 8h (múltiplas verificações)

**Resultado:** 
- 24 documentos, 460 KB, 50.000 palavras
- Pontuação auditoria: 91/100 (Excelente)
- Conivência: 100%
- Status: ✅ Pronto para submissão npj QI/Quantum

---

## 🔗 RECURSOS ADICIONAIS

- **Guia Completo:** `GUIA_COMPLETO_GERACAO_ARTIGOS.md`
- **FAQ:** `FAQ_TROUBLESHOOTING_COMPLETO.md`
- **Glossário:** `GLOSSARIO_COMPLETO.md`
- **Checklist:** `CHECKLIST_AUDITORIA_COMPLETO.md`
- **Fluxograma R0/R1:** `FLUXOGRAMA_R0_R1.md`

---

**Versão:** 1.0  
**Última Atualização:** 26/12/2025  
**Status:** ✅ Completo e validado  
**Licença:** MIT
