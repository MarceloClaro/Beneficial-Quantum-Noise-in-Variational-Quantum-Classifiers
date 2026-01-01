# Cronograma Estimado - Geração de Artigo Científico Qualis A1

## 📅 Visão Geral

Este documento fornece estimativas realistas de tempo para cada fase do processo de geração de artigos científicos com rastreabilidade total.

**Tempo Total Estimado**: 52-78 horas (6-10 dias úteis com dedicação de 8h/dia)


---


## ⏱️ Detalhamento por Fase

### Fase 1: Auditoria Técnica do Código/Dados
**Duração**: 8-12 horas


| Atividade | Tempo | Nível de Automação |
|-----------|-------|-------------------|
| Execução do `enhanced_code_analyzer.py` | 5-10 min | 🤖 Automático |
| Revisão e validação dos componentes extraídos | 2-3h | 👤 Manual |
| Cálculo e verificação de configurações experimentais | 1-2h | 👤 Manual |
| Documentação de metodologia implementada | 2-3h | 👤 Manual |
| Geração de `mapa_execucao.md` | 1-2h | 👤 Manual |
| Identificação de inovações e contribuições | 1-2h | 👤 Manual |
| Quality Gate F1 | 30min | 🤖 Semi-automático |

**Deliverables**:
- ✅ `analise_codigo_inicial.md`
- ✅ `tabela_componentes.md`
- ✅ `mapa_execucao.md`
- ✅ `manifesto_execucao.json`
- ✅ `code_analysis_report.json`


---


### Fase 2: Enquadramento Científico
**Duração**: 4-6 horas


| Atividade | Tempo | Nível de Automação |
|-----------|-------|-------------------|
| Identificação de área e subárea | 1h | 👤 Manual |
| Formulação do problema central | 1-2h | 👤 Manual |
| Identificação de lacuna em 3 dimensões | 1-2h | 👤 Manual |
| Geração de diagrama Mermaid | 30min | 🤖 Semi-automático |
| Alinhamento com objetivos/hipóteses | 30min | 👤 Manual |
| Quality Gate F2 | 30min | 🤖 Semi-automático |

**Deliverables**:
- ✅ `linha_de_pesquisa.md`
- ✅ `diagrama_linha_pesquisa.md`


---


### Fase 3: Curadoria Bibliográfica (R1)
**Duração**: 6-10 horas (varia muito com R0 vs R1)


#### R0 (Referências Travadas): 2-3 horas
| Atividade | Tempo |
|-----------|-------|
| Verificar referências aprovadas | 1h |
| Mapear citações necessárias | 1h |
| Documentar lacunas | 30min |

#### R1 (Referências Expandidas): 8-12 horas
| Atividade | Tempo | Ferramentas |
|-----------|-------|-------------|
| **Categoria 1 - Fundacionais** | 1-2h | Google Scholar, livros-texto |
| **Categoria 2 - Estado da Arte** | 2-3h | arXiv, Web of Science (2022-2024) |
| **Categoria 3 - Metodológicas** | 1-2h | Papers originais de métodos |
| **Categoria 4 - Estatísticas** | 1h | Livros de estatística |
| **Categoria 5 - Frameworks** | 30min | Docs oficiais (PennyLane, Qiskit) |
| **Categoria 6 - Críticas** | 1-2h | Buscar "limitations", "fails" |
| **Categoria 7 - Aplicações** | 1-2h | Casos de uso similares |
| Buscar DOIs para todas as refs | 1h | CrossRef, DOI.org |
| Formatar referências (ABNT/IEEE) | 1h | Zotero, Mendeley |
| Síntese de literatura (não lista!) | 2-3h | Análise crítica |
| Quality Gate F3 | 30min | Verificar DOIs, contrapontos |

**Deliverables**:
- ✅ `referencias_compiladas.md` (35-60 refs)
- ✅ `sintese_literatura.md`
- ✅ `taxonomia_estado_da_arte.md`


**⚠️ Nota**: R1 requer acesso a bases de dados institucionais. Se não disponível, use R0.


---


### Fase 4: Redação do Manuscrito (IMRAD)
**Duração**: 20-30 horas (a mais trabalhosa!)


| Seção | Tempo | Páginas (est.) | Dificuldade |
|-------|-------|----------------|-------------|
| **4.1 Resumo/Abstract** | 2-3h | 0.5 | ⭐⭐⭐ |
| **4.2 Introdução (CARS)** | 3-4h | 2-3 | ⭐⭐⭐⭐ |
| **4.3 Revisão de Literatura** | 3-4h | 3-4 | ⭐⭐⭐⭐ |
| **4.4 Metodologia** | 5-7h | 4-6 | ⭐⭐⭐⭐⭐ |
| - Formulação matemática | 2h | 1-2 | - |
| - Algorithm 1 (LaTeX) | 1h | 0.5 | - |
| - Tabela Código→Método | 2h | 1 | - |
| - Descrição de datasets | 1h | 0.5 | - |
| **4.5 Resultados** | 4-5h | 3-4 | ⭐⭐⭐⭐ |
| **4.6 Discussão** | 3-4h | 2-3 | ⭐⭐⭐⭐ |
| - Interpretação | 1h | 1 | - |
| - Threats to Validity | 1h | 0.5 | - |
| - Comparação com estado da arte | 1h | 0.5 | - |
| **4.7 Conclusão** | 1-2h | 1 | ⭐⭐ |
| **4.8 Seções Editoriais** | 1h | 1 | ⭐ |
| **4.9 Referências (formatação)** | 1-2h | 2-3 | ⭐⭐ |
| Quality Gate F-Redação | 1h | - | ⭐⭐ |

**Total de páginas (estimado)**: 20-30 páginas


**Deliverables**:
- ✅ 9 arquivos `.md` (IMRAD completo)


**⚠️ Dica**: Metodologia é a seção mais técnica. Reserve tempo extra para garantir precisão.


---


### Fase 5: Material Suplementar
**Duração**: 8-12 horas


| Atividade | Tempo | Complexidade |
|-----------|-------|--------------|
| **Tabela S1** (todas as configs) | 2-3h | ⭐⭐⭐⭐ |
| - Gerar CSV a partir de logs | 1h | - |
| - Validar total de linhas | 30min | - |
| - Adicionar metadados | 30min | - |
| **Tabela S2** (comparação estado da arte) | 1-2h | ⭐⭐⭐ |
| **Tabela S3** (hiperparâmetros) | 1h | ⭐⭐ |
| **Tabela S4** (testes post-hoc) | 1-2h | ⭐⭐⭐ |
| **Tabela S5** (configurações de hardware) | 30min | ⭐ |
| **Figuras S1-S8** (descrições) | 3-4h | ⭐⭐⭐ |
| Notas metodológicas adicionais | 1-2h | ⭐⭐ |
| Quality Gate F5 | 30min | ⭐⭐ |

**Deliverables**:
- ✅ `tabelas_suplementares.md`
- ✅ `tabela_s1_configuracoes.csv` (ex: 3.360 linhas)
- ✅ `figuras_suplementares.md`
- ✅ `notas_metodologicas_adicionais.md`


**⚠️ Nota**: S1 é crítica - deve ter exatamente N linhas (onde N = total de configs calculado na Fase 1).


---


### Fase 6: Consolidação e Verificação
**Duração**: 6-8 horas


| Atividade | Tempo | Criticidade |
|-----------|-------|-------------|
| Geração de relatório de consistência | 2-3h | 🔴 Alta |
| - Verificar números vs código | 1h | - |
| - Calcular % de conivência | 30min | - |
| - Identificar discrepâncias | 1h | - |
| **Tabela de Rastreabilidade Completa** | 3-4h | 🔴 Alta |
| - Preencher 50+ entradas | 2h | - |
| - Verificar origens | 1h | - |
| - Testar links arquivo:linha | 30min | - |
| Consolidar documento final | 1h | 🟡 Média |
| Gerar sumário executivo | 30min | 🟢 Baixa |
| Quality Gate Final | 1h | 🔴 Alta |
| - Consistência ≥95%? | 30min | - |
| - Checklist 100 pontos | 30min | - |

**Deliverables**:
- ✅ `relatorio_consistencia.md`
- ✅ `rastreabilidade_completa.md` (50+ entradas)
- ✅ `artigo_abnt_final.md` ou `manuscrito_internacional_final.tex`
- ✅ `sumario_executivo.md`


**🎯 Meta**: Consistência ≥95%, idealmente 100%


---


## 📊 Resumo por Fase

| Fase | Duração | % do Tempo Total | Dificuldade |
|------|---------|------------------|-------------|
| 1. Auditoria | 8-12h | 15% | ⭐⭐⭐ |
| 2. Enquadramento | 4-6h | 8% | ⭐⭐ |
| 3. Bibliografia (R1) | 6-10h | 14% | ⭐⭐⭐ |
| 4. Redação IMRAD | 20-30h | 42% | ⭐⭐⭐⭐⭐ |
| 5. Suplementar | 8-12h | 15% | ⭐⭐⭐⭐ |
| 6. Consolidação | 6-8h | 11% | ⭐⭐⭐⭐ |
| **TOTAL** | **52-78h** | **100%** | - |

---


## 🗓️ Calendários Exemplo

### Cenário 1: Dedicação Integral (8h/dia)

| Dia | Fase | Atividades |
|-----|------|-----------|
| **Dia 1** | Fase 1 | Auditoria completa (8h) |
| **Dia 2** | Fase 2 + Fase 3 início | Enquadramento (4h) + Bibliografia início (4h) |
| **Dia 3** | Fase 3 | Bibliografia R1 completa (8h) |
| **Dia 4** | Fase 4 | Redação: Resumo, Intro, Revisão (8h) |
| **Dia 5** | Fase 4 | Redação: Metodologia (8h) |
| **Dia 6** | Fase 4 | Redação: Resultados, Discussão (8h) |
| **Dia 7** | Fase 4 + Fase 5 | Conclusão, Editoriais (2h) + Suplementar início (6h) |
| **Dia 8** | Fase 5 | Suplementar completo (8h) |
| **Dia 9** | Fase 6 | Consolidação e verificação (8h) |
| **Dia 10** | Revisão Final | Ajustes finais, checklist 100pts (8h) |

**Total**: 10 dias úteis (2 semanas de calendário)


---


### Cenário 2: Dedicação Parcial (4h/dia)

| Semana | Fases | Observações |
|--------|-------|-------------|
| **Semana 1** | Fases 1-3 | Auditoria, enquadramento, bibliografia |
| **Semana 2** | Fase 4 (parte 1) | Resumo, intro, revisão, métodos |
| **Semana 3** | Fase 4 (parte 2) | Resultados, discussão, conclusão |
| **Semana 4** | Fases 5-6 | Suplementar, consolidação, revisão |

**Total**: 4 semanas (1 mês)


---


### Cenário 3: Dedicação Esporádica (2h/dia, fins de semana)

| Período | Fases | Observações |
|---------|-------|-------------|
| **Semanas 1-2** | Fase 1 | Auditoria em partes |
| **Semanas 3-4** | Fases 2-3 | Enquadramento e bibliografia |
| **Semanas 5-8** | Fase 4 | Redação IMRAD (a mais longa) |
| **Semanas 9-10** | Fases 5-6 | Suplementar e consolidação |

**Total**: 10 semanas (2.5 meses)


---


## ⚡ Fatores que Aceleram

1. **Código bem documentado**: -20% tempo (Fase 1)
2. **Seeds fixas e logs existentes**: -30% tempo (Fase 1)
3. **Experiência com LaTeX/ABNT**: -15% tempo (Fase 4)
4. **Acesso a bases de dados**: R1 não é bloqueado (Fase 3)
5. **Uso de Zotero/Mendeley**: -25% tempo (Fase 3)
6. **Templates prontos**: -10% tempo (Fase 4)
7. **Rascunhos parciais existentes**: -20% tempo (Fases 4-5)


**Com todos os fatores**: Pode reduzir para ~35-50h (4-6 dias)


---


## 🐌 Fatores que Retardam

1. **Código sem documentação**: +30% tempo
2. **Sem logs/resultados prontos**: +50% tempo (precisa executar experimentos)
3. **Primeira vez gerando artigo**: +25% tempo (curva de aprendizado)
4. **Restrições de referências (R0 com muitas lacunas)**: +15% tempo
5. **Mudanças de escopo durante processo**: +40% tempo
6. **Metodologia complexa/inovadora**: +20% tempo
7. **Múltiplos coautores com opiniões divergentes**: +50% tempo


**Com vários fatores negativos**: Pode estender para 100-120h (2-3 semanas)


---


## 🎯 Dicas para Otimização de Tempo

### Antes de Começar (Preparação)

1. ✅ **Execute o código completamente** e salve todos os logs
2. ✅ **Fixe seeds** no código antes de iniciar
3. ✅ **Organize referências** (Zotero, Mendeley) antecipadamente
4. ✅ **Defina o periódico-alvo** (MODE_A ou MODE_B)
5. ✅ **Configure access** a Web of Science/Scopus (se R1)


### Durante o Processo

6. ⏰ **Reserve blocos contínuos** de 4h+ (não interrompa)
7. 🎯 **Comece pela Fase 1** (base para tudo)
8. 📝 **Não seja perfeccionista** na primeira passada
9. 🔁 **Use placeholders** `[A COMPLETAR]` e volte depois
10. 🤝 **Delegue quando possível** (coautores podem ajudar na bibliografia)


### Após Completar

11. 🔍 **Quality Gates são críticos** - não pule
12. 📊 **Use o checklist de 100 pontos** para priorizar ajustes
13. 🚀 **Submeta logo** - não espere "perfeição absoluta"


---


## 📈 Benchmarks Reais

**Projeto Beneficial Quantum Noise** (este repositório):
- Fase 1: 6h (código bem estruturado, analyzer automático)
- Fase 2: 3h (lacuna clara)
- Fase 3: 8h (R1, 45 referências)
- Fase 4: 24h (metodologia complexa, 28 páginas)
- Fase 5: 10h (3.360 configurações em S1)
- Fase 6: 7h (consistência de 96%)
- **Total real**: 58h (~7 dias úteis)


---


## 🏁 Checklist de Prontidão

Antes de iniciar, verifique:

- [ ] Código executa sem erros
- [ ] Logs de execução existem ou podem ser gerados
- [ ] Seeds são fixas ou podem ser fixadas
- [ ] Dependências têm versões especificadas
- [ ] Acesso a bases de dados (se R1)
- [ ] Tempo disponível (mínimo 6 dias úteis)
- [ ] Templates baixados
- [ ] `enhanced_code_analyzer.py` testado


**Se todos ✅**: Você está pronto para começar! 🚀


---


**Versão**: 1.0  
**Última Atualização**: 26/12/2025  
**Autor**: Sistema de Geração de Artigos Qualis A1

