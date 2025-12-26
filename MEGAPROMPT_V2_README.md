# MegaPrompt v2.0 — Geração de Artigos Científicos Qualis A1

## 🎯 Visão Geral

Este diretório implementa o **MegaPrompt v2.0**, um framework completo e rigoroso para geração de artigos científicos de alto impacto (Qualis A1) com 100% de rastreabilidade entre código, dados e texto.

**Periódicos-Alvo:**
- Nature family journals
- Science
- Quantum (Verein zur Förderung des Open Access Publizierens)
- Physical Review series (A, Research, Letters)
- npj Quantum Information

## 📖 Estrutura do Framework

### PARTE I: Configuração e Planejamento

#### Configuração
- **`config.json`**: Configuração principal do projeto
  - Modo de saída (LaTeX/ABNT)
  - Política de referências (R0/R1)
  - Perfil editorial
  - Periódicos-alvo
  - Caminhos de entrada

#### Documentação Base
- **`GLOSSARIO.md`**: Termos técnicos e definições
- **`FAQ_TROUBLESHOOTING.md`**: Perguntas frequentes e resolução de problemas
- **`CRONOGRAMA_ESTIMADO.md`**: Estimativa de tempo por fase

### PARTE II: Execução (6 Fases)

#### Fase 1: Auditoria Técnica + Templates
📂 Localização: `artigo_cientifico/fase1_analise/`

**Arquivos:**
- `analise_codigo_inicial.md`: Análise técnica completa
- `tabela_componentes.md`: Resumo executivo técnico
- `mapa_execucao.md`: Passo a passo reprodutível
- `manifesto_execucao.json`: Comandos, seeds, configs

**Quality Gate F1:**
- Cada item tem origem verificável
- Total de configurações calculado
- Ausências explicitadas

#### Fase 2: Bibliografia + Fluxograma R0/R1
📂 Localização: `artigo_cientifico/fase2_bibliografia/`

**Arquivos:**
- `FLUXOGRAMA_R0_R1.md`: Diagrama de decisão para referências
- `referencias_compiladas.md`: Referências com DOI e justificativa
- `sintese_literatura.md`: Síntese crítica
- `taxonomia_estado_da_arte.md`: Tabela comparativa

**Quality Gate F2:**
- Cada técnica central tem referência ou [LACUNA]
- Síntese contém contraste/avaliação

#### Fase 3: Projeto do Artigo + Formal Problem Statement
📂 Localização: `artigo_cientifico/fase3_estrutura/`

**Arquivos:**
- `problema_formal.md`: Formulação matemática do problema
- `titulos_palavras_chave.md`: Opções de título + keywords
- `hipoteses_objetivos.md`: Hipóteses testáveis + objetivos SMART

**Quality Gate F3:**
- Problema formal compatível com código
- Cada hipótese tem teste/métrica correspondente

#### Fase 4: Redação + Algorithm 1 + Tabela Código→Método
📂 Localização: `artigo_cientifico/fase4_secoes/`

**Arquivos:**
- `resumo_abstract.md`: Abstract IMRAD
- `introducao_completa.md`: Introdução CARS
- `revisao_literatura_completa.md`: Related Work
- `metodologia_completa.md`: Methods com Algorithm 1 e tabela
- `resultados_completo.md`: Results com evidências
- `discussao_completa.md`: Discussion crítica
- `conclusao_completa.md`: Conclusion + Future Work
- `secoes_editoriais.md`: Editorial sections
- `agradecimentos_referencias.md`: Acknowledgments + References

**Templates Especiais:**
- `templates/algorithm_latex_template.tex`: Template LaTeX para Algorithm 1
- `templates/tabela_codigo_metodo_template.md`: Mapeamento código→método

**Quality Gate F4:**
- Sem números sem lastro
- R0 respeitado (se aplicável)
- Methods completo (notação + equações + algoritmo + mapa)

#### Fase 5: Material Suplementar + Especificações Detalhadas
📂 Localização: `artigo_cientifico/fase5_suplementar/`

**Arquivos:**
- `tabelas_suplementares.md`: Tabelas S1-S5
- `tabela_s1_configuracoes.csv`: CSV completo de configurações
- `figuras_suplementares.md`: Descrições detalhadas (8 figuras)
- `apendice_suplementar.md`: Derivações matemáticas

**Tabelas Obrigatórias:**
- **S1**: Todas as 2.688 configurações experimentais
- **S2**: Comparação quantitativa com estado da arte
- **S3**: Especificação completa de hiperparâmetros
- **S4**: Testes post-hoc com correção de Bonferroni
- **S5**: Análise de poder estatístico

**Quality Gate F5:**
- S1 confere com total calculado na F1
- Cada tabela/figura aponta script/config/log
- Material core não está apenas no suplemento

#### Fase 6: Consolidação + Tabela de Rastreabilidade + Checklist
📂 Localização: `artigo_cientifico/fase6_consolidacao/`

**Arquivos:**
- `rastreabilidade_completa.md`: Tabela de rastreabilidade
- `relatorio_consistencia.md`: Verificação de conivência
- `manuscrito_internacional_final.tex`: Artigo final LaTeX
- `artigo_abnt_final.md`: Artigo final ABNT (se MODE_B)
- `sumario_executivo.md`: Resumo executivo
- `checklist_auditoria_100pts.md`: Checklist de 0-100 pontos

**Quality Gate Final:**
- Consistência ≥ 95% (meta: 100%)
- Citação↔referência 100%
- Reprodutibilidade completa
- Ameaças à validade explicitadas

### PARTE III: Apoio e Referência

#### Ferramentas de Automação
📂 Localização: `tools/megaprompt_v2/`

**Scripts Python:**
- `generate_s1.py`: Gera Tabela S1 a partir dos logs
- `check_consistency.py`: Verifica conivência código-texto
- `build_paper.sh`: Consolida todas as seções
- `validate_traceability.py`: Valida rastreabilidade

#### Documentação de Apoio
- `EXEMPLOS_PRATICOS.md`: Exemplos do caso Beneficial Quantum Noise
- `FERRAMENTAS_WORKFLOWS.md`: Guia de ferramentas
- `CHECKLIST_AUDITORIA.md`: Checklist expandido (0-100 pontos)

## 🚀 Como Usar

### 1. Configuração Inicial

```bash
# Clone o repositório
git clone https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers.git
cd Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers

# Configure o config.json
nano config.json
```

### 2. Execução das Fases

```bash
# Fase 1: Auditoria Técnica
python tools/megaprompt_v2/fase1_audit.py

# Fase 2: Bibliografia
python tools/megaprompt_v2/fase2_bibliography.py

# Fase 3: Estrutura
python tools/megaprompt_v2/fase3_structure.py

# Fase 4: Redação
python tools/megaprompt_v2/fase4_writing.py

# Fase 5: Suplementar
python tools/megaprompt_v2/fase5_supplementary.py

# Fase 6: Consolidação
python tools/megaprompt_v2/fase6_consolidation.py
```

### 3. Geração do Artigo Final

```bash
# Gera o artigo consolidado
bash tools/megaprompt_v2/build_paper.sh

# Valida rastreabilidade
python tools/megaprompt_v2/validate_traceability.py

# Executa checklist de auditoria
python tools/megaprompt_v2/audit_checklist.py
```

## 📊 Regras Mandatórias de Integridade

1. **NÃO inventar detalhes**: Se algo não estiver em código/dados/logs, usar **[INFORMAÇÃO AUSENTE]**
2. **NÃO inventar números**: Todo valor quantitativo deve ter lastro verificável; caso contrário **[NÃO DISPONÍVEL]**
3. **Se R0**: É proibido alterar o conjunto de referências; quando faltar base, usar **[LACUNA DE CITAÇÃO]**
4. **Reprodutibilidade**: Reportar HW/SW, versões, seeds, configs, scripts e comandos
5. **Auditoria**: Cada seção exige rastreabilidade: **Seção → Evidência → Origem**

## 📈 Métricas de Qualidade

### Checklist de Auditoria (0-100 pontos)

**1. Reprodutibilidade (30 pts)**
- Ambiente documentado (10 pts)
- Seeds fixas e reportadas (10 pts)
- Pipeline executável (10 pts)

**2. Rastreabilidade (30 pts)**
- Tabela de rastreabilidade completa (15 pts)
- Mapa código→método completo (15 pts)

**3. Rigor Estatístico (20 pts)**
- Testes apropriados (5 pts)
- Correção para múltiplas comparações (5 pts)
- Intervalos de confiança (5 pts)
- Tamanhos de efeito (5 pts)

**4. Transparência (20 pts)**
- Código disponível publicamente (10 pts)
- Dados disponíveis publicamente (5 pts)
- Limitações e ameaças à validade discutidas (5 pts)

**Meta**: ≥ 90/100 pontos

## 🎓 Periódicos e Conformidade

### Formato por Periódico

- **MODE_A** (Inglês/LaTeX): Nature, Science, Quantum, Physical Review, npj QI
- **MODE_B** (Português/ABNT): Periódicos brasileiros Qualis A1

### Perfis Editoriais

- **PROFILE_PR_QUANTUM**: Técnico, matemático, focado em física
- **PROFILE_GENERAL**: Narrativo, acessível, amplo público

## 📚 Referências

Este framework implementa boas práticas de:
- **Nature Scientific Data**: "Code availability"
- **PLoS Computational Biology**: "Software and data availability"
- **Quantum Journal**: "Reproducibility statement"
- **ACM**: "Artifact Evaluation"
- **IEEE**: "Open Science"

### Citações Chave

1. Peng, R. D. (2011). "Reproducible research in computational science." *Science*, 334(6060), 1226-1227.
2. Sandve et al. (2013). "Ten simple rules for reproducible computational research." *PLoS Computational Biology*, 9(10), e1003285.
3. Stodden, V. (2014). "The scientific method in practice." *MIT Sloan Research Paper* No. 4773-10.

## ⏱️ Cronograma Estimado

- **Fase 1 (Auditoria)**: 8-12 horas
- **Fase 2 (Bibliografia)**: 6-10 horas
- **Fase 3 (Projeto)**: 4-6 horas
- **Fase 4 (Redação)**: 20-30 horas
- **Fase 5 (Suplementar)**: 8-12 horas
- **Fase 6 (Consolidação)**: 6-8 horas
- **Total**: 52-78 horas (6-10 dias úteis)

## ✅ Status do Projeto

**Versão Atual**: 2.0  
**Framework Base**: v8.0-QAI  
**Status**: ✅ Implementado  
**Última Atualização**: 2025-12-26

## 📧 Suporte

Para questões ou sugestões:
- Abra uma issue no GitHub
- Consulte `FAQ_TROUBLESHOOTING.md`
- Veja exemplos em `EXEMPLOS_PRATICOS.md`

---

**Nota**: Este framework garante a mais alta qualidade científica e conformidade com padrões internacionais de publicação. Use-o com confiança para submissões a periódicos de alto impacto.
