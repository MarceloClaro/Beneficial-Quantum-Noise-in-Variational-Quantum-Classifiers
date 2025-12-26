# 📄 Gerador de Artigo Científico Completo - MODO B + R1

## 🎯 Visão Geral

Este módulo implementa um **sistema completo de geração de artigos científicos** com rastreabilidade total entre código/dados e texto, seguindo:

- **MODO B**: Texto em PORTUGUÊS + normas ABNT (NBR 10520/6023)
- **R1**: Política de referências expandidas (busca permitida com DOI)

O sistema integra com `consultor_metodologico.py` para análise de qualidade metodológica.

---

## 🏗️ Arquitetura do Sistema

### 6 Fases com Quality Gates

```
┌─────────────────────────────────────────────────────────────┐
│ FASE 1: Auditoria Técnica                                  │
│ ├─ analise_codigo_inicial.md                               │
│ ├─ tabela_componentes.md                                   │
│ └─ mapa_execucao.md                                        │
│ Quality Gate F1: Rastreabilidade completa                  │
└─────────────────────────────────────────────────────────────┘
          ↓
┌─────────────────────────────────────────────────────────────┐
│ FASE 2: Enquadramento Científico                           │
│ ├─ linha_de_pesquisa.md                                    │
│ └─ diagrama_linha_pesquisa.md                              │
│ Quality Gate F2: Lacuna operacionalizável                  │
└─────────────────────────────────────────────────────────────┘
          ↓
┌─────────────────────────────────────────────────────────────┐
│ FASE 3: Curadoria Bibliográfica (R1)                       │
│ ├─ referencias_compiladas.md (35-60 refs)                  │
│ └─ sintese_literatura.md                                   │
│ Quality Gate F3: DOIs + contrapontos                       │
└─────────────────────────────────────────────────────────────┘
          ↓
┌─────────────────────────────────────────────────────────────┐
│ FASE 4: Redação IMRAD (PORTUGUÊS)                          │
│ ├─ resumo_abstract.md                                      │
│ ├─ introducao_completa.md                                  │
│ ├─ revisao_literatura_completa.md                          │
│ ├─ metodologia_completa.md                                 │
│ ├─ resultados_completo.md                                  │
│ ├─ discussao_completa.md                                   │
│ ├─ conclusao_completa.md                                   │
│ ├─ secoes_editoriais.md                                    │
│ └─ agradecimentos_referencias.md                           │
│ Quality Gate F-Redação: Números rastreáveis + ABNT         │
└─────────────────────────────────────────────────────────────┘
          ↓
┌─────────────────────────────────────────────────────────────┐
│ FASE 5: Material Suplementar                               │
│ ├─ tabelas_suplementares.md (S1-S5)                        │
│ ├─ tabela_s1_configuracoes.csv                             │
│ ├─ figuras_suplementares.md (S1-S8)                        │
│ └─ notas_metodologicas_adicionais.md                       │
│ Quality Gate F5: Consistência experimental                 │
└─────────────────────────────────────────────────────────────┘
          ↓
┌─────────────────────────────────────────────────────────────┐
│ FASE 6: Consolidação e Verificação                         │
│ ├─ relatorio_consistencia.md                               │
│ ├─ rastreabilidade_completa.md                             │
│ ├─ artigo_abnt_final.md                                    │
│ └─ sumario_executivo.md                                    │
│ Quality Gate Final: Consistência ≥95% + Reprodutibilidade  │
└─────────────────────────────────────────────────────────────┘
```

---

## 🚀 Uso

### Comando Básico

```bash
python gerador_artigo_completo.py --repositorio . --output artigo_gerado
```

### Parâmetros

| Parâmetro | Descrição | Padrão |
|-----------|-----------|--------|
| `--repositorio` | Caminho para o repositório do código | `.` (atual) |
| `--output` | Pasta de saída para arquivos gerados | `artigo_gerado` |
| `--periodico-primario` | Periódico-alvo primário | `[especificar]` |

### Exemplo Completo

```bash
# Gerar artigo do projeto atual
python gerador_artigo_completo.py \
  --repositorio /path/to/meu/projeto \
  --output artigo_meu_projeto \
  --periodico-primario "Nature Communications"
```

---

## 📂 Estrutura de Saída

Após execução, a pasta de saída terá:

```
artigo_gerado/
├── fase1_auditoria/
│   ├── analise_codigo_inicial.md
│   ├── tabela_componentes.md
│   └── mapa_execucao.md
├── fase2_enquadramento/
│   ├── linha_de_pesquisa.md
│   └── diagrama_linha_pesquisa.md
├── fase3_literatura/
│   ├── referencias_compiladas.md
│   └── sintese_literatura.md
├── fase4_redacao/
│   ├── resumo_abstract.md
│   ├── introducao_completa.md
│   ├── revisao_literatura_completa.md
│   ├── metodologia_completa.md
│   ├── resultados_completo.md
│   ├── discussao_completa.md
│   ├── conclusao_completa.md
│   ├── secoes_editoriais.md
│   └── agradecimentos_referencias.md
├── fase5_suplementar/
│   ├── tabelas_suplementares.md
│   ├── tabela_s1_configuracoes.csv
│   ├── figuras_suplementares.md
│   └── notas_metodologicas_adicionais.md
└── fase6_consolidacao/
    ├── relatorio_consistencia.md
    ├── rastreabilidade_completa.md
    ├── artigo_abnt_final.md
    └── sumario_executivo.md
```

**Total**: 24 arquivos gerados automaticamente

---

## 📋 Quality Gates

### Quality Gate F1 (Fase 1)
- ✅ Todas as listas têm origem clara (arquivo/função/linha)
- ✅ Total de configurações calculado e verificável
- ✅ Nenhuma afirmação sem suporte (marcadas [INFORMAÇÃO AUSENTE])

### Quality Gate F2 (Fase 2)
- ✅ Pergunta/objetivos explicitados e alinhados ao código
- ✅ Lacuna é falsificável/operacionalizável

### Quality Gate F3 (Fase 3)
- ✅ Referências têm DOI quando disponível
- ✅ Toda técnica central tem referência
- ✅ Contrapontos (críticos) incluídos (R1)

### Quality Gate F-Redação (Fase 4)
- ✅ Texto não contém números "sem lastro"
- ✅ Seções seguem tom e estrutura MODO B
- ✅ Referências adicionadas rastreáveis (R1)

### Quality Gate F5 (Fase 5)
- ✅ S1 bate com cálculo de configurações
- ✅ Figuras/tabelas têm fonte e script indicados

### Quality Gate Final (Fase 6)
- ✅ Consistência código-texto ≥ 95%
- ✅ Citação↔referência 100% consistente
- ✅ Reprodutibilidade garantida (instruções + ambiente)
- ✅ Limitações explicitadas

---

## 🔍 Rastreabilidade

Cada número ou afirmação no texto tem:

1. **Origem no código**: arquivo, função, linha
2. **Evidência**: log, tabela, figura
3. **Referência**: citação ABNT quando aplicável

Exemplo:
```markdown
A acurácia média foi de 65.83% (IC95%: 63.1-68.5)
[ORIGEM: framework_qiskit.py, linha 245, função calcular_metricas()]
[EVIDÊNCIA: logs/experimento_2025-12-23.log, linha 1847]
[REFERÊNCIA: DU et al., 2021]
```

---

## 📚 Integração com Consultor Metodológico

O gerador pode usar o `consultor_metodologico.py` para:

- **Tarefa A**: Justificativa metodológica (integrada na Metodologia)
- **Tarefa E**: Verificação de elementos da Introdução
- **Tarefa G**: Tabelas comparativas de conceitos (Revisão de Literatura)

Uso integrado:
```bash
# 1. Gerar artigo completo
python gerador_artigo_completo.py --output artigo_gerado

# 2. Revisar metodologia com consultor
python consultor_metodologico.py \
  --insumos artigo_gerado/fase4_redacao/metodologia_completa.md \
  --tarefas A
```

---

## 🎓 Padrões ABNT Implementados

### NBR 10520 (Citações)

**Autor-data:**
```
Du et al. (2021) demonstraram que...
```

**Parentética:**
```
Ruído benéfico melhora VQCs (DU et al., 2021).
```

**Citação direta (com página obrigatória):**
```
"Quantum noise acts as a regularizer" (DU et al., 2021, p. 5).
```

### NBR 6023 (Referências)

```
DU, Y. et al. Quantum noise protects quantum classifiers against 
adversaries and reduces overfitting. Physical Review Applied, v. 15, 
n. 3, p. 034026, 2021. DOI: 10.1103/PhysRevApplied.15.034026.
```

---

## 🛠️ Desenvolvimento Futuro

### Próximas Implementações

1. **Extração automática de código**:
   - Parser AST para identificar modelos/funções
   - Análise de logs para extrair métricas
   - Cálculo automático de configurações experimentais

2. **Busca automática de referências (R1)**:
   - Integração com Crossref API
   - Busca em arXiv, PubMed, Scopus
   - Validação de DOIs

3. **Geração de conteúdo IMRAD**:
   - Templates específicos por área
   - Análise estatística automatizada (ANOVA, post-hoc)
   - Geração de tabelas e figuras

4. **Checklist de rigor**:
   - Validação estatística automática
   - Verificação de pressupostos
   - Análise de poder estatístico

---

## 📊 Exemplo de Saída

### Arquivo: `analise_codigo_inicial.md`

```markdown
# Análise Inicial do Código e Dados

**Data:** 26/12/2025
**Modo:** B (ABNT)
**Política de Referências:** R1 (Expandidas)

## 1. Estrutura Técnica

### 1.1 Estrutura do Projeto

**Pastas principais:**
- artigo_cientifico
- docs
- examples
- figuras
- notebooks
- tests
- tools

**Módulos Python:**
- consultor_metodologico.py
- framework_investigativo_completo.py
- framework_qiskit.py
- gerador_artigo_completo.py

### 1.2 Dependências e Versões

| Biblioteca | Versão |
|-----------|--------|
| numpy | 1.24.3 |
| scipy | 1.11.1 |
| pandas | 2.0.3 |
| matplotlib | 3.7.2 |
| pennylane | 0.38.0 |
| qiskit | 1.0.0 |

[... continua ...]
```

---

## 🤝 Contribuindo

Para melhorar o gerador:

1. Fork o repositório
2. Crie branch (`git checkout -b feature/melhoria`)
3. Commit (`git commit -am 'Adiciona funcionalidade X'`)
4. Push (`git push origin feature/melhoria`)
5. Abra Pull Request

---

## 📄 Licença

MIT License - veja LICENSE para detalhes.

---

## 📞 Suporte

- 🐛 [Reportar bug](https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers/issues)
- 💡 [Sugerir melhoria](https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers/issues)
- 📖 [Documentação completa](GERADOR_ARTIGO_README.md)

---

**Versão:** 1.0.0  
**Última atualização:** 26/12/2025  
**Status:** ✅ Funcional (estrutura base implementada)
