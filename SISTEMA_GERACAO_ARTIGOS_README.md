# Sistema de Geração de Artigos Científicos Qualis A1 com Rastreabilidade Total

## 🎯 Visão Geral

Este sistema implementa um framework completo para geração de artigos científicos de alto impacto (Qualis A1, Nature, Science, Physical Review, Quantum) com 100% de rastreabilidade entre código, dados e texto.

### Características Principais

✅ **6 Fases com Quality Gates**: Desde auditoria técnica até consolidação final  
✅ **Rastreabilidade Total**: Cada afirmação rastreável até o código-fonte  
✅ **Análise Automatizada**: Extração inteligente de componentes do código  
✅ **Templates Profissionais**: LaTeX, Markdown, ABNT, IEEE  
✅ **Checklist de 100 Pontos**: Sistema objetivo de avaliação de qualidade  
✅ **R0/R1 Policies**: Controle rigoroso de referências bibliográficas  

---

## 📂 Estrutura do Sistema

```
/
├── config.json                          # Configuração principal
├── enhanced_code_analyzer.py            # Analisador de código automático
├── gerador_artigo_completo.py          # Gerador principal (6 fases)
│
├── templates/                           # Templates profissionais
│   ├── problema_formal_template.md     # Formulação matemática do problema
│   ├── algorithm_latex_template.tex    # Algoritmos 1, 2, 3 em LaTeX
│   ├── tabela_codigo_metodo_template.md    # Rastreabilidade Código→Método
│   └── rastreabilidade_completa_template.md  # Tabela Seção→Evidência→Origem
│
├── GLOSSARIO.md                         # 90+ termos técnicos definidos
├── FAQ_TROUBLESHOOTING.md               # Perguntas frequentes e soluções
├── FLUXOGRAMA_R0_R1.md                  # Políticas de referências
├── CHECKLIST_AUDITORIA_100PTS.md        # Sistema de avaliação objetivo
│
└── artigo_gerado/                       # Output (criado pelo gerador)
    ├── fase1_auditoria/
    ├── fase2_enquadramento/
    ├── fase3_literatura/
    ├── fase4_redacao/
    ├── fase5_suplementar/
    └── fase6_consolidacao/
```

---

## 🚀 Início Rápido

### 1. Instalação

```bash
# Clone o repositório
git clone https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers.git
cd Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers

# Instale dependências (se necessário)
pip install -r requirements.txt
```

### 2. Configuração

Edite `config.json` com suas preferências:

```json
{
  "output_mode": "MODE_A",              # MODE_A (EN/LaTeX) ou MODE_B (PT/ABNT)
  "reference_policy": "R1",             # R0 (travadas) ou R1 (expandidas)
  "editorial_profile": "PROFILE_PR_QUANTUM",
  "target_journals": {
    "primary": "Quantum",
    "secondary": ["Physical Review A", "Nature Communications"]
  },
  "inputs": {
    "code_path": ".",
    "data_path": "[Gerado pelo código]",
    "artifacts_path": "[Gerado pelo código]"
  }
}
```

### 3. Análise Automática do Código

```bash
# Execute o analisador de código
python enhanced_code_analyzer.py .

# Output: code_analysis_report.json com todos os componentes extraídos
```

**O que é extraído automaticamente:**
- ✅ Ansätze/arquiteturas quânticas
- ✅ Modelos de ruído (Depolarizing, Amplitude Damping, etc.)
- ✅ Datasets (Moons, Circles, Iris, etc.)
- ✅ Métricas (Accuracy, F1-Score, etc.)
- ✅ Schedules de ruído (Constant, Linear, Cosine)
- ✅ Seeds aleatórias
- ✅ **Cálculo automático de configurações** (ex: 3.360 configs)

### 4. Geração do Artigo

```bash
# Gerar artigo completo (6 fases)
python gerador_artigo_completo.py \
    --repositorio . \
    --output artigo_quantum_noise \
    --periodico-primario "Quantum"
```

**Tempo estimado:** 10-30 minutos (execução automatizada)

---

## 📋 As 6 Fases de Geração

### Fase 1: Auditoria Técnica ✅
**Objetivo:** Inventariar todos os componentes do projeto  
**Outputs:**
- `analise_codigo_inicial.md` - Estrutura completa do projeto
- `tabela_componentes.md` - Resumo executivo técnico
- `mapa_execucao.md` - Pipeline reprodutível
- `manifesto_execucao.json` - Comandos, seeds, configs

**Quality Gate F1:** Rastreabilidade completa, sem números inventados

---

### Fase 2: Enquadramento Científico 🔬
**Objetivo:** Posicionar o trabalho na linha de pesquisa  
**Outputs:**
- `linha_de_pesquisa.md` - Área, subárea, problema central
- `diagrama_linha_pesquisa.md` - Fluxograma com Mermaid

**Quality Gate F2:** Lacuna operacionalizável e testável

---

### Fase 3: Curadoria Bibliográfica (R1) 📚
**Objetivo:** Compilar referências em 7 categorias  
**Outputs:**
- `referencias_compiladas.md` - 35-60 refs com DOI
- `sintese_literatura.md` - Análise crítica (não lista)
- `taxonomia_estado_da_arte.md` - Tabela comparativa

**7 Categorias de Busca:**
1. Fundacionais (clássicos, livros-texto)
2. Estado da Arte (últimos 2-3 anos)
3. Metodológicas (algoritmos específicos)
4. Estatísticas (testes, correções)
5. Frameworks (PennyLane, Qiskit)
6. **Críticas/Contrapontos** (limitações, falhas)
7. Aplicações (casos de uso)

**Quality Gate F3:** DOIs presentes, contrapontos incluídos

---

### Fase 4: Redação IMRAD ✍️
**Objetivo:** Redigir todas as seções do artigo  
**Outputs:**
- `resumo_abstract.md`
- `introducao_completa.md` (modelo CARS)
- `revisao_literatura_completa.md`
- `metodologia_completa.md` (com Algorithm 1 e tabela Código→Método)
- `resultados_completo.md`
- `discussao_completa.md` (com Threats to Validity)
- `conclusao_completa.md`
- `secoes_editoriais.md`
- `agradecimentos_referencias.md` (ABNT ou IEEE)

**Quality Gate F4:** Sem números sem lastro, tom consistente (MODE_A/B)

---

### Fase 5: Material Suplementar 📊
**Objetivo:** Gerar tabelas e figuras suplementares  
**Outputs:**
- `tabelas_suplementares.md` (S1-S5)
- `tabela_s1_configuracoes.csv` (todas as N configurações)
- `figuras_suplementares.md` (descrições S1-S8)
- `notas_metodologicas_adicionais.md`

**Quality Gate F5:** S1 confere com cálculo da Fase 1

---

### Fase 6: Consolidação e Auditoria ✅
**Objetivo:** Verificar consistência e rastreabilidade  
**Outputs:**
- `relatorio_consistencia.md` (% de conivência código-texto)
- `rastreabilidade_completa.md` (tabela com 50+ entradas)
- `artigo_abnt_final.md` ou `manuscrito_internacional_final.tex`
- `sumario_executivo.md`

**Quality Gate Final:** Consistência ≥95%, reprodutibilidade completa

---

## 🛠️ Ferramentas e Scripts

### 1. Analisador de Código (`enhanced_code_analyzer.py`)

```bash
python enhanced_code_analyzer.py /path/to/repo

# Gera: code_analysis_report.json
```

**Funcionalidades:**
- Extração automática de componentes
- Cálculo de configurações experimentais
- Identificação de seeds, métricas, datasets

---

### 2. Gerador Principal (`gerador_artigo_completo.py`)

```bash
# Uso básico
python gerador_artigo_completo.py --repositorio . --output artigo

# Opções avançadas
python gerador_artigo_completo.py \
    --repositorio /path/to/project \
    --output artigo_customizado \
    --periodico-primario "Nature Communications" \
    --mode B \                    # PT/ABNT
    --policy R0                   # Refs travadas
```

**Parâmetros:**
- `--repositorio`: Caminho do projeto a analisar
- `--output`: Pasta de saída
- `--periodico-primario`: Periódico-alvo
- `--mode`: A (Internacional/LaTeX) ou B (ABNT/Markdown)
- `--policy`: R0 (travadas) ou R1 (expandidas)

---

### 3. Verificador de Rastreabilidade

```bash
# Verificar consistência do artigo
python tools/verificar_rastreabilidade.py \
    artigo_gerado/fase6_consolidacao/rastreabilidade_completa.md \
    code_analysis_report.json
```

---

## 📖 Documentação Completa

### Para Iniciantes

1. 📘 **[GLOSSARIO.md](GLOSSARIO.md)**: Entenda todos os termos técnicos
2. ❓ **[FAQ_TROUBLESHOOTING.md](FAQ_TROUBLESHOOTING.md)**: Problemas comuns e soluções

### Para Usuários Avançados

3. 🔀 **[FLUXOGRAMA_R0_R1.md](FLUXOGRAMA_R0_R1.md)**: Políticas de referências
4. ✅ **[CHECKLIST_AUDITORIA_100PTS.md](CHECKLIST_AUDITORIA_100PTS.md)**: Sistema de avaliação

### Templates

5. 📐 **[templates/problema_formal_template.md](templates/problema_formal_template.md)**: Formulação matemática
6. 🔢 **[templates/algorithm_latex_template.tex](templates/algorithm_latex_template.tex)**: Algoritmos em LaTeX
7. 🗺️ **[templates/tabela_codigo_metodo_template.md](templates/tabela_codigo_metodo_template.md)**: Rastreabilidade
8. 📊 **[templates/rastreabilidade_completa_template.md](templates/rastreabilidade_completa_template.md)**: Auditoria completa

---

## 📊 Checklist de 100 Pontos

### Categorias

| Categoria | Peso | Descrição |
|-----------|------|-----------|
| **1. Reprodutibilidade** | 30 pts | Ambiente, seeds, pipeline executável |
| **2. Rastreabilidade** | 30 pts | Tabelas completas, mapa código→método |
| **3. Rigor Estatístico** | 20 pts | Testes, correções, IC, effect sizes |
| **4. Transparência** | 20 pts | Código/dados públicos, limitações discutidas |

### Interpretação

| Pontuação | Classificação | Ação |
|-----------|---------------|------|
| 90-100 | 🥇 Excelente | Pronto para Nature, Science, Physical Review |
| 75-89 | 🥈 Muito Bom | Pronto para Qualis A1/A2 com pequenos ajustes |
| 60-74 | 🥉 Bom | Melhorias moderadas necessárias (2-4 semanas) |
| 40-59 | ⚠️ Insuficiente | Trabalho substancial necessário (1-2 meses) |
| 0-39 | ❌ Inadequado | Repensar abordagem (3+ meses) |

---

## 🎓 Exemplos de Uso

### Caso 1: Beneficial Quantum Noise

```bash
# 1. Analisar código do projeto
python enhanced_code_analyzer.py .

# Output:
# ✅ 7 ansätze encontrados
# ✅ 3 modelos de ruído encontrados
# ✅ 3.360 configurações calculadas

# 2. Gerar artigo MODE A (internacional)
python gerador_artigo_completo.py \
    --repositorio . \
    --output artigo_quantum_v1 \
    --periodico-primario "Quantum" \
    --mode A \
    --policy R1

# 3. Verificar qualidade
python tools/avaliar_qualidade.py artigo_quantum_v1/

# Output: Pontuação: 87/100 (🥈 Muito Bom)
```

### Caso 2: Projeto de ML Clássico

```bash
# Adaptar para domínio de ML clássico
python gerador_artigo_completo.py \
    --repositorio /path/to/ml_project \
    --output artigo_ml \
    --periodico-primario "Journal of Machine Learning Research" \
    --mode A \
    --policy R1

# O sistema adapta automaticamente:
# - Ansätze → Arquiteturas de rede
# - Qubits → Neurônios/camadas
# - Ruído quântico → Dropout/regularização
```

---

## 🔧 Troubleshooting

### Problema: "Código não tem seeds fixas"

**Solução:**
1. Documente como `[INFORMAÇÃO AUSENTE]`
2. Execute 10-30 vezes e reporte média ± desvio padrão
3. Adicione à seção "Threats to Validity"

### Problema: "Muitas lacunas de citação (R0)"

**Solução:**
```bash
# Mudar para R1
sed -i 's/"reference_policy": "R0"/"reference_policy": "R1"/g' config.json

# Regenerar Fase 3
python gerador_artigo_completo.py --regenerate fase3_literatura
```

### Problema: "Quality Gate falhou"

**Checklist:**
- [ ] Todos os números têm evidência?
- [ ] Tabelas estão completas?
- [ ] DOIs presentes?
- [ ] Contrapontos incluídos?

Ver `FAQ_TROUBLESHOOTING.md` para mais soluções.

---

## 📈 Estatísticas do Sistema

- **Templates**: 8 arquivos profissionais
- **Documentação**: 4 guias completos (>40K palavras)
- **Glossário**: 90+ termos definidos
- **Checklist**: 100 pontos objetivos
- **Fases**: 6 com quality gates
- **Outputs**: 24+ arquivos gerados
- **Análise automática**: 7 componentes extraídos

---

## 🤝 Contribuindo

Este é um framework vivo. Contribuições são bem-vindas:

1. **Issues**: Reporte bugs ou sugira melhorias
2. **Pull Requests**: Envie templates, ferramentas, melhorias
3. **Documentação**: Ajude a melhorar guias e exemplos

---

## 📄 Licença

MIT License - Ver arquivo LICENSE

---

## 📧 Contato

**Issues**: https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers/issues  
**Documentação**: Ver arquivos `.md` na raiz do repositório

---

## 🎉 Agradecimentos

Este sistema foi desenvolvido para facilitar a geração de artigos científicos de alto impacto com rigor metodológico e transparência total.

**Versão**: 1.0  
**Última Atualização**: 26/12/2025  
**Status**: ✅ Pronto para uso

---

**Próximos Passos:**
1. Leia o [GLOSSARIO.md](GLOSSARIO.md) para familiarizar-se com os termos
2. Execute `python enhanced_code_analyzer.py .` no seu projeto
3. Edite `config.json` com suas preferências
4. Execute `python gerador_artigo_completo.py --repositorio . --output meu_artigo`
5. Revise os arquivos gerados em `meu_artigo/`
6. Use o [CHECKLIST_AUDITORIA_100PTS.md](CHECKLIST_AUDITORIA_100PTS.md) para avaliar qualidade

**Boa sorte na submissão!** 🚀
