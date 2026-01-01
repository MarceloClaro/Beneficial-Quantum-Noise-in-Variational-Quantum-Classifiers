# Guia de Uso do Framework de Geração de Artigos Científicos QUALIS A1

**Versão:** 1.0  
**Data:** Dezembro 2025  
**Framework:** Beneficial Quantum Noise in Variational Quantum Classifiers


---


## 🎯 Visão Geral

Este framework permite gerar artigos científicos completos e rigorosos, prontos para submissão a periódicos de alto impacto (Nature, Science, Quantum, Physical Review, Qualis A1), com **100% de conivência entre código/dados e texto**.

---


## 📋 Pré-requisitos

### Ambiente Python

```bash

# Criar ambiente virtual
python -m venv venv
source venv/bin/activate  # Linux/Mac

# ou
venv\Scripts\activate  # Windows

# Instalar dependências
pip install -r requirements.txt

```text

### Estrutura do Repositório

```

.
├── framework_investigativo_completo.py  # Código principal
├── artigo_cientifico/                   # Artigo gerado (6 fases)
├── tools/                               # Ferramentas de validação
├── MEGA_PROMPT_QUALIS_A1.md            # Especificação completa
└── WORKFLOW_ARTIGO.md                   # Este arquivo

```text

---


## 🚀 Workflows

### Workflow 1: Validar Artigo Existente

**Objetivo:** Verificar conformidade do artigo já gerado com critérios QUALIS A1.


```bash

# 1. Validar conformidade QUALIS A1
python tools/validate_qualis_a1.py \

    --article artigo_cientifico/ \
    --report VALIDATION_REPORT.md


# 2. Verificar conivência código-texto
python tools/verify_code_text_congruence.py \

    --code framework_investigativo_completo.py \
    --article artigo_cientifico/ \
    --output CONGRUENCE_REPORT.md


# 3. Revisar relatórios
cat VALIDATION_REPORT.md
cat CONGRUENCE_REPORT.md

```text

#### Output Esperado:
- Pontuação QUALIS A1: 80-100 (Bom a Excelente)
- Congruência: ≥95%


---


### Workflow 2: Gerar Novo Artigo do Zero

**Objetivo:** Gerar artigo completo a partir de código/dados.


```bash

# 1. Executar gerador completo
python gerador_artigo_completo.py \

    --repositorio . \
    --output artigo_gerado \
    --periodico-primario "Nature Communications"


# 2. Validar artigo gerado
python tools/validate_qualis_a1.py \

    --article artigo_gerado/ \
    --report artigo_gerado/VALIDATION_REPORT.md


# 3. Verificar congruência
python tools/verify_code_text_congruence.py \

    --code framework_investigativo_completo.py \
    --article artigo_gerado/ \
    --output artigo_gerado/CONGRUENCE_REPORT.md

```text

**Fases Executadas:**


1. **Fase 1:** Análise Inicial e Planejamento (3 horas)
   - `analise_codigo_inicial.md`
   - `linha_de_pesquisa.md`


2. **Fase 2:** Pesquisa Bibliográfica Profunda (4 horas)
   - `referencias_compiladas.md`
   - `sintese_literatura.md`


3. **Fase 3:** Elaboração da Estrutura (2 horas)
   - `titulos_palavras_chave.md`
   - `hipoteses_objetivos.md`


4. **Fase 4:** Redação das Seções (30-40 horas)
   - `resumo_abstract.md`
   - `introducao_completa.md`
   - `revisao_literatura_completa.md`
   - `metodologia_completa.md`
   - `resultados_completo.md`
   - `discussao_completa.md`
   - `conclusao_completa.md`
   - `agradecimentos_referencias.md`


5. **Fase 5:** Material Suplementar (5-8 horas)
   - `tabelas_suplementares.md`
   - `figuras_suplementares.md`
   - `notas_metodologicas_adicionais.md`


6. **Fase 6:** Consolidação e Verificação (3-5 horas)
   - `relatorio_conivencia.md`
   - `artigo_completo_final.md`
   - `sumario_executivo.md`


---


### Workflow 3: Atualizar Artigo com Novos Resultados

**Objetivo:** Incorporar novos experimentos ao artigo existente.


```bash

# 1. Executar novos experimentos
python framework_investigativo_completo.py \

    --resultados resultados_novos


# 2. Atualizar seção de resultados
python atualizar_artigos_com_resultados.py \

    --resultados resultados_novos \
    --artigo artigo_cientifico/fase4_secoes/resultados_completo.md


# 3. Validar atualização
python tools/validate_qualis_a1.py \

    --article artigo_cientifico/ \
    --report VALIDATION_REPORT_UPDATED.md


# 4. Verificar congruência
python tools/verify_code_text_congruence.py \

    --code framework_investigativo_completo.py \
    --article artigo_cientifico/ \
    --output CONGRUENCE_REPORT_UPDATED.md

```text

---


### Workflow 4: Preparar para Submissão

**Objetivo:** Finalizar artigo para submissão a periódico específico.


```bash

# 1. Consolidar artigo completo
cd artigo_cientifico/fase6_consolidacao
cat ../fase4_secoes/*.md > artigo_completo_consolidado.md

# 2. Converter para LaTeX (se necessário)
# Para npj Quantum Information (recomendado)
pandoc artigo_completo_consolidado.md \

    -o npj_qi_submission.tex \
    --template=../latex_template/npj_qi_submission.tex \
    --bibliography=../fase4_secoes/referencias.bib


# 3. Gerar figuras finais
python ../../visualize_results.py \

    --output figuras_finais/


# 4. Verificar checklist final
python ../../tools/validate_qualis_a1.py \

    --article ../ \
    --report FINAL_VALIDATION.md


# 5. Criar pacote de submissão
mkdir submission_package
cp artigo_completo_consolidado.md submission_package/
cp -r ../fase5_suplementar submission_package/
cp figuras_finais/* submission_package/
cp FINAL_VALIDATION.md submission_package/

```text

---


## 🔍 Interpretando os Relatórios

### Relatório de Validação QUALIS A1

#### Pontuação:
- **90-100:** 🥇 EXCELENTE - Pronto para Nature/Science/Quantum
- **80-89:** 🥈 BOM - Pronto para periódicos Qualis A1
- **70-79:** 🥉 REGULAR - Necessita ajustes
- **<70:** ❌ INSUFICIENTE - Revisão substancial necessária


**Critérios Avaliados:**


1. **Estruturais:**
   - Número de referências (35-50)
   - Cobertura DOI (≥80%)
   - Hipóteses (≥3)
   - Objetivos SMART (≥3)


2. **Extensão (palavras):**
   - Abstract: 250-300
   - Introdução: 3.000-4.000
   - Revisão: 4.000-5.000
   - Metodologia: 4.000-5.000
   - Resultados: 3.000-4.000
   - Discussão: 4.000-5.000
   - Conclusão: 1.000-1.500


3. **Qualidade:**
   - Tabelas (≥5)
   - Equações LaTeX (≥10)
   - Conivência código-texto (≥95%)


### Relatório de Congruência Código-Texto

#### Congruência Geral:
- **≥95%:** ✅ EXCELENTE - Conivência total
- **85-94%:** ✅ BOA - Conivência adequada
- **70-84%:** ⚠️ REGULAR - Algumas inconsistências
- **<70%:** ❌ INSUFICIENTE - Revisão necessária


#### Componentes Verificados:
- Classes implementadas vs. mencionadas
- Datasets utilizados vs. descritos
- Modelos de ruído implementados vs. documentados
- Ansätze quânticos implementados vs. descritos
- Métricas de avaliação utilizadas vs. reportadas
- Bibliotecas e versões
- Seeds de reprodutibilidade


---


## 🛠️ Ferramentas Disponíveis

### 1. Validador QUALIS A1

**Arquivo:** `tools/validate_qualis_a1.py`


**Uso:**

```bash
python tools/validate_qualis_a1.py \

    --article artigo_cientifico/ \
    --report VALIDATION_REPORT.md

```text

#### Saída:
- Pontuação geral (0-100)
- Tabela de conformidade por critério
- Recomendações de melhoria


### 2. Verificador de Congruência

**Arquivo:** `tools/verify_code_text_congruence.py`


**Uso:**

```bash
python tools/verify_code_text_congruence.py \

    --code framework_investigativo_completo.py \
    --article artigo_cientifico/ \
    --output CONGRUENCE_REPORT.md

```text

#### Saída:
- Percentual de congruência geral
- Análise detalhada por componente
- Inconsistências detectadas


### 3. Gerador de Artigo Completo

**Arquivo:** `gerador_artigo_completo.py`


**Uso:**

```bash
python gerador_artigo_completo.py \

    --repositorio . \
    --output artigo_gerado \
    --periodico-primario "Nature Communications"

```text

#### Saída:
- Estrutura completa de 6 fases
- 20 documentos markdown
- Relatórios de qualidade


### 4. Consultor Metodológico

**Arquivo:** `consultor_metodologico.py`


**Uso:**

```bash
python consultor_metodologico.py \

    --insumos exemplo_insumos_consultor.json \
    --output parecer_metodologico.md

```text

#### Saída:
- Parecer metodológico detalhado
- Sugestões de melhoria
- Análise estatística


---


## 📚 Documentação Adicional

### Arquivos de Referência

1. **MEGA_PROMPT_QUALIS_A1.md**
   - Especificação completa do framework
   - Todos os critérios QUALIS A1
   - Exemplos de uso


2. **GUIA_COMPLETO_GERACAO_ARTIGOS.md**
   - Tutorial passo a passo
   - Exemplos práticos
   - Troubleshooting


3. **GERADOR_ARTIGO_README.md**
   - Uso do gerador automático
   - Parâmetros e configurações


4. **CONSULTOR_METODOLOGICO_README.md**
   - Análise de qualidade metodológica
   - Pareceres e recomendações


### Estrutura do Artigo Gerado

```

artigo_cientifico/
├── fase1_analise/
│   ├── analise_codigo_inicial.md
│   └── linha_de_pesquisa.md
├── fase2_bibliografia/
│   ├── referencias_compiladas.md
│   └── sintese_literatura.md
├── fase3_estrutura/
│   ├── titulos_palavras_chave.md
│   └── hipoteses_objetivos.md
├── fase4_secoes/
│   ├── resumo_abstract.md
│   ├── introducao_completa.md
│   ├── revisao_literatura_completa.md
│   ├── metodologia_completa.md
│   ├── resultados_completo.md
│   ├── discussao_completa.md
│   ├── conclusao_completa.md
│   └── agradecimentos_referencias.md
├── fase5_suplementar/
│   ├── tabelas_suplementares.md
│   ├── figuras_suplementares.md
│   └── notas_metodologicas_adicionais.md
└── fase6_consolidacao/
    ├── relatorio_conivencia.md
    ├── rastreabilidade_completa.md
    ├── artigo_completo_final.md
    └── sumario_executivo.md

```text

---


## ✅ Checklist de Validação Final

### Antes da Submissão

- [ ] Pontuação QUALIS A1 ≥ 80
- [ ] Congruência código-texto ≥ 95%
- [ ] Todas as hipóteses testadas
- [ ] Todos os objetivos atingidos
- [ ] Referências completas (35-50)
- [ ] DOI para ≥80% das referências
- [ ] Limitações honestamente discutidas
- [ ] Trabalhos futuros propostos
- [ ] Código e dados públicos (GitHub)
- [ ] Material suplementar completo
- [ ] Figuras em alta resolução (300 DPI)
- [ ] Formatação LaTeX (se aplicável)
- [ ] Revisão por pares interna
- [ ] English language editing (se internacional)


### Documentos para Submissão

- [ ] Artigo principal (PDF/LaTeX)
- [ ] Material suplementar
- [ ] Cover letter
- [ ] Competing interests statement
- [ ] Author contributions
- [ ] Data availability statement
- [ ] Code availability statement
- [ ] Figuras individuais (alta resolução)
- [ ] Tabelas suplementares (CSV/Excel)


---


## 🎯 Periódicos-Alvo Recomendados

### Para Este Framework (Computação Quântica)

1. **npj Quantum Information** ⭐⭐⭐
   - Impact Factor: 7.6
   - Open Access
   - Nature Portfolio
   - **Mais recomendado** para este trabalho


2. **Nature Communications** ⭐⭐⭐
   - Impact Factor: 14.9
   - Open Access
   - Multidisciplinar
   - Alta visibilidade


3. **Quantum** ⭐⭐
   - Impact Factor: 5.1
   - Open Access
   - Especializado em computação quântica
   - Revisão rápida


4. **Physical Review A** ⭐⭐
   - Impact Factor: 2.9
   - Rigor técnico alto
   - Comunidade de física quântica


5. **Physical Review Research** ⭐
   - Impact Factor: 4.2
   - Open Access
   - Multidisciplinar


---


## 🆘 Troubleshooting

### Problema: Pontuação QUALIS A1 Baixa

**Causa:** Seções com extensão fora do esperado


**Solução:**

```bash

# Identificar seções problemáticas
python tools/validate_qualis_a1.py --article artigo_cientifico/ --report report.md
cat report.md | grep "⚠️\|❌"

# Expandir seções curtas ou condensar seções longas
# Editar manualmente os arquivos .md em fase4_secoes/

```text

### Problema: Congruência Código-Texto Baixa

**Causa:** Componentes do código não mencionados no artigo


**Solução:**

```bash

# Identificar inconsistências
python tools/verify_code_text_congruence.py \

    --code framework_investigativo_completo.py \
    --article artigo_cientifico/ \
    --output congruence.md


# Atualizar seções relevantes (principalmente Metodologia e Resultados)
# para incluir todos os componentes implementados

```text

### Problema: Referências Insuficientes

**Causa:** Menos de 35 referências


**Solução:**

```bash

# Buscar referências adicionais por categoria:
# - Trabalhos fundacionais (5-8)
# - Estado da arte recente (8-12)
# - Metodologia (6-10)
# - Análise estatística (4-6)
# - Frameworks (3-5)
# - Trabalhos críticos (3-5)
# - Aplicações (3-5)

# Atualizar fase2_bibliografia/referencias_compiladas.md

```text

### Problema: Equações LaTeX Não Renderizadas

**Causa:** Sintaxe LaTeX incorreta


**Solução:**

```bash

# Validar sintaxe LaTeX
grep -n '\$\$' artigo_cientifico/fase4_secoes/metodologia_completa.md

# Formatos válidos:
# - Inline: $equação$
# - Display: $$equação$$
# - Numbered: \begin{equation} equação \end{equation}

```text

---


## 📞 Suporte

### Documentação

- GitHub: <https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers>
- Issues: <https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers/issues>


### Citação

Se usar este framework, cite:

```bibtex
@software{claro2025beneficial,
  title={Beneficial Quantum Noise in Variational Quantum Classifiers:
         A Framework for Qualis A1 Scientific Article Generation},
  author={Claro, Marcelo},
  year={2025},
  url={<https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers}>
}

```

---


**Versão:** 1.0  
**Data:** Dezembro 2025  
**Próxima Revisão:** Após submissão a periódico

