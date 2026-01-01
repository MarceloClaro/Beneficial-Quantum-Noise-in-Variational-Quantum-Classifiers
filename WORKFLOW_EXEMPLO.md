# MegaPrompt v2.0 - Exemplo de Workflow Completo

Este documento demonstra um workflow completo de geração de artigo científico Qualis A1 usando o MegaPrompt v2.0.

## 📋 Pré-requisitos

- Python 3.7+
- Git
- Repositório com código, dados e experimentos
- 6-10 dias úteis disponíveis


## 🚀 Workflow Passo a Passo

### Passo 0: Configuração Inicial (30 min)

```bash

# Clone o repositório
git clone <https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers.git>
cd Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers

# Instale dependências
pip install -r requirements.txt

# Configure o projeto
cp config.json config.json.backup  # Backup se já existir
nano config.json

```text

**Edite `config.json`:**


```json
{
  "output_mode": "MODE_A",
  "reference_policy": "R1",
  "editorial_profile": "PROFILE_PR_QUANTUM",
  "target_journals": {
    "primary": "Quantum",
    "secondary": ["Physical Review A", "npj Quantum Information"]
  },
  "user_inputs": {
    "research_question": "Can quantum noise be beneficial for VQCs?",
    "hypotheses": [
      "H₀: Optimal noise level p* > 0 improves generalization"
    ]
  }
}

```text

### Fase 1: Auditoria Técnica (8-12h)

**Objetivo**: Inventariar todo o código e experimentos


```bash

# Examine o código
python enhanced_code_analyzer.py --output artigo_cientifico/fase1_analise/

# Gere manifesto de execução
python -c "
from qualis_a1_modules.reproducibility import gerar_manifesto_execucao, configurar_seeds_reprodutiveis

seed_info = configurar_seeds_reprodutiveis(seed=42)
config = {'version': '8.0-QAI', 'default_seed': 42}
gerar_manifesto_execucao(config, 'artigo_cientifico/fase1_analise/', seed_info)
"

# Calcule total de configurações
python tools/megaprompt_v2/generate_s1.py --verify-only

```text

**Arquivos Gerados**:
- ✅ `analise_codigo_inicial.md`
- ✅ `manifesto_execucao.json`
- ✅ Contagem de configurações verificada


**Quality Gate F1**:
- [ ] Cada item tem origem
- [ ] Total calculado e conferido
- [ ] Ausências explicitadas


### Fase 2: Bibliografia (6-10h)

**Objetivo**: Compilar referências e síntese da literatura


```bash

# Revise referências existentes
cat artigo_cientifico/fase2_bibliografia/referencias_compiladas.md

# Se R1, busque novas referências em:
# - Google Scholar
# - arXiv
# - Web of Science
# - Scopus

# Organize por categorias:
# 1. Quantum noise models
# 2. VQC architectures  
# 3. Optimization methods
# 4. Statistical techniques
# 5. Related work
# 6. Datasets
# 7. Applications

```text

**Arquivos Gerados**:
- ✅ `referencias_compiladas.md` (com DOI)
- ✅ `sintese_literatura.md` (crítica)
- ✅ `taxonomia_estado_da_arte.md`


**Quality Gate F2**:
- [ ] Cada técnica tem referência ou [LACUNA]
- [ ] Síntese contém contraste


### Fase 3: Estrutura do Artigo (4-6h)

**Objetivo**: Definir problema formal e hipóteses


```bash

# Revise templates
cat templates/problema_formal_template.md

# Edite formulação matemática
nano artigo_cientifico/fase3_estrutura/problema_formal.md

```text

**Conteúdo do `problema_formal.md`**:


```markdown

# Formal Problem Statement

#### Seja:
- $\mathcal{D} = \{(x_i, y_i)\}_{i=1}^N$ dataset com $N$ amostras
- $U(\theta)$ circuito quântico parametrizado
- $\mathcal{N}_p(\cdot)$ canal de ruído com parâmetro $p$
- $L(\theta, p)$ função de custo


**O problema é:**

$$
(\theta^*, p^*) = \arg\min_{\theta, p} \mathbb{E}_{(x,y) \sim \mathcal{D}_{test}} [L(y, f(x; \theta, p))]
$$

**Hipótese Principal (H₀):**

$$
\exists p^* > 0 : \text{Acc}(p^*) > \text{Acc}(0)
$$

```text

**Arquivos Gerados**:
- ✅ `problema_formal.md`
- ✅ `titulos_palavras_chave.md`
- ✅ `hipoteses_objetivos.md`


**Quality Gate F3**:
- [ ] Problema formal compatível com código
- [ ] Cada hipótese tem teste correspondente


### Fase 4: Redação (20-30h)

**Objetivo**: Escrever todas as seções do artigo


```bash

# Edite seções principais
nano artigo_cientifico/fase4_secoes/metodologia_completa.md

```text

**Adicione Algorithm 1**:


```latex
\begin{algorithm}[H]
\caption{Experimental Pipeline}
\begin{algorithmic}[1]
\REQUIRE Datasets, Configurations, Seeds
\STATE Initialize results table
\FOR{each configuration}
  \STATE Train model with noise
  \STATE Evaluate on test set
  \STATE Record metrics
\ENDFOR
\RETURN Results table
\end{algorithmic}
\end{algorithm}

```text

**Crie Tabela Código→Método**:


```markdown
| Componente | Arquivo | Função | Linha | Parâmetros |
|-----------|---------|--------|-------|------------|
| Depolarizing Channel | framework.py | RuidoDepolarizante | 1523-1548 | nivel: float |
| QNG Optimizer | framework.py | ClassificadorVQC | 2341-2398 | lr: float |

```text

**Arquivos Gerados**:
- ✅ `resumo_abstract.md`
- ✅ `introducao_completa.md`
- ✅ `revisao_literatura_completa.md`
- ✅ `metodologia_completa.md` (com Algorithm 1)
- ✅ `resultados_completo.md`
- ✅ `discussao_completa.md`
- ✅ `conclusao_completa.md`
- ✅ `agradecimentos_referencias.md`


**Quality Gate F4**:
- [ ] Sem números sem lastro
- [ ] Methods completo


### Fase 5: Material Suplementar (8-12h)

**Objetivo**: Criar tabelas e figuras suplementares


```bash

# Gere Tabela S1
python tools/megaprompt_v2/generate_s1.py \

  --config config.json \
  --output artigo_cientifico/fase5_suplementar/tabela_s1_configuracoes.csv


# Verifique geração
wc -l artigo_cientifico/fase5_suplementar/tabela_s1_configuracoes.csv

# Deve mostrar número correto de linhas + 1 (header)

```text

**Crie outras tabelas**:


```markdown

## Tabela S2: Comparação com Estado da Arte

| Método | Acurácia | Tempo (s) | Framework | Referência |
|--------|----------|-----------|-----------|------------|
| Nosso (melhor) | 65.83% | 234 | PennyLane | - |
| Du et al. (2021) | 62.10% | [NÃO DISPONÍVEL] | Qiskit | (Du et al., 2021) |
| Smith et al. (2022) | 61.50% | 180 | Cirq | (Smith et al., 2022) |

```text

**Arquivos Gerados**:
- ✅ `tabelas_suplementares.md`
- ✅ `tabela_s1_configuracoes.csv`
- ✅ `figuras_suplementares.md`
- ✅ `apendice_suplementar.md`


**Quality Gate F5**:
- [ ] S1 confere com total calculado
- [ ] Cada tabela/figura rastreável


### Fase 6: Consolidação (6-8h)

**Objetivo**: Consolidar e validar


```bash

# 1. Consolide o manuscrito
bash tools/megaprompt_v2/build_paper.sh

# Saída esperada:
# ✅ Manuscript: artigo_cientifico/fase6_consolidacao/manuscrito_internacional_final.md
# ✅ Summary: artigo_cientifico/fase6_consolidacao/sumario_executivo.md
# Words: 22620
# Lines: 2433

# 2. Verifique consistência
python tools/megaprompt_v2/check_consistency.py

# Meta: ≥ 95% consistency
# Se < 95%, revisar issues no relatório

# 3. Execute auditoria
python tools/megaprompt_v2/audit_checklist.py

# Meta: ≥ 90/100 pontos
# Categorias: Reprodutibilidade, Rastreabilidade, Estatística, Transparência

# 4. Gere tabela de rastreabilidade
nano artigo_cientifico/fase6_consolidacao/rastreabilidade_completa.md

```text

**Tabela de Rastreabilidade**:


```markdown
| Seção | Afirmação/Número | Evidência | Referência |
|-------|------------------|-----------|------------|
| Abstract | 65.83% accuracy | resultados/melhor_config.json:accuracy | - |
| Methods | Lindblad equation | framework.py:L1523 | (Lindblad, 1976) |
| Results | p < 0.001 | resultados/anova.json:p_value | - |

```text

**Arquivos Gerados**:
- ✅ `rastreabilidade_completa.md`
- ✅ `relatorio_consistencia.md`
- ✅ `checklist_auditoria_100pts.md`
- ✅ `manuscrito_internacional_final.md`
- ✅ `sumario_executivo.md`


**Quality Gate Final**:
- [ ] Consistência ≥ 95%
- [ ] Auditoria ≥ 90/100
- [ ] Citação↔referência 100%
- [ ] Reprodutibilidade completa


## 📊 Verificação Final

```bash

# Verifique arquivos gerados
tree artigo_cientifico/

# Esperado:
# artigo_cientifico/
# ├── fase1_analise/ (2 arquivos)
# ├── fase2_bibliografia/ (3 arquivos)
# ├── fase3_estrutura/ (3 arquivos)
# ├── fase4_secoes/ (8 arquivos)
# ├── fase5_suplementar/ (4 arquivos)
# └── fase6_consolidacao/ (5 arquivos)
#
# Total: ~25 arquivos

```text

## ✅ Checklist Pré-Submissão

### Documentação
- [ ] Manuscrito final revisado
- [ ] Todas as seções completas
- [ ] Figuras e tabelas numeradas
- [ ] Referências formatadas corretamente


### Qualidade
- [ ] Consistência código-texto ≥ 95%
- [ ] Auditoria ≥ 90/100 pontos
- [ ] Todos os marcadores [INFORMAÇÃO AUSENTE] revisados
- [ ] Tabela de rastreabilidade completa


### Reprodutibilidade
- [ ] `requirements.txt` atualizado
- [ ] Seeds documentadas
- [ ] Comandos de execução incluídos
- [ ] Manifesto de execução gerado


### Material Suplementar
- [ ] Tabela S1 com todas as configurações
- [ ] Tabelas S2-S5 completas
- [ ] Figuras S1-S8 descritas
- [ ] Código disponível publicamente


### Submissão
- [ ] Periódico-alvo escolhido
- [ ] Formato específico verificado
- [ ] Carta de apresentação preparada
- [ ] Todos os arquivos empacotados


## 🎯 Comandos Rápidos

```bash

# Workflow completo automatizado
make all  # Se houver Makefile

# Ou manualmente:
python tools/megaprompt_v2/generate_s1.py && \
bash tools/megaprompt_v2/build_paper.sh && \
python tools/megaprompt_v2/check_consistency.py && \
python tools/megaprompt_v2/audit_checklist.py

# Ver relatórios
cat artigo_cientifico/fase6_consolidacao/sumario_executivo.md
cat artigo_cientifico/fase6_consolidacao/checklist_auditoria_100pts.md

```text

## 📈 Métricas de Sucesso

Ao final do workflow, você deve ter:

- ✅ **100%** de rastreabilidade (Seção → Evidência → Origem)
- ✅ **≥95%** de consistência código-texto
- ✅ **≥90/100** pontos no audit checklist
- ✅ **22,000+** palavras no manuscrito
- ✅ **45+** referências compiladas
- ✅ **14+** tabelas (9 principais + 5 suplementares)
- ✅ **8+** figuras suplementares
- ✅ **25+** arquivos de documentação


## 🐛 Troubleshooting Comum

### Problema: Consistência < 95%

**Solução:**

```bash

# Identifique problemas específicos
cat artigo_cientifico/fase6_consolidacao/relatorio_consistencia.md

# Corrija:
# - Adicione citações para números
# - Verifique arquivos referenciados
# - Resolva marcadores [INFORMAÇÃO AUSENTE]

# Re-execute
python tools/megaprompt_v2/check_consistency.py

```text

### Problema: Auditoria < 90 pontos

**Solução:**

```bash

# Veja detalhamento por categoria
cat artigo_cientifico/fase6_consolidacao/checklist_auditoria_100pts.md

# Foque em categorias com < 70%
# Reprodutibilidade: Adicione seeds, ambiente, comandos
# Rastreabilidade: Complete tabelas
# Estatística: Adicione testes, correções, IC, effect sizes
# Transparência: Documente limitações, disponibilize dados

```text

### Problema: Build falha

**Solução:**

```bash

# Verifique se todas as seções existem
ls -la artigo_cientifico/fase4_secoes/

# Se faltando, crie esqueleto:
touch artigo_cientifico/fase4_secoes/{resumo_abstract,introducao_completa,metodologia_completa,resultados_completo,discussao_completa,conclusao_completa,agradecimentos_referencias}.md

# Re-execute
bash tools/megaprompt_v2/build_paper.sh

```

## 📚 Recursos Adicionais

- 📖 [MegaPrompt v2.0 README](MEGAPROMPT_V2_README.md)
- 💡 [Exemplos Práticos](EXEMPLOS_PRATICOS.md)
- ❓ [FAQ e Troubleshooting](FAQ_TROUBLESHOOTING.md)
- 📚 [Glossário de Termos](GLOSSARIO.md)
- 🔧 [Documentação de Ferramentas](tools/megaprompt_v2/README.md)


## 🎓 Dicas Finais

1. **Seja Paciente**: 6-10 dias é normal para um artigo Qualis A1 de qualidade
2. **Use Quality Gates**: Não pule fases, valide cada etapa
3. **Documente Tudo**: Melhor excesso que falta de documentação
4. **Peça Revisão**: Tenha alguém revisando seu trabalho
5. **Mantenha Backup**: `git commit` frequentemente


---


**Boa sorte com sua submissão! 🚀**


*MegaPrompt v2.0 - Framework para Artigos Científicos Qualis A1*

