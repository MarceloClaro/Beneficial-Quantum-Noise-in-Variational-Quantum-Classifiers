# Instruções para Reproduzir os Resultados QUALIS A1

**Data deste Guia:** 23 de dezembro de 2025  
**Framework Version:** 7.2


---


## 🎯 Objetivo

Este documento fornece instruções passo-a-passo para reproduzir os resultados e gráficos atualizados do framework, conforme documentado em `RESULTADOS_ATUALIZADOS_QUALIS_A1.md`.

---


## ⚡ Execução Rápida (15-30 minutos)

### Modo Quick Bayesian (Validação)

Esta foi a configuração usada para gerar os resultados documentados:

```bash

# 1. Configurar ambiente
export VQC_QUICK=1
export VQC_BAYESIAN=1

# 2. Executar framework
python framework_investigativo_completo.py \

    --bayes \
    --trials 5 \
    --dataset-bayes moons

```text

**Tempo estimado:** 4-5 minutos  
**Trials executados:** 5  
**Épocas por trial:** 5  
**Melhor acurácia esperada:** ~80%


### Resultados Gerados

Após a execução, você encontrará:

```

resultados_YYYY-MM-DD_HH-MM-SS/
├── otimizacao_bayesiana/
│   ├── resultado_otimizacao.json      # Configuração ótima
│   ├── historico_trials.csv           # Histórico completo
│   └── README_otimizacao.md           # Documentação
├── figura2_beneficial_noise.html      # Visualização interativa
├── analises_estatisticas_completo.csv # Análises estatísticas
├── comparacao_baselines.csv           # Comparação VQC vs SVM/RF
└── metadata.json                      # Metadados da execução

```text

---


## 🚀 Execução Completa Para QUALIS A1

### Opção 1: Modo Bayesian Completo (Recomendado - 1-2 horas)

```bash
python framework_investigativo_completo.py \

    --bayes \
    --trials 200 \
    --dataset-bayes all \
    --epocas-bayes 15

```text

#### Vantagens:
- 10-20x mais rápido que Grid Search
- Encontra configurações ótimas inteligentemente
- Análise de importância de hiperparâmetros
- Ideal para pesquisa exploratória


**Tempo estimado:** 1-2 horas  
**Trials:** 200  
**Datasets:** Todos (moons, circles, iris, breast_cancer, wine)  
**Épocas:** 15 por trial  


### Opção 2: Modo Grid Search Completo (Rigor Máximo - 15-20 horas)

```bash
python framework_investigativo_completo.py

```text

#### Vantagens:
- Cobertura exhaustiva do espaço de hiperparâmetros
- 8,280 configurações testadas
- Máxima reprodutibilidade científica
- Ideal para artigos QUALIS A1


**Tempo estimado:** 15-20 horas  
**Configurações:** 8,280  
**Épocas:** 15 por configuração  


### Opção 3: Modo Híbrido (Máxima Precisão - 20-25 horas)

```bash
python framework_investigativo_completo.py \

    --run-both \
    --trials 150

```text

#### Vantagens:
- Combina exploração (Grid) com refinamento (Bayesiano)
- Melhor de ambos os mundos
- Recomendado para trabalhos definitivos


---


## 📊 Verificação dos Resultados

### Validar Execução

```bash

# 1. Verificar se o diretório de resultados foi criado
ls -la resultados_*/

# 2. Verificar arquivo de configuração ótima
cat resultados_*/otimizacao_bayesiana/resultado_otimizacao.json

# 3. Ver histórico de trials
head -20 resultados_*/otimizacao_bayesiana/historico_trials.csv

# 4. Abrir visualização interativa
firefox resultados_*/figura2_beneficial_noise.html

```text

### Métricas Esperadas

Se a execução foi bem-sucedida, você deve ver:

#### Para Quick Bayesian (5 trials):
- Melhor acurácia: ~75-85%
- Trials completos: 5/5
- Arquitetura ótima: Strongly Entangling ou Hardware Efficient
- Tipo de ruído: Depolarizante
- Nível ótimo: ~0.001-0.005


#### Para Bayesian Completo (200 trials):
- Melhor acurácia: ~85-95%
- Trials completos: 180-200 (alguns podem ser podados)
- Convergência clara para configuração ótima
- Análise de importância bem definida


#### Para Grid Search Completo:
- Total de experimentos: 8,280
- CSVs individuais gerados: 8,280
- Todas as 9 figuras geradas
- Análises estatísticas completas


---


## 🔧 Solução de Problemas

### Erro: "ModuleNotFoundError: No module named 'numpy'"

```bash
pip install -r requirements.txt

```text

### Erro: "Image export using kaleido engine requires..."

```bash
pip install --upgrade kaleido

```text

### Execução muito lenta

```bash

# Use modo Quick ou reduza trials
export VQC_QUICK=1
python framework_investigativo_completo.py --bayes --trials 10

```text

### Memória insuficiente

```bash

# Execute para um dataset de cada vez
python framework_investigativo_completo.py \

    --bayes \
    --trials 100 \
    --dataset-bayes moons

```text

---


## 📁 Estrutura de Resultados Esperada

Após execução completa:

```

resultados_YYYY-MM-DD_HH-MM-SS/
├── README.md
├── metadata.json
├── metadata_orchestrator.json
│
├── # Resultados Principais
├── resultados_completos_artigo.csv
├── comparacao_baselines.csv
├── analise_comparacao_inicializacoes.csv
├── analises_estatisticas_completo.csv
│
├── # Visualizações (9 figuras)
├── figura2_beneficial_noise.html
├── figura2b_beneficial_noise_ic95.html
├── figura3_noise_types.html
├── figura3b_noise_types_ic95.html
├── figura4_initialization.html
├── figura5_architecture_tradeoffs.html
├── figura6_effect_sizes.html
├── figura7_overfitting.html
├── figura_correlacao.html
│
├── # Análises Bayesianas
├── otimizacao_bayesiana/
│   ├── resultado_otimizacao.json
│   ├── historico_trials.csv
│   └── README_otimizacao.md
│
├── # Artefatos
├── circuitos/                    # Diagramas de circuitos quânticos
├── barren_plateaus/              # Análise de gradientes
├── experimentos_individuais/     # CSVs por experimento (8,280)
├── analises_individuais/         # Análises granulares
└── visualizacoes_individuais/    # Gráficos individuais

```text

---


## 📈 Análise dos Resultados

### Explorar Resultados Interativamente

```python
import pandas as pd
import json

# Carregar resultados
with open('resultados_*/otimizacao_bayesiana/resultado_otimizacao.json') as f:
    resultados = json.load(f)

# Melhor configuração
print("Melhor Acurácia:", resultados['melhor_acuracia'])
print("Melhor Configuração:", resultados['melhor_params'])

# Importância dos hiperparâmetros
print("\nImportância:")
for param, imp in sorted(resultados['importancias'].items(),
                         key=lambda x: x[1], reverse=True):
    print(f"  {param}: {imp:.1%}")

# Carregar histórico
historico = pd.read_csv('resultados_*/otimizacao_bayesiana/historico_trials.csv')
print("\nResumo do Histórico:")
print(historico[['trial', 'acuracia']].describe())

```text

### Gerar Relatório Customizado

```python
from framework_investigativo_completo import consolidar_e_gerar_metadados

# Consolidar todos os resultados
resultado = consolidar_e_gerar_metadados(
    pasta_resultados='resultados_2025-12-23_14-05-56',
    verbose=True
)

print("Arquivos consolidados:", resultado['consolidacao']['status'])
print("Total de experimentos:", resultado['consolidacao']['rows_consolidated'])

```text

---


## 🎓 Checklist QUALIS A1

Após execução, verifique:

### Dados
- [ ] `resultados_completos_artigo.csv` gerado
- [ ] 8,280 experimentos registrados (Grid) ou 200+ (Bayesian)
- [ ] Metadados completos em `metadata.json`
- [ ] Seeds fixas utilizadas (42-46)


### Visualizações
- [ ] 9 figuras HTML interativas geradas
- [ ] Figuras exportadas em PNG/PDF/SVG 300 DPI
- [ ] Intervalos de confiança (IC95%) incluídos
- [ ] Labels em LaTeX (Mathjax)


### Análises Estatísticas
- [ ] ANOVA 2-way e 3-way executadas
- [ ] Effect sizes calculados (Cohen's d, Glass's Δ, Hedges' g)
- [ ] Testes post-hoc realizados (Tukey HSD, Bonferroni)
- [ ] Comparação com baselines clássicos (SVM, RF)


### Reprodutibilidade
- [ ] Código versionado (Git)
- [ ] Environment especificado (`requirements.txt`)
- [ ] Seeds documentadas
- [ ] Tempo de execução registrado
- [ ] Ambiente computacional documentado (CPU, RAM, OS)


---


## 🚀 Próximos Passos

### Para Submissão QUALIS A1

1. **Execute modo completo:**

   ```bash
   python framework_investigativo_completo.py --bayes --trials 200 --dataset-bayes all
   ```text

2. **Gere todas as figuras em alta resolução:**

   ```bash
   pip install --upgrade kaleido

   # Figuras serão exportadas automaticamente em PNG/PDF/SVG
   ```

3. **Upload no Zenodo:**
   - Crie conta em <https://zenodo.org/>
   - Faça upload de `resultados_completos_artigo.csv`
   - Faça upload da pasta `experimentos_individuais/`
   - Obtenha DOI permanente


4. **Submeta preprint no arXiv:**
   - Categoria: quant-ph (Quantum Physics)
   - Inclua link para repositório GitHub
   - Inclua DOI do Zenodo


5. **Atualize documentação:**
   - Substitua placeholders de DOI/arXiv
   - Atualize `README.md` com resultados finais
   - Verifique todas as referências


---


## 📞 Suporte

**Dúvidas sobre execução?**
1. Consulte `GUIA_EXECUCAO.md`
2. Verifique `docs/AUTOMACAO_FRAMEWORK.md`
3. Abra issue no GitHub
4. Contato: marceloclaro@gmail.com


**Problemas técnicos?**
1. Verifique `requirements.txt` está instalado
2. Confirme Python 3.9+ está instalado
3. Veja seção "Troubleshooting" em `GUIA_EXECUCAO.md`


---


## 📚 Documentação Relacionada

- [RESULTADOS_ATUALIZADOS_QUALIS_A1.md](RESULTADOS_ATUALIZADOS_QUALIS_A1.md) - Resultados da execução 23/12/2025
- [GRAFICOS_ATUALIZADOS.md](GRAFICOS_ATUALIZADOS.md) - Documentação de visualizações
- [ANALISE_QUALIS_A1.md](ANALISE_QUALIS_A1.md) - Análise de adequação para QUALIS A1
- [GUIA_EXECUCAO.md](GUIA_EXECUCAO.md) - Guia completo de execução
- [README.md](README.md) - Visão geral do projeto


---


**Última Atualização:** 23 de dezembro de 2025  
**Framework Version:** 7.2  
**Status:** Instruções validadas ✅  


🌟 *Siga essas instruções para reproduzir os resultados publicados!*
