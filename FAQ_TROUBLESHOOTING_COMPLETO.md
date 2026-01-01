# FAQ e Troubleshooting - Geração de Artigos Científicos QUALIS A1

**Versão:** 1.0  
**Data:** 26/12/2025  
**Contexto:** Framework de Geração de Artigos com Rastreabilidade Total


---


## 📋 ÍNDICE

1. [Perguntas Gerais](#perguntas-gerais)
2. [Reprodutibilidade](#reprodutibilidade)
3. [Referências e Citações](#referencias-e-citacoes)
4. [Dados e Evidências](#dados-e-evidencias)
5. [Estatística e Análise](#estatistica-e-analise)
6. [Formatação e Estilo](#formatacao-e-estilo)
7. [Problemas Técnicos](#problemas-tecnicos)
8. [Adaptação para Outras Áreas](#adaptacao-outras-areas)


---


<a name="perguntas-gerais"></a>

## 1. PERGUNTAS GERAIS

### Q1.1: Quando devo usar MODE_A vs MODE_B?

#### R:
- **MODE_A (Inglês/LaTeX):** Para submissão a periódicos internacionais (Nature, Science, Physical Review, npj QI, Quantum)
- **MODE_B (Português/ABNT):** Para submissão a periódicos brasileiros ou teses/dissertações em português


**Arquivo de configuração:** `config_artigo.json`

```json
{
  "output_mode": "MODE_A"  // ou "MODE_B"
}

```text

---


### Q1.2: Quando devo usar política R0 vs R1?

#### R:
- **R0 (Referências Travadas):** Quando a lista de referências já foi aprovada pelo orientador ou quando submetendo revisão de artigo e não pode adicionar novas referências
- **R1 (Referências Expandidas):** Durante escrita inicial, quando pode buscar e adicionar novas citações


**Fluxograma de Decisão:**

```

Posso adicionar novas referências?
├─ SIM → Use R1
└─ NÃO → Use R0 (marcar lacunas com [LACUNA DE CITAÇÃO])

```text

**Ver também:** `FLUXOGRAMA_R0_R1_COMPLETO.md`


---


### Q1.3: Quanto tempo leva cada fase?

**R:** Estimativas baseadas em projeto real:
- **Fase 1 (Auditoria):** 8-12 horas
- **Fase 2 (Bibliografia):** 6-10 horas
- **Fase 3 (Projeto):** 4-6 horas
- **Fase 4 (Redação):** 20-30 horas
- **Fase 5 (Suplementar):** 8-12 horas
- **Fase 6 (Consolidação):** 6-8 horas
- **Total:** 52-78 horas (6-10 dias úteis)


**Aceleração possível:** Com IA assistência (GPT-4/Claude), pode reduzir 30-40%.


**Ver também:** `CRONOGRAMA_ESTIMADO_COMPLETO.md`


---


<a name="reprodutibilidade"></a>

## 2. REPRODUTIBILIDADE

### Q2.1: O que fazer se o código não tem seeds fixas?

**R:**
1. **Documentar como [INFORMAÇÃO AUSENTE]**
   - No texto: "Seeds aleatórias: [INFORMAÇÃO AUSENTE] - código não especifica valores fixos"

   
2. **Adicionar em Threats to Validity:**

   ```markdown

   **Ameaça à Reprodutibilidade Estocástica:** O código original não utiliza

   seeds fixas, resultando em variabilidade entre execuções. Para mitigar,
   executamos cada configuração N=10 vezes e reportamos média ± desvio padrão.
   ```text

3. **Executar múltiplas vezes:**

   ```python
   results = []
   for run in range(10):

       # Não definir seed - deixar aleatório
       result = execute_experiment()
       results.append(result)
   
   mean = np.mean(results)
   std = np.std(results)
   print(f"Acurácia: {mean:.2%} ± {std:.2%}")
   ```text

4. **Reportar variabilidade:**
   - "Acurácia média: 63.5% ± 2.1% (N=10 execuções, seeds aleatórias)"


**Importante:** Sempre priorize adicionar seeds fixas ao código quando possível!


---


### Q2.2: Como proceder se não há logs de execução?

**R:**


**Opção A - Executar e Gerar Logs (RECOMENDADO):**

```bash

# Executar pipeline completo com logging
python framework_investigativo_completo.py --log-file resultados.log \

  --save-config config_exec.json \
  --output-dir resultados_$(date +%Y%m%d_%H%M%S)

```text

**Opção B - Execução Inviável (limitações computacionais):**
1. Marcar resultados como **[NÃO DISPONÍVEL]**
2. Focar o artigo na **metodologia** e **plano experimental**
3. Adicionar seção "Experimental Design" detalhada
4. Propor análise teórica/simulação em menor escala


**Exemplo de texto:**

```markdown
Devido a limitações computacionais (tempo de execução estimado:
720 horas em hardware disponível), os resultados experimentais completos
estão [NÃO DISPONÍVEL]. Este trabalho apresenta:

1. Framework metodológico completo e validado em subset reduzido (N=100)
2. Análise teórica fundamentada em trabalhos anteriores
3. Plano experimental detalhado para execução futura

```text

---


### Q2.3: Como documentar o ambiente de execução?

**R:** Criar arquivo `ENVIRONMENT.md`:


```markdown

# Ambiente de Execução

## Hardware
- **CPU:** Intel Xeon E5-2680 v4 @ 2.40GHz (28 cores)
- **RAM:** 128 GB DDR4
- **GPU:** NVIDIA Tesla V100 32GB (não utilizada neste projeto)
- **Armazenamento:** 2TB SSD NVMe


## Software
- **OS:** Ubuntu 22.04 LTS (kernel 5.15.0)
- **Python:** 3.9.18
- **PennyLane:** 0.38.0
- **Qiskit:** 1.0.2
- **Cirq:** 1.4.0
- **NumPy:** 1.24.3
- **SciPy:** 1.10.1
- **Matplotlib:** 3.7.1
- **Scikit-learn:** 1.3.0


## Instalação

```bash
pip install -r requirements.txt

```text

## Reprodução

```bash
python framework_investigativo_completo.py --seed 42

```text

## Tempo de Execução
- **Configuração completa:** 48-72 horas
- **Modo rápido (--bayes --trials 100):** 1-2 horas

```

**Automatização:**

```python
import platform
import sys

def generate_environment_info():
    info = {
        "os": platform.system(),
        "python": sys.version,
        "packages": {}
    }

    # Coletar versões de pacotes
    import pkg_resources
    for pkg in pkg_resources.working_set:
        info["packages"][pkg.key] = pkg.version
    return info

```text

---


<a name="referencias-e-citacoes"></a>

## 3. REFERÊNCIAS E CITAÇÕES

### Q3.1: Quando usar [INFORMAÇÃO AUSENTE] vs [NÃO DISPONÍVEL] vs [LACUNA DE CITAÇÃO]?

**R:**


| Marcador | Quando Usar | Exemplo |
|----------|-------------|---------|
| **[INFORMAÇÃO AUSENTE]** | Info deveria existir mas não foi encontrada | "Versão do PyTorch: [INFORMAÇÃO AUSENTE]" |
| **[NÃO DISPONÍVEL]** | Info não pode ser gerada/obtida | "Resultados em hardware IBM real: [NÃO DISPONÍVEL]" |
| **[LACUNA DE CITAÇÃO]** | Falta referência (apenas em R0) | "Quantum supremacy [LACUNA DE CITAÇÃO] foi demonstrado" |

**Fluxograma:**

```

A informação existe no código/logs?
├─ SIM → Extrair e reportar
├─ NÃO (mas deveria existir) → [INFORMAÇÃO AUSENTE]
├─ NÃO (e não pode obter) → [NÃO DISPONÍVEL]
└─ Falta citação (modo R0) → [LACUNA DE CITAÇÃO]

```text

---


### Q3.2: Como buscar referências em modo R1?

**R:** Seguir as **7 categorias** obrigatórias:


#### 1. Fundamentos Teóricos
- Nielsen & Chuang (2010) - Quantum Computing
- Preskill (2018) - NISQ era
- Teoria de informação quântica


**Onde buscar:** Livros-texto clássicos, reviews em Rev. Mod. Phys.


#### 2. Estado da Arte
- Últimos 5 anos (2019-2024)
- Trabalhos em periódicos QUALIS A1


**Onde buscar:** Google Scholar, arXiv (quant-ph, cs.LG)


#### 3. Metodologia
- Referências que justificam escolhas metodológicas
- Protocolos, testes estatísticos


**Onde buscar:** Methodological papers, software documentation


#### 4. Benchmarks
- Datasets padrão (Iris, Wine, Breast Cancer)
- Baselines estabelecidos


**Onde buscar:** UCI ML Repository, Kaggle, papers de benchmark


#### 5. Frameworks e Ferramentas
- PennyLane, Qiskit, Cirq documentation
- Software libraries


**Onde buscar:** Official docs, JOSS (Journal of Open Source Software)


#### 6. Aplicações
- Casos de uso práticos
- Estudos de viabilidade


**Onde buscar:** Application-focused journals


#### 7. Surveys e Reviews
- Revisões abrangentes
- Taxonomias


**Onde buscar:** ACM Computing Surveys, IEEE Access


**Ferramenta automatizada:**

```python
from scholarly import scholarly

def buscar_referencias(termo, num_resultados=10):
    search = scholarly.search_pubs(termo)
    refs = []
    for i, pub in enumerate(search):
        if i >= num_resultados:
            break
        refs.append({
            "title": pub['bib']['title'],
            "author": pub['bib']['author'],
            "year": pub['bib'].get('pub_year', 'N/A'),
            "doi": pub.get('doi', 'N/A')
        })
    return refs

# Exemplo
refs = buscar_referencias("beneficial quantum noise variational")

```text

---


### Q3.3: Como formatar referências ABNT?

**R:** Template padrão:


**Artigo de periódico:**

```

AUTOR, Nome. Título do artigo. **Nome da Revista**, v. X, n. Y, p. Z-W, ano.
DOI: <https://doi.org/10.xxxx/yyyyy>

```text

**Exemplo real:**

```

DU, Yuxuan et al. Quantum noise can help quantum sensing. **Physical Review Letters**,
v. 128, n. 8, p. 080506, 2021. DOI: <https://doi.org/10.1103/PhysRevLett.128.080506>

```text

**Livro:**

```

AUTOR, Nome. **Título do livro**. Edição. Cidade: Editora, ano.

```text

**Exemplo:**

```

NIELSEN, Michael A.; CHUANG, Isaac L. **Quantum computation and quantum information**.
10th ed. Cambridge: Cambridge University Press, 2010.

```text

**Preprint arXiv:**

```

AUTOR, Nome. Título. **arXiv preprint** arXiv:XXXX.YYYYY, ano.

```text

**Ferramenta:** Use BibTeX + `bibtex2abnt.py` para conversão automática


---


<a name="dados-e-evidencias"></a>

## 4. DADOS E EVIDÊNCIAS

### Q4.1: Como calcular o total de configurações?

**R:**


**Fórmula Geral:**

```

Total = Datasets × Ansätze × Ruídos × Schedules × Inicializações × Seeds × Outras_dims

```text

**Exemplo deste projeto:**

```

Total = 4 × 7 × 6 × 2 × 8 × 2 × 1 = 2.688 configurações (grid search)

```text

Ou com otimização Bayesiana:

```

Total_executado = 8.280 configurações (100 trials × ~83 combinações relevantes)

```text

**Código para verificação:**

```python
import itertools

configs = {
    "datasets": ["iris", "wine", "breast_cancer", "moons"],
    "ansatze": ["basic", "strongly_entangling", "random", ...],  # 7 total
    "noises": ["depolarizing", "phase_damping", ...],  # 6 total
    "schedules": ["static", "cosine"],  # 2 total
    "inits": list(range(8)),  # 8 total
    "seeds": [42, 43]  # 2 total
}

# Cartesian product
all_configs = list(itertools.product(
    configs["datasets"],
    configs["ansatze"],
    configs["noises"],
    configs["schedules"],
    configs["inits"],
    configs["seeds"]
))

print(f"Total de configurações: {len(all_configs)}")

# Output: Total de configurações: 2688

```text

**Registrar em:** `fase1_analise/analise_codigo_inicial.md`


---


### Q4.2: Como validar conivência (código-texto)?

**R:** Processo em 3 etapas:


#### Etapa 1: Extrair Valores do Código

```python
import re

def extract_values_from_code(file_path):
    with open(file_path, 'r') as f:
        code = f.read()
    
    # Extrair constantes
    constants = re.findall(r'(\w+)\s*=\s*([0-9.]+)', code)
    return dict(constants)

values_code = extract_values_from_code("framework_investigativo_completo.py")

```text

#### Etapa 2: Extrair Valores do Texto

```python
def extract_values_from_text(md_file):
    with open(md_file, 'r') as f:
        text = f.read()
    
    # Extrair valores numéricos com contexto
    values = re.findall(r'(\w+)[:\s]+([0-9.]+)%?', text)
    return dict(values)

values_text = extract_values_from_text("artigo_cientifico/fase4_secoes/resultados_completo.md")

```text

#### Etapa 3: Comparar

```python
def verificar_conivencia(values_code, values_text):
    discrepancias = []
    
    for key in values_text:
        if key in values_code:
            if float(values_code[key]) != float(values_text[key]):
                discrepancias.append({
                    "variavel": key,
                    "codigo": values_code[key],
                    "texto": values_text[key]
                })
    
    conivencia = 100 * (1 - len(discrepancias) / len(values_text))
    return conivencia, discrepancias

conivencia, discrep = verificar_conivencia(values_code, values_text)
print(f"Conivência: {conivencia:.1f}%")

```text

**Gerar relatório:** `fase6_consolidacao/relatorio_conivencia.md`


---


### Q4.3: Como criar tabela de rastreabilidade?

**R:** Template obrigatório:


```markdown
| Seção | Afirmação/Número | Evidência (Arquivo:Linha) | Referência |
|-------|------------------|---------------------------|------------|
| 4.5 Results | Acurácia 65.83% | `resultados.csv:row_1523:col_acc` | - |
| 4.4 Methods | Eq. Lindblad | `noise.py:L15` | (Lindblad, 1976) |
| 4.6 Discussion | Ruído como regularizador | `sintese_literatura.md` | (Du et al., 2021) |

```text

**Automatização:**

```python
import pandas as pd

def create_traceability_table():
    table = []
    
    # Extrair de resultados
    df = pd.read_csv("resultados.csv")
    max_acc = df["accuracy"].max()
    idx = df["accuracy"].idxmax()
    
    table.append({
        "Seção": "4.5 Results",
        "Afirmação": f"Acurácia máxima {max_acc:.2%}",
        "Evidência": f"resultados.csv:row_{idx}:col_accuracy",
        "Referência": "-"
    })
    
    # Extrair de código
    # ... (similar para outros valores)
    
    df_table = pd.DataFrame(table)
    df_table.to_markdown("rastreabilidade_completa.md", index=False)

```text

**Ver também:** `templates/rastreabilidade_completa_template.md`


---


<a name="estatistica-e-analise"></a>

## 5. ESTATÍSTICA E ANÁLISE

### Q5.1: Que testes estatísticos devo usar?

**R:** Depende do tipo de comparação:


#### Comparação 2 Grupos (Com ruído vs Sem ruído)

```python
from scipy.stats import ttest_ind, mannwhitneyu

# Dados normalmente distribuídos → t-test
t_stat, p_value = ttest_ind(com_ruido, sem_ruido)

# Dados não-normais → Mann-Whitney U
u_stat, p_value = mannwhitneyu(com_ruido, sem_ruido, alternative='greater')

```text

#### Comparação 3+ Grupos (Múltiplos tipos de ruído)

```python
from scipy.stats import f_oneway, kruskal

# Dados normalmente distribuídos → ANOVA
f_stat, p_value = f_oneway(depolarizing, phase_damping, amplitude_damping)

# Dados não-normais → Kruskal-Wallis
h_stat, p_value = kruskal(depolarizing, phase_damping, amplitude_damping)

```text

#### Análise Multifatorial (Interações)

```python
import statsmodels.api as sm
from statsmodels.formula.api import ols

# ANOVA 3-way (ruído × ansatz × dataset)
model = ols('accuracy ~ C(noise_type) * C(ansatz) * C(dataset)', data=df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)

```text

#### Post-hoc (após ANOVA significativa)

```python
from scipy.stats import tukey_hsd

# Tukey HSD
res = tukey_hsd(depolarizing, phase_damping, amplitude_damping)
print(res)

```text

**Importante:** Sempre reportar:
1. Estatística de teste
2. p-valor
3. Tamanho de efeito (Cohen's d, η², etc.)
4. Intervalo de confiança


---


### Q5.2: Como calcular tamanhos de efeito?

**R:**


#### Cohen's d (diferença entre médias)

```python
import numpy as np

def cohens_d(group1, group2):
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    
    # Pooled standard deviation
    pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    
    d = (np.mean(group1) - np.mean(group2)) / pooled_std
    return d

d = cohens_d(com_ruido, sem_ruido)
print(f"Cohen's d = {d:.2f}")

# Interpretação
if abs(d) < 0.2:
    print("Efeito pequeno")
elif abs(d) < 0.5:
    print("Efeito médio")
elif abs(d) < 0.8:
    print("Efeito grande")
else:
    print("Efeito muito grande")

```text

#### Eta-squared (η² - para ANOVA)

```python
def eta_squared(anova_table):
    ss_between = anova_table["sum_sq"][0]  # Variância entre grupos
    ss_total = anova_table["sum_sq"].sum()  # Variância total
    
    eta_sq = ss_between / ss_total
    return eta_sq

eta_sq = eta_squared(anova_table)
print(f"η² = {eta_sq:.3f}")  # 0.01 (pequeno), 0.06 (médio), 0.14 (grande)

```text

#### Cliff's Delta (alternativa não-paramétrica)

```python
def cliffs_delta(group1, group2):
    n1, n2 = len(group1), len(group2)
    
    # Contar pares onde group1 > group2
    dominance = sum(x > y for x in group1 for y in group2)
    
    delta = (dominance - (n1*n2 - dominance)) / (n1 * n2)
    return delta

delta = cliffs_delta(com_ruido, sem_ruido)
print(f"Cliff's Δ = {delta:.2f}")  # |Δ| < 0.147 (negligível), < 0.33 (pequeno), < 0.474 (médio), ≥ 0.474 (grande)

```text

**Reportar sempre:**

```markdown
Phase Damping superou Depolarizing em +3.75 pontos percentuais
(t(98) = 4.12, p < 0.001, Cohen's d = 0.83 [IC 95%: 0.41, 1.25], grande).

```text

---


### Q5.3: Como corrigir para múltiplas comparações?

**R:**


#### Bonferroni Correction (conservadora)

```python
from scipy.stats import ttest_ind

p_values = []
comparisons = [
    ("depolarizing", "phase_damping"),
    ("depolarizing", "amplitude_damping"),
    ("phase_damping", "amplitude_damping")
]

for group1_name, group2_name in comparisons:
    t_stat, p = ttest_ind(data[group1_name], data[group2_name])
    p_values.append(p)

# Bonferroni
alpha = 0.05
bonferroni_alpha = alpha / len(p_values)
print(f"α corrigido: {bonferroni_alpha:.4f}")

significant = [p < bonferroni_alpha for p in p_values]

```text

#### Benjamini-Hochberg (FDR - menos conservadora)

```python
from statsmodels.stats.multitest import multipletests

# FDR correction
reject, p_adjusted, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

for i, (comp, p_raw, p_adj, sig) in enumerate(zip(comparisons, p_values, p_adjusted, reject)):
    print(f"{comp}: p_raw={p_raw:.4f}, p_adj={p_adj:.4f}, sig={sig}")

```text

#### Escolha do método:
- Bonferroni: Quando poucas comparações (≤10) ou custo de erro tipo I alto
- Benjamini-Hochberg: Quando muitas comparações (>10) ou exploratório


---


<a name="formatacao-e-estilo"></a>

## 6. FORMATAÇÃO E ESTILO

### Q6.1: Quantas palavras deve ter cada seção?

**R:** Guideline para periódicos QUALIS A1:


| Seção | Palavras | Caracteres | Parágrafos |
|-------|----------|------------|------------|
| **Abstract** | 250-300 | 1.500-2.000 | 1 (estruturado IMRAD) |
| **Introduction** | 1.500-2.500 | 9.000-15.000 | 8-12 |
| **Related Work** | 2.000-3.000 | 12.000-18.000 | 10-15 |
| **Methods** | 2.500-4.000 | 15.000-24.000 | 12-20 |
| **Results** | 2.000-3.000 | 12.000-18.000 | 10-15 |
| **Discussion** | 2.500-4.000 | 15.000-24.000 | 12-20 |
| **Conclusion** | 500-800 | 3.000-5.000 | 4-6 |
| **Total** | 12.000-18.000 | 72.000-108.000 | 60-100 |

**Verificar contagem:**

```bash

# Palavras
wc -w introducao_completa.md

# Caracteres
wc -c introducao_completa.md

# Parágrafos (linhas vazias + 1)
grep -c '^$' introducao_completa.md

```text

---


### Q6.2: Como estruturar parágrafos para QUALIS A1?

**R:**


**Anatomia do Parágrafo Ideal (5-8 frases):**


1. **Topic sentence:** Ideia principal
2. **Supporting evidence:** Dados/citações
3. **Analysis:** Interpretação
4. **Connection:** Link com próxima ideia


**Exemplo Ruim (muito curto):**

```markdown
Ruído quântico degrada fidelidade. Isso é um problema.

```text

**Exemplo Bom:**

```markdown
Ruído quântico degrada a fidelidade de operações em dispositivos NISQ,
representando o principal obstáculo para vantagem quântica prática [1].
Taxas de erro típicas de 10⁻³ por porta resultam em fidelidade exponencialmente
decrescente com profundidade do circuito [2]. No entanto, trabalhos recentes
sugerem que ruído controlado pode agir como regularizador, análogo ao dropout
em redes neurais clássicas [3,4]. Esta mudança de paradigma — de obstáculo a
oportunidade — motiva nossa investigação sistemática dos regimes benéficos de ruído.

```text

**Transições entre parágrafos:**

```markdown
<!-- Fim parágrafo 1 -->
... resultando em melhor generalização.

<!-- Transição -->
No entanto, esta abordagem apresenta limitações quando...

<!-- Início parágrafo 2 -->

```text

---


### Q6.3: Como formatar equações LaTeX?

**R:**


**Inline (dentro do texto):**

```markdown
A fidelidade é dada por $F = \langle\psi|\rho|\psi\rangle$.

```text

**Display (centralizada, numerada):**

```markdown
A evolução do sistema quântico aberto segue a equação mestra de Lindblad:

$$
\frac{d\rho}{dt} = -i[H, \rho] + \sum_k \gamma_k \left(
  L_k \rho L_k^\dagger - \frac{1}{2}\{L_k^\dagger L_k, \rho\}
\right) \tag{1}
$$

onde $\gamma_k$ são taxas de dissipação e $L_k$ são operadores de Lindblad [1].

```text

**Alinhamento multi-linha:**

```markdown
$$
\begin{align}
\mathcal{L}(\theta) &= \sum_{i=1}^N \ell(y_i, f(x_i; \theta)) \tag{2a} \\
&= -\sum_{i=1}^N y_i \log p_i + (1-y_i)\log(1-p_i) \tag{2b}
\end{align}
$$

```text

**Sempre incluir:**
1. Número da equação (tag)
2. Parágrafo explicativo após a equação
3. Definição de todas as variáveis


**Exemplo completo:**

```markdown
O circuito variacional é definido por:

$$
U(\theta) = \prod_{l=1}^L U_l(\theta_l) \tag{3}
$$

onde $U_l$ representa a camada $l$, $\theta_l \in \mathbb{R}^{P_l}$ são parâmetros
treináveis da camada, $L$ é a profundidade total, e $P_l$ é a dimensionalidade
dos parâmetros em cada camada. Para nossa implementação, utilizamos $L=2$ camadas
e $P_l=16$ parâmetros por camada, totalizando 32 parâmetros variacionais.

```text

---


<a name="problemas-tecnicos"></a>

## 7. PROBLEMAS TÉCNICOS

### Q7.1: Erro ao importar bibliotecas

**R:**


**Problema:**

```python
ModuleNotFoundError: No module named 'pennylane'

```text

**Solução:**

```bash

# Instalar bibliotecas
pip install pennylane qiskit cirq numpy scipy scikit-learn matplotlib

# Verificar instalação
python -c "import pennylane as qml; print(qml.__version__)"

# Criar requirements.txt
pip freeze > requirements.txt

```text

**requirements.txt mínimo:**

```

pennylane==0.38.0
qiskit==1.0.2
cirq==1.4.0
numpy>=1.24.0
scipy>=1.10.0
scikit-learn>=1.3.0
matplotlib>=3.7.0
pandas>=2.0.0
seaborn>=0.12.0

```text

---


### Q7.2: Execução muito lenta

**R:**


**Diagnóstico:**

```python
import time

start = time.time()

# ... seu código ...
end = time.time()

print(f"Tempo: {end-start:.2f} segundos")

```text

**Otimizações:**


#### 1. Reduzir configurações (Bayesian Optimization)

```python

# Em vez de grid search completo (2.688 configs)
from skopt import gp_minimize

def objective(params):
    gamma, lr, depth = params

    # ... executar experimento ...
    return -accuracy  # Minimizar negativo = maximizar acurácia

# Buscar apenas 100 configurações promissoras
result = gp_minimize(
    objective,
    dimensions=[(0, 0.01), (0.001, 0.1), (1, 5)],
    n_calls=100,
    random_state=42
)

```text

#### 2. Paralelização

```python
from joblib import Parallel, delayed

results = Parallel(n_jobs=4)(
    delayed(run_experiment)(config)
    for config in configurations
)

```text

#### 3. Caching de resultados

```python
import pickle

def run_with_cache(config, cache_file="cache.pkl"):

    # Tentar carregar do cache
    try:
        with open(cache_file, 'rb') as f:
            cache = pickle.load(f)
        if config in cache:
            return cache[config]
    except FileNotFoundError:
        cache = {}
    
    # Executar experimento
    result = run_experiment(config)
    
    # Salvar no cache
    cache[config] = result
    with open(cache_file, 'wb') as f:
        pickle.dump(cache, f)
    
    return result

```text

---


### Q7.3: Resultados inconsistentes entre execuções

**R:**


**Problema:** Mesma config, resultados diferentes


**Causa:** Seeds não fixas


**Solução:**

```python
import random
import numpy as np
import pennylane as qml

def set_all_seeds(seed=42):
    """Fixar todas as seeds para reprodutibilidade"""
    random.seed(seed)
    np.random.seed(seed)

    # PennyLane usa NumPy, então np.random.seed suficiente
    
    # Para PyTorch (se usado)
    try:
        import torch
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass
    
    # Para TensorFlow (se usado)
    try:
        import tensorflow as tf
        tf.random.set_seed(seed)
    except ImportError:
        pass

# Usar no início do script
set_all_seeds(42)

```text

**Validação:**

```python

# Executar 3 vezes com mesma seed
for i in range(3):
    set_all_seeds(42)
    result = run_experiment(config)
    print(f"Run {i+1}: {result}")

# Deve imprimir valores idênticos

```text

---


<a name="adaptacao-outras-areas"></a>

## 8. ADAPTAÇÃO PARA OUTRAS ÁREAS

### Q8.1: Como adaptar para Machine Learning Clássico?

**R:**


**Substituições:**


| Quantum ML | Classical ML |
|------------|--------------|
| Ansatz | Network architecture (CNN, RNN, Transformer) |
| Quantum noise | Dropout, Data augmentation |
| Qubits | Hidden units, Layers |
| Phase Damping | Specific dropout variants |
| γ (noise intensity) | dropout_rate |
| Lindblad equation | Backpropagation |

**Exemplo adaptado:**

```python

# Original (Quantum)
qml.RY(params[0], wires=0)
qml.DepolarizingChannel(gamma, wires=0)

# Adaptado (Classical)
x = Dense(64, activation='relu')(x)
x = Dropout(dropout_rate)(x, training=True)  # training=True para inferência com dropout

```text

---


### Q8.2: Como adaptar para Bioinformática?

**R:**


**Substituições:**


| Quantum Context | Bioinformatics Context |
|----------------|----------------------|
| Circuit depth | Sequence length |
| Qubits | Genes, Proteins |
| Entanglement | Gene co-expression, Protein-protein interaction |
| Noise types | Sequencing errors, Batch effects |
| VQC | Classification (disease vs healthy) |

**Exemplo de Hipótese Adaptada:**

```markdown

**H₀ (Original):**

Existe um γ > 0 tal que ruído quântico melhora VQC.

**H₀ (Adaptado para Bio):**

Existe um nível de ruído de sequenciamento onde algoritmos de classificação
de doenças apresentam melhor generalização (via regularização estocástica).

```

---


### Q8.3: Como adaptar para Ciências Sociais?

**R:**


**Substituições:**


| Quantum/CS | Social Sciences |
|------------|----------------|
| Experimental configurations | Survey conditions, Treatment groups |
| Accuracy metric | Response rate, Effect size |
| Noise | Measurement error, Response bias |
| Hyperparameters | Interview protocols, Questionnaire versions |
| Code-text congruence | Data-claims correspondence |

**Adaptações metodológicas:**


1. **Reprodutibilidade:**
   - Seeds fixas → Pre-registration (osf.io)
   - Version control → Protocol versioning
   - Logs → Field notes


2. **Estatística:**
   - ANOVA → ANCOVA (covariates)
   - Effect sizes → Same (Cohen's d, η²)
   - Multiple comparisons → Same (Bonferroni, FDR)


3. **Auditoria:**
   - Código → Interview transcripts
   - Dados → Survey responses (anonymized)
   - Rastreabilidade → Chain of evidence


**Exemplo de Tabela Rastreabilidade Adaptada:**


| Seção | Afirmação | Evidência | Referência |
|-------|-----------|-----------|------------|
| Results | 67% concordam | `survey_data.csv:Q15:row_1-500` | - |
| Methods | Escala Likert 5 pontos | `protocol_v2.pdf:p.8` | (Likert, 1932) |
| Discussion | Efeito moderado | Calculado de survey_data.csv | (Cohen, 1988) |

---


## 🆘 HELP! Meu Problema Não Está Aqui

### Procedimento:

1. **Verificar documentação principal:**
   - `README.md`
   - `GUIA_COMPLETO_GERACAO_ARTIGOS.md`
   - `artigo_cientifico/README.md`


2. **Verificar templates:**
   - `templates/` directory
   - `artigo_cientifico/fase*/` directories


3. **Consultar exemplos:**
   - Arquivos `exemplo_*` no repositório
   - `artigo_cientifico/fase4_secoes/` (seções completas como referência)


4. **Criar issue no GitHub:**
   - Descrever problema em detalhe
   - Incluir mensagem de erro completa
   - Anexar código reproduzível mínimo


5. **Contato direto:**
   - Email: [inserir email do mantenedor]
   - Twitter/X: [inserir handle]


---


## 📚 Referências Úteis

1. **Metodologia Científica:**
   - Creswell, J. W. (2014). Research Design. 4th ed. SAGE.
   - Field, A. (2013). Discovering Statistics Using IBM SPSS. 4th ed. SAGE.


2. **Escrita Acadêmica:**
   - Swales, J. M. (1990). Genre Analysis. Cambridge University Press.
   - Belcher, W. L. (2019). Writing Your Journal Article in Twelve Weeks. 2nd ed. Chicago.


3. **Reprodutibilidade:**
   - Wilkinson et al. (2016). The FAIR Guiding Principles. Scientific Data.
   - Baker, M. (2016). 1,500 scientists lift the lid on reproducibility. Nature.


4. **Quantum ML:**
   - Cerezo et al. (2021). Variational Quantum Algorithms. Nature Reviews Physics.
   - Benedetti et al. (2019). Quantum Machine Learning. PRX Quantum.


---


**Última atualização:** 26/12/2025  
**Versão:** 1.0  
**Mantenedor:** Framework Geração Artigos QUALIS A1  
**Status:** ✅ Ativo e mantido

