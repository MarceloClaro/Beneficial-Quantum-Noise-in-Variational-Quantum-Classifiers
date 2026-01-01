# Automação do Framework Investigativo v7.2

## 📋 Visão Geral

O Framework Investigativo agora inclui **consolidação e orquestração automática integrada**, eliminando a necessidade de scripts externos para pós-processamento de resultados.

## 🆕 Novos Recursos (v7.2)

### 1. Consolidação Automática
- **Função**: `consolidar_resultados_individuais(pasta_resultados, verbose=True)`
- **O que faz**:
  - Localiza e lê todos os CSVs na pasta `experimentos_individuais/`
  - Consolida em um único arquivo `resultados_completos_artigo.csv`
  - Retorna sumário com estatísticas (linhas, colunas, status)


### 2. Comparação de Baselines
- **Função**: `gerar_comparacao_baselines(df_all, pasta_resultados, verbose=True)`
- **O que faz**:
  - Compara o melhor resultado VQC vs. SVM e Random Forest por dataset
  - Calcula deltas (vantagem/desvantagem do VQC)
  - Salva em `comparacao_baselines.csv`


### 3. Geração de Metadados
- **Função**: `gerar_metadata_orchestrator(pasta_resultados, consolidacao_info, verbose=True)`
- **O que faz**:
  - Cria inventário completo de arquivos e subpastas
  - Registra estatísticas da consolidação
  - Salva em `metadata_orchestrator.json`


### 4. Função Unificada
- **Função**: `consolidar_e_gerar_metadados(pasta_resultados, verbose=True)`
- **O que faz**:
  - Executa todas as etapas acima em sequência
  - Chamada automaticamente ao final do `main()`
  - Retorna dict com sumário completo das operações


## 🚀 Modos de Execução

### Execução Completa (Grid Search)

```bash
python framework_investigativo_completo.py

```text

- Executa grid search completo (8,280 configurações × 5 seeds)
- Gera todas as visualizações e análises
- **NOVO**: Consolida automaticamente os resultados ao final


### Modo Rápido

```bash

# Windows PowerShell
$env:VQC_QUICK="1"; python framework_investigativo_completo.py

# Linux/Mac
VQC_QUICK=1 python framework_investigativo_completo.py

```text

- Reduz épocas para 5 (vs. 15 padrão)
- Útil para testes e validação rápida


### Otimização Bayesiana

```bash

# Windows PowerShell
$env:VQC_BAYESIAN="1"; python framework_investigativo_completo.py

# Linux/Mac
VQC_BAYESIAN=1 python framework_investigativo_completo.py

```text

- Usa Optuna (TPE) para busca inteligente
- 100-200 trials (vs. 8,280 do grid)
- 10-20x mais rápido que grid completo


### Bayesiano com Parâmetros Customizados

```bash
python framework_investigativo_completo.py --bayes --trials 150 --epocas-bayes 20 --dataset-bayes all

```text

- `--trials`: número de trials Optuna (padrão: 100-200)
- `--epocas-bayes`: épocas por trial (padrão: 15)
- `--dataset-bayes`: dataset alvo (`moons`, `circles`, `iris`, `breast_cancer`, `wine`, `all`)


### Grid + Bayesiano (Híbrido)

```bash
python framework_investigativo_completo.py --bayes-after-grid

# ou
python framework_investigativo_completo.py --run-both

```text

- Executa grid search completo primeiro
- Depois refina com otimização Bayesiana
- Melhor para artigos Qualis A1 (exploração + refinamento)


## 📂 Estrutura de Resultados

Após a execução, o diretório `resultados_YYYY-MM-DD_HH-MM-SS/` contém:

```

resultados_2025-10-28_17-09-33/
├── README.md                              # Descrição do experimento
├── metadata.json                          # Metadados raiz (gerado pelo main)
├── metadata_orchestrator.json             # ✨ NOVO: Metadados de consolidação
├── resultados_completos_artigo.csv        # ✨ Consolidado automaticamente
├── comparacao_baselines.csv               # ✨ NOVO: VQC vs. SVM/RF
├── figura2_beneficial_noise.html          # Visualizações interativas
├── figura2b_beneficial_noise_ic95.html    # (com intervalos de confiança)
├── figura3_noise_types.html
├── figura3b_noise_types_ic95.html
├── figura4_initialization.html
├── figura5_architecture_tradeoffs.html
├── figura6_effect_sizes.html
├── figura7_overfitting.html
├── experimentos_individuais/              # CSVs de cada configuração
│   ├── exp_00001.csv
│   ├── exp_00002.csv
│   └── ...
├── circuitos/                             # Circuitos exportados (PNG/QASM)
├── barren_plateaus/                       # Análises de Barren Plateaus
├── analises_individuais/                  # Análises estatísticas (ANOVA, etc.)
└── otimizacao_bayesiana/                  # Resultados Optuna (se --bayes)
    ├── resultado_otimizacao.json
    ├── historico_trials.csv
    └── README_otimizacao.md

```text

## 📊 Formato dos Arquivos Gerados

### resultados_completos_artigo.csv
Consolidação de todos os experimentos individuais:

```csv
dataset,arquitetura,estrategia_init,tipo_ruido,nivel_ruido,n_qubits,n_camadas,acuracia_treino,acuracia_teste,gap_treino_teste,tempo_segundos,custo_final,cm_tn,cm_fp,cm_fn,cm_tp,seed
moons,basic_entangler,matematico,sem_ruido,0.0,4,2,0.585714,0.550000,0.035714,83.41,0.929987,44,16,38,22,42
...

```text

### comparacao_baselines.csv
Comparação por dataset:

```csv
dataset,vqc_melhor,vqc_sem_ruido_media,svm,rf,delta_vqc_svm,delta_vqc_rf
moons,0.9250,0.8800,0.8950,0.9100,0.0300,-0.0150
circles,0.9500,0.9200,0.9350,0.9400,0.0150,0.0100
...

```text

### metadata_orchestrator.json
Inventário e estatísticas:

```json
{
  "tipo": "metadata_orchestrator",
  "versao_framework": "7.2",
  "timestamp": "2025-10-28 17:45:32",
  "pasta_resultados": "c:\\...\\resultados_2025-10-28_17-09-33",
  "arquivos_raiz": [
    "comparacao_baselines.csv",
    "metadata_orchestrator.json",
    "resultados_completos_artigo.csv",
    ...
  ],
  "subpastas": [
    "analises_individuais",
    "barren_plateaus",
    "circuitos",
    "experimentos_individuais",
    "otimizacao_bayesiana"
  ],
  "consolidacao": {
    "status": "ok",
    "num_csvs_individuais": 71,
    "rows_consolidated": 71,
    "columns": ["dataset", "arquitetura", ...]
  }
}

```text

## 🛠️ Uso Programático

### Consolidar Resultados Existentes

```python
from framework_investigativo_completo import consolidar_e_gerar_metadados

# Pós-processar um diretório existente
resultado = consolidar_e_gerar_metadados(
    pasta_resultados="resultados_2025-10-28_17-09-33",
    verbose=True
)

print(f"Consolidados: {resultado['consolidacao']['rows_consolidated']} experimentos")
print(f"Metadados em: {resultado['metadata_orchestrator']}")

```text

### Apenas Consolidação

```python
from framework_investigativo_completo import consolidar_resultados_individuais

info = consolidar_resultados_individuais(
    pasta_resultados="resultados_2025-10-28_17-09-33",
    verbose=True
)

if info['status'] == 'ok':
    print(f"CSV consolidado: {info['consolidated_path']}")
    print(f"Total de linhas: {info['rows_consolidated']}")

```text

### Apenas Comparação de Baselines

```python
import pandas as pd
from framework_investigativo_completo import gerar_comparacao_baselines

df = pd.read_csv("resultados_2025-10-28_17-09-33/resultados_completos_artigo.csv")
comp_path = gerar_comparacao_baselines(
    df_all=df,
    pasta_resultados="resultados_2025-10-28_17-09-33",
    verbose=True
)

print(f"Comparação salva em: {comp_path}")

```text

## 🔄 Fluxo de Execução

```

┌─────────────────────────────────────────────────────────────────┐
│ 1. Carregar Datasets (moons, circles, iris, breast_cancer, wine)│
└────────────────────────┬────────────────────────────────────────┘
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│ 2. Executar Grid Search / Bayesiano / Ambos                     │
│    - Salva CSVs individuais em experimentos_individuais/        │
│    - Salva circuitos em circuitos/                              │
│    - Salva análises de Barren Plateaus em barren_plateaus/      │
└────────────────────────┬────────────────────────────────────────┘
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│ 3. Análises Estatísticas (ANOVA, post-hoc, effect sizes)        │
│    - Salva em analises_individuais/                             │
└────────────────────────┬────────────────────────────────────────┘
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│ 4. Gerar Visualizações (9 figuras interativas HTML)             │
└────────────────────────┬────────────────────────────────────────┘
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│ 5. Análises Profundas (PCA, clustering, correlação)             │
└────────────────────────┬────────────────────────────────────────┘
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│ 6. Resumo Final (estatísticas, melhores configs, ruídos benéficos)│
└────────────────────────┬────────────────────────────────────────┘
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│ 7. ✨ CONSOLIDAÇÃO AUTOMÁTICA (v7.2)                            │
│    - Consolida experimentos_individuais/*.csv                   │
│    - Gera resultados_completos_artigo.csv                       │
│    - Gera comparacao_baselines.csv                              │
│    - Gera metadata_orchestrator.json                            │
└─────────────────────────────────────────────────────────────────┘

```text

## ⚡ Comparação de Performance

| Modo | Configurações | Tempo Estimado | Uso |
|------|--------------|----------------|-----|
| **Grid Completo** | 8,280 × 5 seeds = 41,400 | ~48-72h | Publicações Qualis A1 |
| **Grid Rápido** (VQC_QUICK=1) | 8,280 × 5 seeds = 41,400 | ~16-24h | Testes extensivos |
| **Bayesiano** | 100-200 trials | ~2-4h | Prototipagem rápida |
| **Bayesiano Rápido** | 100 trials × 5 épocas | ~1-2h | Validação de conceito |
| **Híbrido** (Grid + Bayes) | 41,400 + 600 | ~50-76h | Exploração + Refinamento |

## 🎯 Recomendações de Uso

### Para Artigo Qualis A1

```bash

# Execução completa com todos os recursos
python framework_investigativo_completo.py --run-both --trials 150 --dataset-bayes all

```text

- Grid search completo para exploração exaustiva
- Bayesiano adicional para refinamento fino
- Consolida tudo automaticamente
- Tempo: ~50-76 horas


### Para Prototipagem Rápida

```bash

# Windows
$env:VQC_BAYESIAN="1"; $env:VQC_QUICK="1"; python framework_investigativo_completo.py

# Linux/Mac
VQC_BAYESIAN=1 VQC_QUICK=1 python framework_investigativo_completo.py

```text

- Otimização Bayesiana rápida
- Resultados em 1-2 horas
- Útil para validar conceitos


### Para Validação Intermediária

```bash
python framework_investigativo_completo.py --bayes --trials 150 --epocas-bayes 20 --dataset-bayes moons

```text

- Foco em um dataset específico
- Refinamento com mais épocas
- Tempo: ~3-5 horas


## 📝 Notas Técnicas

### Consolidação Robusta
- Lê todos os CSVs individuais, mesmo que alguns estejam corrompidos
- Ignora arquivos com erro e continua processamento
- Valida colunas necessárias antes de gerar comparações


### Detecção Automática de Ambiente
- Suporta execução local, Google Colab e Google Drive
- Variável de ambiente `RESULTS_BASE_DIR` para customizar diretório
- Argumento CLI `--resultados <dir>` para override


### Compatibilidade
- Python 3.9+
- PennyLane 0.38.x
- Optuna (opcional, para modo Bayesiano)
- Todas as dependências em `requirements.txt`


## 🐛 Troubleshooting

### Consolidação não executou
**Problema**: Nenhum CSV individual encontrado.
**Solução**: Verifique se a pasta `experimentos_individuais/` existe e contém arquivos `.csv`.


### Comparação de baselines não gerada
**Problema**: Colunas necessárias ausentes no DataFrame.
**Solução**: Certifique-se de que os CSVs individuais contêm: `dataset`, `arquitetura`, `tipo_ruido`, `nivel_ruido`, `acuracia_teste`.


### Optuna não disponível
**Problema**: Modo Bayesiano falha.
**Solução**:

```bash
pip install optuna

```

## 📚 Referências

- **Framework original**: v7.1 (Beneficial Quantum Noise in VQCs)
- **Automação integrada**: v7.2 (este documento)
- **Artigo**: "From Obstacle to Opportunity: Harnessing Beneficial Quantum Noise in Variational Classifiers"


## 🆘 Suporte

Para dúvidas, issues ou contribuições:

1. Verifique este documento primeiro
2. Consulte `README.md` na raiz do projeto
3. Abra uma issue no repositório GitHub

