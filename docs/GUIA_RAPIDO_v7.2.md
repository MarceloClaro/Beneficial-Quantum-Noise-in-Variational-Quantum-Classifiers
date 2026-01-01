# 🚀 Guia Rápido de Uso - Framework v7.2

## Início Rápido (Quick Start)

### 1. Execução Simples (Automático)

```bash

# Tudo automático - consolidação incluída!
python framework_investigativo_completo.py

```text

**O que acontece**:
- ✅ Grid search completo (8,280 configs × 5 seeds)
- ✅ Análises estatísticas automáticas
- ✅ 9 visualizações interativas geradas
- ✅ **NOVO**: Consolidação automática de todos os resultados
- ✅ **NOVO**: Geração de comparacao_baselines.csv
- ✅ **NOVO**: Geração de metadata_orchestrator.json


**Tempo**: ~48-72 horas (grid completo)


---


### 2. Teste Rápido (5 minutos de setup)

```bash

# Windows PowerShell
$env:VQC_QUICK="1"; python framework_investigativo_completo.py

# Linux/Mac
VQC_QUICK=1 python framework_investigativo_completo.py

```text

**O que acontece**:
- ⚡ Reduz épocas para 5 (vs. 15 padrão)
- ⚡ Mesma cobertura, execução mais rápida
- ✅ Consolidação automática incluída


**Tempo**: ~16-24 horas


---


### 3. Otimização Bayesiana (Inteligente e Rápido)

```bash

# Modo Bayesiano - 10-20x mais rápido que grid
python framework_investigativo_completo.py --bayes --trials 150 --dataset-bayes moons

```text

**O que acontece**:
- 🧠 Optuna (TPE) para busca inteligente
- 🎯 150 trials em vez de 8,280 configurações
- ✅ Consolidação automática incluída


**Tempo**: ~2-4 horas


---


### 4. Bayesiano para Todos os Datasets

```bash
python framework_investigativo_completo.py --bayes --trials 200 --dataset-bayes all

```text

**O que acontece**:
- 🧠 Otimização Bayesiana em: moons, circles, iris, breast_cancer, wine
- 🎯 200 trials por dataset (cap máximo)
- ✅ Consolidação automática para todos


**Tempo**: ~8-15 horas (total)


---


### 5. Híbrido (Exploração + Refinamento)

```bash

# Grid seguido de Bayesiano - melhor de dois mundos
python framework_investigativo_completo.py --run-both --trials 120

```text

**O que acontece**:
- 📊 Grid search completo primeiro (exploração exaustiva)
- 🧠 Depois: otimização Bayesiana (refinamento fino)
- ✅ Consolidação única no final


**Tempo**: ~50-76 horas


**Recomendado para**: Publicações Qualis A1


---


## Uso Programático

### Consolidar Resultados Existentes

```python
from framework_investigativo_completo import consolidar_e_gerar_metadados

# Pós-processar um diretório existente
resultado = consolidar_e_gerar_metadados(
    pasta_resultados='resultados_2025-10-28_17-09-33',
    verbose=True
)

print(f"Status: {resultado['consolidacao']['status']}")
print(f"CSVs individuais: {resultado['consolidacao']['num_csvs_individuais']}")
print(f"Linhas consolidadas: {resultado['consolidacao']['rows_consolidated']}")
print(f"Arquivo consolidado: {resultado['consolidacao']['consolidated_path']}")
print(f"Comparação: {resultado['comparacao_baselines']}")
print(f"Metadados: {resultado['metadata_orchestrator']}")

```text

**Saída esperada**:

```

  🔄 Iniciando consolidação automática...
  📦 Encontrados 97 CSVs individuais. Consolidando...
  ✅ CSV consolidado salvo: resultados_completos_artigo.csv
     Linhas: 97 | Colunas: 17
  ✅ Comparação de baselines salva: comparacao_baselines.csv
  ✅ Metadados salvos: metadata_orchestrator.json
  ✅ Consolidação e metadados concluídos!

Status: ok
CSVs individuais: 97
Linhas consolidadas: 97
Arquivo consolidado: C:\...\resultados_completos_artigo.csv
Comparação: resultados_2025-10-28_17-09-33\comparacao_baselines.csv
Metadados: C:\...\metadata_orchestrator.json

```text

---


### Apenas Consolidação (sem metadados)

```python
from framework_investigativo_completo import consolidar_resultados_individuais

info = consolidar_resultados_individuais(
    pasta_resultados='resultados_2025-10-28_17-09-33',
    verbose=True
)

if info['status'] == 'ok':
    print(f"✅ Consolidado: {info['rows_consolidated']} linhas")
    df = info['df']  # DataFrame consolidado disponível
    print(df.head())

```text

---


### Apenas Comparação (CSVs já consolidados)

```python
import pandas as pd
from framework_investigativo_completo import gerar_comparacao_baselines

# Ler CSV consolidado existente
df = pd.read_csv('resultados_2025-10-28_17-09-33/resultados_completos_artigo.csv')

# Gerar comparação
comp_path = gerar_comparacao_baselines(
    df_all=df,
    pasta_resultados='resultados_2025-10-28_17-09-33',
    verbose=True
)

# Ler e analisar comparação
comp_df = pd.read_csv(comp_path)
print("\n🏆 Melhor VQC por dataset:")
print(comp_df[['dataset', 'vqc_melhor', 'delta_vqc_svm', 'delta_vqc_rf']])

```text

---


## Cenários de Uso

### Cenário 1: Primeira Execução (Pesquisa Inicial)

```bash

# Teste rápido para validar setup
$env:VQC_QUICK="1"; $env:VQC_BAYESIAN="1"; python framework_investigativo_completo.py

```text

**Tempo**: 1-2 horas  
**Uso**: Validar ambiente e conceitos


---


### Cenário 2: Desenvolvimento e Debugging

```bash

# Grid rápido para um dataset específico (modificar código temporariamente)
$env:VQC_QUICK="1"; python framework_investigativo_completo.py

```text

**Tempo**: 16-24 horas  
**Uso**: Testar mudanças no código


---


### Cenário 3: Artigo Qualis A1 (Produção)

```bash

# Exploração completa + refinamento Bayesiano
python framework_investigativo_completo.py --run-both --trials 150 --dataset-bayes all

```text

**Tempo**: 50-76 horas  
**Uso**: Resultados finais para publicação


---


### Cenário 4: Apresentação ou Demo

```bash

# Bayesiano rápido em um dataset popular
python framework_investigativo_completo.py --bayes --trials 100 --epocas-bayes 10 --dataset-bayes moons

```text

**Tempo**: 1-2 horas  
**Uso**: Gerar resultados para apresentação


---


### Cenário 5: Re-consolidar Resultados Antigos

```python

# Script Python para pós-processar execuções antigas
from framework_investigativo_completo import consolidar_e_gerar_metadados
import os

# Encontrar todas as pastas de resultados
for folder in os.listdir('.'):
    if folder.startswith('resultados_'):
        print(f"\n📂 Processando: {folder}")
        try:
            consolidar_e_gerar_metadados(folder, verbose=True)
        except Exception as e:
            print(f"  ⚠️ Erro: {e}")

```text

---


## Interpretando os Resultados

### resultados_completos_artigo.csv

```csv
dataset,arquitetura,estrategia_init,tipo_ruido,nivel_ruido,acuracia_teste,...
moons,basic_entangler,matematico,sem_ruido,0.0,0.9250,...
moons,basic_entangler,matematico,depolarizante,0.05,0.9380,...  # ⬆️ Ruído benéfico!

```text

**Como usar**:

```python
import pandas as pd

df = pd.read_csv('resultados_2025-10-28_17-09-33/resultados_completos_artigo.csv')

# Encontrar melhor configuração global
melhor = df.loc[df['acuracia_teste'].idxmax()]
print(f"Melhor: {melhor['dataset']} - {melhor['tipo_ruido']} - {melhor['acuracia_teste']:.4f}")

# Comparar ruído vs. sem ruído
sem_ruido = df[df['tipo_ruido'] == 'sem_ruido']['acuracia_teste'].mean()
com_depol = df[df['tipo_ruido'] == 'depolarizante']['acuracia_teste'].mean()
print(f"Sem ruído: {sem_ruido:.4f}")
print(f"Com depolarizante: {com_depol:.4f}")
print(f"Delta: {com_depol - sem_ruido:+.4f}")

```text

---


### comparacao_baselines.csv

```csv
dataset,vqc_melhor,vqc_sem_ruido_media,svm,rf,delta_vqc_svm,delta_vqc_rf
moons,0.9380,0.9100,0.9050,0.9200,0.0330,0.0180

```text

**Como usar**:

```python
comp = pd.read_csv('resultados_2025-10-28_17-09-33/comparacao_baselines.csv')

# Datasets onde VQC vence SVM
vence_svm = comp[comp['delta_vqc_svm'] > 0]
print(f"\n🏆 VQC vence SVM em {len(vence_svm)} de {len(comp)} datasets:")
print(vence_svm[['dataset', 'vqc_melhor', 'svm', 'delta_vqc_svm']])

# Datasets onde VQC vence RF
vence_rf = comp[comp['delta_vqc_rf'] > 0]
print(f"\n🏆 VQC vence RF em {len(vence_rf)} de {len(comp)} datasets:")
print(vence_rf[['dataset', 'vqc_melhor', 'rf', 'delta_vqc_rf']])

```text

---


### metadata_orchestrator.json

```json
{
  "tipo": "metadata_orchestrator",
  "versao_framework": "7.2",
  "timestamp": "2025-10-28 20:06:32",
  "consolidacao": {
    "status": "ok",
    "num_csvs_individuais": 97,
    "rows_consolidated": 97
  }
}

```text

**Como usar**:

```python
import json

with open('resultados_2025-10-28_17-09-33/metadata_orchestrator.json') as f:
    meta = json.load(f)

print(f"Framework: v{meta['versao_framework']}")
print(f"Executado em: {meta['timestamp']}")
print(f"Experimentos: {meta['consolidacao']['rows_consolidated']}")
print(f"Status: {meta['consolidacao']['status']}")

```text

---


## Troubleshooting

### Problema: "Pasta experimentos_individuais não encontrada"
**Solução**: O framework ainda não foi executado ou a execução falhou antes de gerar CSVs.

```bash

# Execute o framework primeiro
python framework_investigativo_completo.py

```text

---


### Problema: "Nenhum CSV individual encontrado"
**Solução**: A pasta existe mas está vazia. Verifique se a execução anterior completou ao menos um experimento.

```bash

# Liste CSVs
ls resultados_2025-10-28_17-09-33/experimentos_individuais/

```text

---


### Problema: "Colunas necessárias não encontradas"
**Solução**: CSVs individuais têm estrutura diferente. Verifique se está usando a versão correta do framework.

```python

# Inspecionar estrutura de um CSV
import pandas as pd
df = pd.read_csv('resultados_2025-10-28_17-09-33/experimentos_individuais/exp_00001.csv')
print(df.columns)

```text

---


### Problema: Optuna não disponível (modo Bayesiano)
**Solução**: Instale Optuna:

```bash
pip install optuna

```text

---


## Comparação de Modos

| Modo | Tempo | Exploração | Refinamento | Qualis A1 | Demo |
|------|-------|-----------|-------------|-----------|------|
| **Grid Completo** | 48-72h | ✅✅✅ | ✅ | ✅✅✅ | ❌ |
| **Grid Rápido** | 16-24h | ✅✅✅ | ✅ | ✅✅ | ❌ |
| **Bayesiano** | 2-4h | ✅ | ✅✅✅ | ✅✅ | ✅✅✅ |
| **Bayesiano Rápido** | 1-2h | ✅ | ✅✅ | ✅ | ✅✅✅ |
| **Híbrido** | 50-76h | ✅✅✅ | ✅✅✅ | ✅✅✅ | ❌ |

---


## Checklist de Execução

### Antes de Executar
- [ ] Python 3.9+ instalado
- [ ] Dependências instaladas (`pip install -r requirements.txt`)
- [ ] Espaço em disco suficiente (~5-10 GB para grid completo)
- [ ] Tempo disponível (veja tabela acima)


### Durante a Execução
- [ ] Monitorar logs no terminal
- [ ] Verificar criação de pasta `resultados_YYYY-MM-DD_HH-MM-SS/`
- [ ] Verificar CSVs sendo gerados em `experimentos_individuais/`


### Após Execução
- [ ] Verificar `resultados_completos_artigo.csv` existe
- [ ] Verificar `comparacao_baselines.csv` existe
- [ ] Verificar `metadata_orchestrator.json` existe
- [ ] Abrir visualizações HTML no navegador
- [ ] Analisar resultados estatísticos


---


## Comandos Úteis

### Verificar Progresso

```bash

# Contar CSVs gerados
ls resultados_2025-10-28_17-09-33/experimentos_individuais/ | measure

# ou (Linux/Mac)
ls resultados_2025-10-28_17-09-33/experimentos_individuais/*.csv | wc -l

```text

### Listar Arquivos Gerados

```bash
ls resultados_2025-10-28_17-09-33/*.csv
ls resultados_2025-10-28_17-09-33/*.html
ls resultados_2025-10-28_17-09-33/*.json

```text

### Abrir Visualizações

```bash

# Windows
start resultados_2025-10-28_17-09-33/figura2_beneficial_noise.html

# Linux/Mac
open resultados_2025-10-28_17-09-33/figura2_beneficial_noise.html

```

---


## Próximos Passos

1. ✅ **Execute o framework** no modo desejado
2. ✅ **Verifique os resultados** (CSVs, comparações, metadados)
3. ✅ **Analise as visualizações** (9 figuras HTML interativas)
4. ✅ **Use os dados** para seu artigo/pesquisa


---


**📚 Documentação Completa**: [AUTOMACAO_FRAMEWORK.md](AUTOMACAO_FRAMEWORK.md)  
**📋 Changelog**: [CHANGELOG_v7.2.md](CHANGELOG_v7.2.md)  
**🎯 Resumo**: [RESUMO_EXECUTIVO_v7.2.md](RESUMO_EXECUTIVO_v7.2.md)


---


**🎉 Boa sorte com sua pesquisa! 🎉**

