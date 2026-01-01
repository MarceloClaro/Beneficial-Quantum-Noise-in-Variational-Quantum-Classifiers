# 📖 Playbook: Como Reproduzir a Auditoria QAOA

**Objetivo**: Reproduzir exatamente a auditoria QAOA com rastreabilidade completa  
**Público**: Pesquisadores, revisores QUALIS A1, auditores de código  
**Tempo estimado**: 5-10 minutos (incluindo leitura)

---

## 🎯 Pré-requisitos

### 1. Ambiente Python 3.13+ com Qiskit

```powershell
# Verificar versão
python --version
# Esperado: Python 3.13.x
```

### 2. Dependências instaladas

```powershell
# Verificar instalação (activate venv primeiro)
& ".\.venv\Scripts\Activate.ps1"
pip list | grep -E "qiskit|pandas|numpy|matplotlib|plotly"
```

**Pacotes requeridos**:
- qiskit >= 1.0
- qiskit-aer >= 0.13
- pandas >= 2.0
- numpy >= 1.24
- matplotlib >= 3.7
- plotly >= 5.0
- scipy >= 1.10

### 3. Diretórios esperados

```
projeto/
├── experimento_qaoa_otimizado.py
├── enriquecer_resultados_qaoa.py
├── auditoria_qaoa_resultados.py
├── validar_auditoria_qaoa.py
├── calculador_hashes_qaoa.py
├── framework_qaoa_100qubits.py
├── test_pipeline_verificacao.py
├── resultados_qaoa_otimizado/          (será criado)
└── auditoria_qaoa/                     (será criado)
```

---

## 🚀 Reprodução Passo a Passo

### Passo 1: Ativar Ambiente Python

```powershell
cd "seu\caminho\do\projeto"
& ".\.venv\Scripts\Activate.ps1"
```

### Passo 2: Verificar Integridade de Código

**Antes de começar**, verifique os hashes SHA-256:

```powershell
python calculador_hashes_qaoa.py
```

**Saída esperada**:
```
📋 Manifest de Hashes SHA-256 (Rastreabilidade de Código)
======================================================================
experimento_qaoa_otimizado.py         | 146892e1c9f4d54f...
enriquecer_resultados_qaoa.py         | 9069b979c7ca335...
auditoria_qaoa_resultados.py          | 2c36bbb229275f5b...
validar_auditoria_qaoa.py             | aae51347bd1f3d38...
framework_qaoa_100qubits.py           | c9a8ba174ad4ca0a...
✅ Manifest salvo em: .../manifest_codigo.json
```

**Ação**: Anote os hashes. Se forem diferentes do seu projeto, significa que o código foi modificado.

### Passo 3: Executar Experimento QAOA

```powershell
python experimento_qaoa_otimizado.py
```

**Saída esperada**:
```
🔬 QAOA com Ruído Quântico - Análise Benéfico Prejudicial
Configuração: 6 qubits, 1 P-layer, 512 shots, seed=42
...
✅ Experimento sem ruído: Energia Final = -6.50
✅ Experimento com Depolarizing: Energia Final = -6.52
✅ Experimento com Phase Damping: Energia Final = -6.51
✅ Experimento com Amplitude Damping: Energia Final = -6.50

✅ Resultados salvos em: .../resultados_20251228_HHMMSS.csv
✅ JSON salvo em: .../resultados_20251228_HHMMSS.json
```

**Tempo**: ~2-3 minutos  
**Artefatos gerados**:
- `resultados_qaoa_otimizado/resultados_YYYYMMDD_HHMMSS.csv`
- `resultados_qaoa_otimizado/resultados_YYYYMMDD_HHMMSS.json`

### Passo 4: Enriquecer CSV com Metadados e Métricas

```powershell
python enriquecer_resultados_qaoa.py
```

**Saída esperada**:
```
📊 Enriquecendo resultados com metadados e métricas...
✅ Adicionadas colunas: Run ID, Timestamp, Qubits, P-layers, Shots, Max Iter, Seed
✅ Adicionadas métricas: Energia Normalizada, Melhora vs Sem Ruído, Classificação, AUE, TREX
✅ Resumo adicionado com tempo_total e num_experimentos
✅ Hashes SHA-256 inclusos

CSV atualizado: .../resultados_20251228_HHMMSS.csv
```

**Tempo**: ~30 segundos  
**Modificações**:
- 20+ colunas adicionadas ao CSV
- Linha de resumo inserida
- Hashes SHA-256 capturados

### Passo 5: Consolidar e Gerar Gráficos

```powershell
python auditoria_qaoa_resultados.py
```

**Saída esperada**:
```
📋 Consolidando resultados de múltiplos runs...
✅ Consolidado salvo em: .../auditoria_qaoa_master.csv
✅ Gráficos salvos em: .../auditoria_qaoa/
✅ Ambiente de execução documentado em: .../ambiente_execucao.json
✅ Hashes SHA-256 salvos em: .../manifest_codigo.json
```

**Tempo**: ~1-2 minutos  
**Artefatos gerados**:
- `auditoria_qaoa/auditoria_qaoa_master.csv` (tabela consolidada)
- `auditoria_qaoa/auditoria_qaoa_energia.png` (gráfico)
- `auditoria_qaoa/auditoria_qaoa_energia.html` (interativo)
- `auditoria_qaoa/auditoria_qaoa_tempo.png` (gráfico)
- `auditoria_qaoa/auditoria_qaoa_tempo.html` (interativo)
- `auditoria_qaoa/ambiente_execucao.json` (metadata)
- `auditoria_qaoa/manifest_codigo.json` (hashes)

### Passo 6: Validar Conformidade QUALIS A1

```powershell
python validar_auditoria_qaoa.py
```

**Saída esperada**:
```
[PASS] resultados_20251228_HHMMSS.csv possui linha de resumo
[PASS] resultados_20251228_HHMMSS.csv Energia Normalizada máxima ≈ 100.00
[PASS] Master consolidado contém o run atual e as colunas essenciais

✅ Validação QUALIS A1: PASS
```

**Tempo**: ~10 segundos  
**Validações**:
1. ✅ CSV enriquecido com colunas
2. ✅ Linha de resumo presente
3. ✅ Energia normalizada ~100%
4. ✅ Classificação consistente

### Passo 7: Verificar Toda a Pipeline

```powershell
python test_pipeline_verificacao.py
```

**Saída esperada**:
```
✅ VERIFICAÇÃO PIPELINE QAOA
...
Total............................................. 5/5

✅ PIPELINE FUNCIONANDO CORRETAMENTE!
Todos os artefatos foram gerados com sucesso.

Próximos passos:
  1. Revisar auditoria_qaoa/auditoria_qaoa_master.csv
  2. Abrir auditoria_qaoa/auditoria_qaoa_energia.html (interativo)
  3. Verificar manifest_codigo.json para hashes SHA-256
  4. Consultar PIPELINE_AUTOMATIZADA_QAOA.md para detalhes
```

**Tempo**: ~30 segundos  
**Verificações**: 5/5 PASS

---

## 📋 Checklist de Reprodução

Após completar todos os passos, verifique:

- [ ] **Passo 2**: Hashes SHA-256 anotados (compare com commit anterior)
- [ ] **Passo 3**: CSV de resultado gerado e contém 4 experimentos
- [ ] **Passo 4**: CSV agora tem ~20 colunas + linha de resumo
- [ ] **Passo 5**: Diretório `auditoria_qaoa/` contém 8 arquivos
- [ ] **Passo 6**: Validação retorna ✅ PASS
- [ ] **Passo 7**: Teste retorna 5/5 verificações

---

## 🔍 Verificações Críticas

### Verificar CSV Resultado

```powershell
# Ver primeiras linhas
python -c "import pandas as pd; df = pd.read_csv('resultados_qaoa_otimizado/resultados_*.csv'); print(df.head()); print(f'Colunas: {len(df.columns)}')"
```

**Esperado**:
- 4 linhas de experimentos (Sem Ruído, Depolarizing, Phase Damping, Amplitude Damping)
- 1 linha de resumo (RESUMO)
- ~20+ colunas incluindo Energia, Tempo, Run ID, AUE, TREX, Classificação

### Verificar Master CSV

```powershell
cat auditoria_qaoa/auditoria_qaoa_master.csv | head -5
```

**Esperado**: Consolidação de múltiplos runs (se houver)

### Verificar Manifest de Hashes

```powershell
cat auditoria_qaoa/manifest_codigo.json
```

**Esperado**:
```json
{
  "experimento_qaoa_otimizado.py": "146892e1c9f4d54f...",
  "enriquecer_resultados_qaoa.py": "9069b979c7ca335...",
  ...
}
```

### Verificar Ambiente Documentado

```powershell
python -c "import json; d=json.load(open('auditoria_qaoa/ambiente_execucao.json')); print(d['python_version']); print(d['sistema_operacional'])"
```

**Esperado**:
- Python 3.13.x
- Windows-10 ou similar

---

## 🛑 Troubleshooting

| Problema | Solução |
| --- | --- |
| `No columns to parse from file` | CSV corrompido; deletar e re-executar experimento |
| `Qiskit not found` | `pip install qiskit qiskit-aer` |
| Hashes divergem | Código foi modificado; compare com git diff |
| Validação FAIL | Verificar se CSV tem todas as colunas; re-enriquecer |
| Gráficos não aparecem | Normal (Agg backend); verifique arquivos .html |

---

## 📊 Exemplo de Saída Esperada

### Resumo final do CSV enriquecido

```
Experimento,Energia Final,Tempo (s),Run ID,...,Classificação,AUE,TREX
Sem Ruído,-6.50,0.85,run_20251228_145335,...,negligenciável,0.00,7.65
Depolarizing,-6.52,0.88,run_20251228_145335,...,negligenciável,0.00,7.41
Phase Damping,-6.51,0.87,run_20251228_145335,...,negligenciável,0.00,7.48
Amplitude Damping,-6.50,0.86,run_20251228_145335,...,negligenciável,0.00,7.56
RESUMO,,,(tempo_total),run_20251228_145335,...,(resumo),...,Hashes: 5/5
```

### Gráfico Energia.html

Abrir `auditoria_qaoa/auditoria_qaoa_energia.html` no navegador para ver gráfico interativo com:
- Eixo X: Experimento (Sem Ruído, Depolarizing, Phase Damping, Amplitude Damping)
- Eixo Y: Energia Final
- Hover: Tempo, Classificação, AUE, TREX

---

## ✅ Sucesso!

Se todos os 7 passos executarem sem erros e as verificações passarem, você tem uma **auditoria completa e rastreável** pronta para:

1. **Inclusão em artigo QUALIS A1**
2. **Reproducibilidade garantida** (SHA-256)
3. **Auditoria independente** (hashes + metadados)
4. **Conformidade científica** (validação automática)

---

## 📞 Referências

- [PIPELINE_AUTOMATIZADA_QAOA.md](PIPELINE_AUTOMATIZADA_QAOA.md) — Detalhes técnicos
- [VISAO_GERAL_PIPELINE.md](VISAO_GERAL_PIPELINE.md) — Diagrama e arquitetura
- [auditoria_qaoa/README_AUDITORIA_QAOA.md](auditoria_qaoa/README_AUDITORIA_QAOA.md) — Metodologia
- [RESUMO_EXECUTIVO_PIPELINE_QAOA.md](RESUMO_EXECUTIVO_PIPELINE_QAOA.md) — Status geral

**Status**: ✅ **PRONTA PARA REPRODUCIBILIDADE E AUDITORIA QUALIS A1**
