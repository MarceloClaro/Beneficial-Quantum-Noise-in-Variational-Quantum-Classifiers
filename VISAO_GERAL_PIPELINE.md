# Visão Geral: Pipeline Automatizada QAOA

## 📊 Diagrama da Pipeline

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                     PIPELINE AUTOMATIZADA QAOA - QUALIS A1                  │
└─────────────────────────────────────────────────────────────────────────────┘

  [1️⃣  EXPERIMENTO]
         ↓
    ┌─────────────────────────┐
    │ experimento_qaoa_        │
    │ otimizado.py            │
    │                         │
    │ • 6-qubit QAOA          │
    │ • 4 variantes de ruído  │
    │ • Tempo: ~2-3 min       │
    │                         │
    │ OUTPUT:                 │
    │ ├─ resultados_YYYYMMDD_ │
    │ │  HHMMSS.csv           │
    │ └─ resultados_YYYYMMDD_ │
    │    HHMMSS.json          │
    └─────────────────────────┘
         ↓
  [2️⃣  ENRIQUECIMENTO]
         ↓
    ┌──────────────────────────┐
    │ enriquecer_resultados_   │
    │ qaoa.py                  │
    │                          │
    │ • Metadados (Run ID,     │
    │   Timestamp, Config)     │
    │ • Métricas derivadas     │
    │   (AUE, TREX, Class)     │
    │ • Resumo + hashes        │
    │                          │
    │ OUTPUT:                  │
    │ └─ CSV enriquecido com   │
    │    20+ colunas + resumo  │
    └──────────────────────────┘
         ↓
  [3️⃣  CONSOLIDAÇÃO]
         ↓
    ┌─────────────────────────┐
    │ auditoria_qaoa_         │
    │ resultados.py           │
    │                         │
    │ • Mescla múltiplos runs │
    │ • Gera gráficos         │
    │   (Energia, Tempo)      │
    │ • Captura ambiente      │
    │ • Exporta hashes SHA256 │
    │                         │
    │ OUTPUT:                 │
    │ ├─ auditoria_qaoa_      │
    │ │  master.csv           │
    │ ├─ *.png *.html         │
    │ ├─ ambiente_execucao    │
    │ │  .json                │
    │ └─ manifest_codigo.json │
    └─────────────────────────┘
         ↓
  [4️⃣  VALIDAÇÃO]
         ↓
    ┌──────────────────────────┐
    │ validar_auditoria_qaoa   │
    │ .py                      │
    │                          │
    │ ✅ QUALIS A1 Compliance  │
    │    • Colunas presentes   │
    │    • Resumo presente     │
    │    • Energia ~100%       │
    │    • Classificação OK    │
    │                          │
    │ OUTPUT:                  │
    │ └─ ✅ PASS ou ❌ FAIL    │
    └──────────────────────────┘
         ↓
  [5️⃣  RASTREABILIDADE]
         ↓
    ┌──────────────────────────┐
    │ calculador_hashes_qaoa   │
    │ .py                      │
    │                          │
    │ SHA-256 dos 5 scripts:   │
    │ • experimento_qaoa_      │
    │   otimizado.py           │
    │ • enriquecer_resultados_ │
    │   qaoa.py                │
    │ • auditoria_qaoa_        │
    │   resultados.py          │
    │ • validar_auditoria_     │
    │   qaoa.py                │
    │ • framework_qaoa_        │
    │   100qubits.py           │
    │                          │
    │ OUTPUT:                  │
    │ └─ manifest_codigo.json  │
    │    (com 5 hashes)        │
    └──────────────────────────┘
         ↓
    ✅ PIPELINE COMPLETA

```

---

## 📁 Arquivos Gerados

### Estrutura de Saída

```
resultados_qaoa_otimizado/
└── resultados_YYYYMMDD_HHMMSS.csv      [CSV enriquecido: metadados + métricas + resumo]
    resultados_YYYYMMDD_HHMMSS.json     [JSON com manifest_codigo]

auditoria_qaoa/
├── auditoria_qaoa_master.csv           [Tabela mestre consolidada]
├── auditoria_qaoa_energia.png          [Gráfico: Energia (plotly)]
├── auditoria_qaoa_energia.html         [Gráfico interativo]
├── auditoria_qaoa_tempo.png            [Gráfico: Tempo (plotly)]
├── auditoria_qaoa_tempo.html           [Gráfico interativo]
├── ambiente_execucao.json              [Python, SO, pacotes]
├── manifest_codigo.json                [SHA-256 dos 5 scripts]
└── README_AUDITORIA_QAOA.md            [Este guia de auditoria]
```

---

## ⏱️ Tempo de Execução Esperado

| Etapa | Tempo | Observações |
| --- | --- | --- |
| 1️⃣ Experimento | 2-3 min | QAOA 6 qubits × 4 variantes |
| 2️⃣ Enriquecimento | 30 seg | Cálculos rápidos em Pandas |
| 3️⃣ Consolidação | 1-2 min | Gráficos, JSON, ambiente |
| 4️⃣ Validação | 10 seg | Checks simples |
| 5️⃣ Rastreabilidade | 5 seg | SHA-256 dos 5 scripts |
| **Total** | **4-7 min** | **Tempo total end-to-end** |

---

## 🎯 Checklist: Como Usar

### ✅ Passo 1: Configurar Ambiente

```powershell
python -m venv .venv
& ".\.venv\Scripts\Activate.ps1"
pip install -r requirements.txt
```

### ✅ Passo 2: Rodar Pipeline Automática

**Opção A: Manualmente (ver cada etapa)**
```powershell
& ".\.venv\Scripts\python.exe" "experimento_qaoa_otimizado.py"
& ".\.venv\Scripts\python.exe" "enriquecer_resultados_qaoa.py"
& ".\.venv\Scripts\python.exe" "auditoria_qaoa_resultados.py"
& ".\.venv\Scripts\python.exe" "validar_auditoria_qaoa.py"
& ".\.venv\Scripts\python.exe" "calculador_hashes_qaoa.py"
```

**Opção B: Tudo junto (mais rápido)**
```powershell
& ".\.venv\Scripts\python.exe" "test_pipeline_automatizada.py"
```

### ✅ Passo 3: Verificar Resultados

**Tabela mestre consolidada**:
```powershell
cat "auditoria_qaoa\auditoria_qaoa_master.csv" | head -20
```

**Manifest de hashes**:
```powershell
cat "auditoria_qaoa\manifest_codigo.json" | python -m json.tool
```

**Validação QUALIS A1**:
```powershell
& ".\.venv\Scripts\python.exe" "validar_auditoria_qaoa.py"
# Esperado: ✅ Validação QUALIS A1: PASS
```

---

## 🔍 O que Cada Script Faz

| Script | Função | Chave de Sucesso |
| --- | --- | --- |
| `experimento_qaoa_otimizado.py` | Executa QAOA 6 qubits × 4 ruídos | CSV + JSON gerados |
| `enriquecer_resultados_qaoa.py` | Adiciona 20+ colunas + resumo | CSV com Hashes SHA-256 |
| `auditoria_qaoa_resultados.py` | Consolida + gera gráficos | Master CSV + PNGs + manifest |
| `validar_auditoria_qaoa.py` | Valida QUALIS A1 (4 checks) | ✅ PASS |
| `calculador_hashes_qaoa.py` | Calcula SHA-256 dos 5 scripts | manifest_codigo.json |
| `test_pipeline_automatizada.py` | Testa toda a pipeline | 5/5 etapas ✅ |

---

## 🛡️ Rastreabilidade de Código

Cada run captura **hashes SHA-256** dos scripts:

```json
{
  "timestamp": "2025-12-28T14:55:00",
  "scripts": {
    "experimento_qaoa_otimizado.py": "146892e1c9f4d54ffcb5d4caf07114f992e691746b7f236df85680cfb7b2fa30",
    "enriquecer_resultados_qaoa.py": "9069b979c7ca335010a838133445558fe62fe7f5d631e883196a069346d8285f",
    "auditoria_qaoa_resultados.py": "2c36bbb229275f5be79b3cf970232fada91c514303af98583f55bde4c5d96148",
    "validar_auditoria_qaoa.py": "aae51347bd1f3d382790fea9be57017b3d4ed0bd62889590d3d6088a50227294",
    "framework_qaoa_100qubits.py": "c9a8ba174ad4ca0a68b2258f52f6985390358fae242ef66fe276e51fc909880d"
  }
}
```

### Verificar Integridade

Se qualquer script for alterado:
```powershell
& ".\.venv\Scripts\python.exe" "calculador_hashes_qaoa.py"
# Compare os hashes com manifest_codigo.json anterior
```

---

## 📊 Métricas Geradas

### Enriquecimento CSV

Cada experimento tem:
- **Energias**: Energia Final, Energia Máxima
- **Normalizado**: Energia Normalizada (%)
- **Comparação**: Melhora vs Sem Ruído (%)
- **Classificação**: "Ruído benéfico" / "prejudicial" / "negligenciável"
- **Heurísticas**: AUE, TREX
- **Metadados**: Run ID, Timestamp, Qubits, P-layers, Shots, Max Iter, Seed

### Exemplo de Linha

```
Experimento,Energia Final,Tempo (s),Run ID,Qubits,Energia Normalizada (%),Melhora vs Sem Ruído (%),Classificação,AUE,TREX
Sem Ruído,-6.50,0.85,run_20251228_145335,6,100.00,0.00,negligenciável,0.00,7.65
Depolarizing,-6.52,0.88,run_20251228_145335,6,100.31,0.31,negligenciável,0.00,7.41
...
RESUMO,,(tempo_total),run_20251228_145335,6,99.50,--,--,--,Hashes: 5/5
```

---

## ✅ Validação QUALIS A1

A pipeline valida **4 requisitos**:

1. **CSV enriquecido** presente com todas colunas
2. **Linha de resumo** com tempo_total e num_experimentos
3. **Energia Normalizada** ≈ 100% (baseline sem ruído)
4. **Classificação** consistente (thresholds ±1%)

**Resultado**: ✅ **PASS** (todos os 4 checks passam)

---

## 🚀 Próximos Passos

### Para Maior Escala

Se quiser aumentar o número de qubits:

1. Editar `experimento_qaoa_otimizado.py`:
```python
CONFIG = ConfigQAOA(
    n_qubits=8,          # ↑ de 6 para 8
    p_layers=1,
    shots=1024,          # ↑ de 512
    max_iter=100,        # ↑ de 50
    seed=42
)
```

2. Re-executar pipeline (será consolidado automaticamente)

### Integração com Artigo QUALIS A1

- Cite `PIPELINE_AUTOMATIZADA_QAOA.md` na seção Methods
- Inclua `auditoria_qaoa_master.csv` como Material Suplementar (S1)
- Referencie `manifest_codigo.json` para reproducibilidade
- Use gráficos `auditoria_qaoa_energia.html` para interatividade no paper

---

## 📞 Troubleshooting

| Problema | Solução |
| --- | --- |
| CSV vazio | Re-executar `experimento_qaoa_otimizado.py` |
| Hashes divergem | Código foi editado; novo hash será capturado |
| Validação FAIL | Verificar CSV tem todas as colunas enriquecidas |
| Gráficos não aparecem | Normal (Agg backend); verifique .html (interativo) |

---

**Status da Pipeline**: ✅ **Pronta para Produção QUALIS A1**
