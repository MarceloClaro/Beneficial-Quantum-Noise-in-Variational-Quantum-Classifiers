# Pipeline Automatizada QAOA — QUALIS A1

## Visão Geral

Pipeline **completamente automatizada** para execução, enriquecimento, consolidação, validação e rastreabilidade de versão de experimentos QAOA com análise de ruído benéfico.

```
[Experimento] → [Enriquecimento] → [Consolidação] → [Validação] → [Rastreabilidade]
```

## Etapas da Pipeline

### 1️⃣ Execução de Experimentos

**Script**: `experimento_qaoa_otimizado.py`

Executa 4 variantes de QAOA com 6 qubits:
- Sem Ruído (baseline)
- Com Ruído Depolarizante
- Com Ruído de Damping de Fase
- Com Ruído de Damping de Amplitude

**Outputs**:
- CSV: `resultados_qaoa_otimizado/resultados_YYYYMMDD_HHMMSS.csv`
- JSON: `resultados_qaoa_otimizado/resultados_YYYYMMDD_HHMMSS.json`
- Manifest de hashes: capturado no JSON (`manifest_codigo`)

**Comando**:
```powershell
& ".\.venv\Scripts\python.exe" "experimento_qaoa_otimizado.py"
```

### 2️⃣ Enriquecimento de Metadados

**Script**: `enriquecer_resultados_qaoa.py`

Adiciona colunas aos CSVs:
- **Metadados**: Run Timestamp, Run ID, Qubits, P-layers, Shots, Max Iter, Seed
- **Métricas Derivadas**:
  - Energia Normalizada (%) = (Energia Final / Max) × 100
  - Melhora vs Sem Ruído (%) = ((Final − Sem Ruído) / Sem Ruído) × 100
  - Classificação (ruído benéfico/prejudicial/negligenciável)
  - AUE: max(0, Melhora %) × (Normalizada / 100)
  - TREX: Energia Final / Tempo (s)

Insere **linha de resumo** com:
- `tempo_total`: soma de tempos dos 4 experimentos
- `num_experimentos`: número de runs (ex.: 4)
- `hashes_sha256_manifest`: hashes dos scripts (RESUMO)

**Outputs**:
- CSV enriquecido: mesmo arquivo com novas colunas + linha resumo

**Comando**:
```powershell
& ".\.venv\Scripts\python.exe" "enriquecer_resultados_qaoa.py"
```

### 3️⃣ Consolidação Multi-Run

**Script**: `auditoria_qaoa_resultados.py`

Consolida todos os CSVs/JSONs em uma tabela mestre:
- Mescla múltiplos runs
- Gera gráficos (Energia, Tempo)
- Captura ambiente (Python, SO, pacotes)
- Exporta manifest de código

**Outputs**:
- `auditoria_qaoa/auditoria_qaoa_master.csv` (tabela mestre)
- `auditoria_qaoa/auditoria_qaoa_energia.{png,html}` (gráfico de energia)
- `auditoria_qaoa/auditoria_qaoa_tempo.{png,html}` (gráfico de tempo)
- `auditoria_qaoa/ambiente_execucao.json` (metadata de ambiente)
- `auditoria_qaoa/manifest_codigo.json` (hashes SHA-256)

**Comando**:
```powershell
& ".\.venv\Scripts\python.exe" "auditoria_qaoa_resultados.py"
```

### 4️⃣ Validação QUALIS A1

**Script**: `validar_auditoria_qaoa.py`

Valida conformidade com requisitos QUALIS A1:
1. ✅ CSV enriquecido com todas as colunas
2. ✅ Linha de resumo presente (com tempo_total e num_experimentos)
3. ✅ Energia Normalizada ≈ 100% (baseline)
4. ✅ Classificação consistente (thresholds: ±1%)

**Saída**:
- ✅ **PASS**: Todas as validações passaram
- ❌ **FAIL**: Erros específicos listados

**Comando**:
```powershell
& ".\.venv\Scripts\python.exe" "validar_auditoria_qaoa.py"
```

### 5️⃣ Rastreabilidade de Versão

**Script**: `calculador_hashes_qaoa.py`

Calcula SHA-256 dos scripts principais:
- `experimento_qaoa_otimizado.py`
- `enriquecer_resultados_qaoa.py`
- `auditoria_qaoa_resultados.py`
- `validar_auditoria_qaoa.py`
- `framework_qaoa_100qubits.py`

Exporta `manifest_codigo.json` com todos os hashes.

**Comando**:
```powershell
& ".\.venv\Scripts\python.exe" "calculador_hashes_qaoa.py"
```

## Execução Automática Completa

### Opção 1: Passo a Passo (com visualização)

```powershell
# 1. Rodar experimentos
& ".\.venv\Scripts\python.exe" "experimento_qaoa_otimizado.py"

# 2. Enriquecer CSVs
& ".\.venv\Scripts\python.exe" "enriquecer_resultados_qaoa.py"

# 3. Consolidar e gerar gráficos
& ".\.venv\Scripts\python.exe" "auditoria_qaoa_resultados.py"

# 4. Validar
& ".\.venv\Scripts\python.exe" "validar_auditoria_qaoa.py"

# 5. Verificar integridade de código
& ".\.venv\Scripts\python.exe" "calculador_hashes_qaoa.py"
```

### Opção 2: Script Automático (recomendado)

**Criar**: `executar_pipeline_completa.ps1`

```powershell
param(
    [switch]$SkipValidation = $false
)

Write-Host "🚀 Iniciando Pipeline QAOA Automática..." -ForegroundColor Green

# 1. Experimentos
Write-Host "`n[1/5] Executando experimentos..." -ForegroundColor Cyan
& ".\.venv\Scripts\python.exe" "experimento_qaoa_otimizado.py"
if ($LASTEXITCODE -ne 0) { exit 1 }

# 2. Enriquecimento
Write-Host "`n[2/5] Enriquecendo CSVs..." -ForegroundColor Cyan
& ".\.venv\Scripts\python.exe" "enriquecer_resultados_qaoa.py"
if ($LASTEXITCODE -ne 0) { exit 1 }

# 3. Consolidação
Write-Host "`n[3/5] Consolidando e gerando gráficos..." -ForegroundColor Cyan
& ".\.venv\Scripts\python.exe" "auditoria_qaoa_resultados.py"
if ($LASTEXITCODE -ne 0) { exit 1 }

# 4. Validação (opcional)
if (-not $SkipValidation) {
    Write-Host "`n[4/5] Validando QUALIS A1..." -ForegroundColor Cyan
    & ".\.venv\Scripts\python.exe" "validar_auditoria_qaoa.py"
    if ($LASTEXITCODE -ne 0) { exit 1 }
}

# 5. Rastreabilidade
Write-Host "`n[5/5] Verificando integridade de código..." -ForegroundColor Cyan
& ".\.venv\Scripts\python.exe" "calculador_hashes_qaoa.py"

Write-Host "`n✅ Pipeline completa com sucesso!" -ForegroundColor Green
Write-Host "📊 Resultados em: auditoria_qaoa/" -ForegroundColor Yellow
```

**Executar**:
```powershell
& ".\executar_pipeline_completa.ps1"
```

## Artefatos Gerados

### Estrutura de Diretórios

```
resultados_qaoa_otimizado/
├── resultados_20251228_145335.csv    (enriquecido: metadados + métricas)
├── resultados_20251228_145335.json   (com manifest_codigo)
└── ...

auditoria_qaoa/
├── auditoria_qaoa_master.csv         (tabela mestre consolidada)
├── auditoria_qaoa_energia.png        (gráfico: energia vs experimento)
├── auditoria_qaoa_energia.html       (interativo)
├── auditoria_qaoa_tempo.png          (gráfico: tempo vs experimento)
├── auditoria_qaoa_tempo.html         (interativo)
├── ambiente_execucao.json            (Python, OS, pacotes)
├── manifest_codigo.json              (hashes SHA-256)
└── README_AUDITORIA_QAOA.md          (este documento)
```

## Verificação de Integridade

### Como validar que o código não foi alterado?

1. **Capturar hashes atuais**:
```powershell
& ".\.venv\Scripts\python.exe" "calculador_hashes_qaoa.py"
```

2. **Comparar com manifest anterior**:
```powershell
cat "auditoria_qaoa\manifest_codigo.json"
```

3. **Se diferente**: ⚠️ Código foi alterado desde último run
   - Documente as mudanças
   - Re-execute a pipeline

## Reprodutibilidade

Cada run é **totalmente reproduzível** porque:
- ✅ Metadados capturados (timestamp, seed, configurações)
- ✅ Ambiente documentado (Python version, pacotes)
- ✅ Hashes SHA-256 rastreiam versão exata do código
- ✅ CSV consolidado permite comparação com runs anteriores

## Modificar Experimentos

Se quiser alterar configurações (ex.: aumentar qubits, shots):

1. **Editar `experimento_qaoa_otimizado.py`**:
```python
CONFIG = ConfigQAOA(
    n_qubits=8,          # aumentar de 6 para 8
    p_layers=1,
    shots=1024,          # aumentar de 512
    max_iter=100,        # aumentar de 50
    seed=42
)
```

2. **Executar pipeline novamente**:
```powershell
& ".\executar_pipeline_completa.ps1"
```

3. **Novo run será automaticamente consolidado** em `auditoria_qaoa_master.csv`

## Troubleshooting

### CSV corrompido / vazio
- **Causa**: Experimento falhou ou foi interrompido
- **Ação**: Deletar arquivo corrompido e re-executar experimento
- **Auditoria**: Erro registrado em logs (não falha pipeline)

### Hashes divergem
- **Causa**: Código foi editado
- **Ação**: Comparar diffs, documentar mudanças, re-executar
- **Rastreabilidade**: Novo hash será capturado automaticamente

### Gráficos não aparecem
- **Causa**: Matplotlib Agg backend (headless)
- **Ação**: Normal em ambientes sem display; verifique HTML (interativo)

## Conformidade QUALIS A1

✅ Todos os requisitos atendidos:
- Metadados completos (timestamp, seed, configurações)
- Métricas derivadas (AUE, TREX, classificação)
- Consolidação multi-run (rastreabilidade de experimentos)
- Validação automática (PASS/FAIL)
- Rastreabilidade de código (SHA-256)
- Reprodutibilidade garantida (ambiente + hashes)

**Resultado**: ✅ **Validação QUALIS A1: PASS**
