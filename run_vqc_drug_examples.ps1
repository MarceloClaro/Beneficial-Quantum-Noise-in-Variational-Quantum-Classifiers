# ============================================================================
# run_vqc_drug_examples.ps1
# Execute multiple VQC-Molecular experiments (Windows PowerShell)
# ============================================================================

$ErrorActionPreference = "Stop"

Write-Host "╔════════════════════════════════════════════════════════════════════╗" -ForegroundColor Cyan
Write-Host "║     VQC-Molecular v8.0 - Drug Screening Pipeline Examples          ║" -ForegroundColor Cyan
Write-Host "║  Quantum-Enhanced QSAR with Automatic Hyper-parameter Tuning       ║" -ForegroundColor Cyan
Write-Host "╚════════════════════════════════════════════════════════════════════╝" -ForegroundColor Cyan

# Create directories
@("results_vqc_drug", "qsar_cache", "logs") | ForEach-Object {
    if (-not (Test-Path $_)) {
        New-Item -ItemType Directory -Path $_ -Force | Out-Null
        Write-Host "✓ Criado: $_" -ForegroundColor Green
    }
}

$TIMESTAMP = Get-Date -Format "yyyy-MM-dd_HH-mm-ss"
$LOG_FILE = "logs/vqc_drug_${TIMESTAMP}.log"

Write-Host ""
Write-Host "$(Get-Date) ========== Verificando Dependências ==========" -ForegroundColor Yellow

try {
    $version = python -c "import pennylane as qml; print(qml.__version__)" 2>$null
    Write-Host "✓ PennyLane $version" -ForegroundColor Green
}
catch {
    Write-Host "✗ PennyLane não encontrado!" -ForegroundColor Red
    exit 1
}

try {
    $version = python -c "import optuna; print(optuna.__version__)" 2>$null
    Write-Host "✓ Optuna $version" -ForegroundColor Green
}
catch {
    Write-Host "✗ Optuna não encontrado!" -ForegroundColor Red
}

try {
    $version = python -c "import deepchem as dc; print(dc.__version__)" 2>$null
    Write-Host "✓ DeepChem $version" -ForegroundColor Green
}
catch {
    Write-Host "✗ DeepChem não encontrado!" -ForegroundColor Red
}

try {
    python -c "import rdkit" 2>$null
    Write-Host "✓ RDKit OK" -ForegroundColor Green
}
catch {
    Write-Host "✗ RDKit não encontrado!" -ForegroundColor Red
}

Write-Host ""

# ============================================================================
# EXEMPLO 1: EGFR (Piloto)
# ============================================================================
Write-Host "$(Get-Date) ========== EXEMPLO 1: EGFR Kinase (Piloto) ==========" -ForegroundColor Yellow
Write-Host "Target: EGFR (6,847 moléculas)" 
Write-Host "Configuração: 20 qubits máx, 300 trials Optuna"
Write-Host "Tempo estimado: 45 min (GPU), 120 min (CPU)"
Write-Host ""
$response = Read-Host "Executar EGFR? (s/n)"
if ($response -eq "s" -or $response -eq "S") {
    Write-Host "Iniciando EGFR..." -ForegroundColor Cyan
    python vqc_drug_tuner.py `
        --target EGFR `
        --max-qubits 20 `
        --trials 300 `
        --seed 42 `
        --out-dir results_vqc_drug `
        | Tee-Object -FilePath $LOG_FILE -Append
    Write-Host "✅ EGFR completo" -ForegroundColor Green
}

Write-Host ""

# ============================================================================
# EXEMPLO 2: HIV (Produção)
# ============================================================================
Write-Host "$(Get-Date) ========== EXEMPLO 2: HIV (Produção) ==========" -ForegroundColor Yellow
Write-Host "Target: HIV (41,913 moléculas - grande!)"
Write-Host "Configuração: 16 qubits máx, 200 trials"
Write-Host "Tempo estimado: 90 min (GPU), 240 min (CPU)"
Write-Host ""
$response = Read-Host "Executar HIV? (s/n)"
if ($response -eq "s" -or $response -eq "S") {
    Write-Host "Iniciando HIV..." -ForegroundColor Cyan
    python vqc_drug_tuner.py `
        --target HIV `
        --max-qubits 16 `
        --trials 200 `
        --seed 42 `
        --out-dir results_vqc_drug `
        | Tee-Object -FilePath $LOG_FILE -Append
    Write-Host "✅ HIV completo" -ForegroundColor Green
}

Write-Host ""

# ============================================================================
# EXEMPLO 3: Malaria (Rápido)
# ============================================================================
Write-Host "$(Get-Date) ========== EXEMPLO 3: Malaria (Rápido) ==========" -ForegroundColor Yellow
Write-Host "Target: Malaria (13,281 moléculas)"
Write-Host "Configuração: 12 qubits máx, 150 trials"
Write-Host "Tempo estimado: 30 min (GPU), 80 min (CPU)"
Write-Host ""
$response = Read-Host "Executar Malaria? (s/n)"
if ($response -eq "s" -or $response -eq "S") {
    Write-Host "Iniciando Malaria..." -ForegroundColor Cyan
    python vqc_drug_tuner.py `
        --target Malaria `
        --max-qubits 12 `
        --trials 150 `
        --seed 42 `
        --out-dir results_vqc_drug `
        | Tee-Object -FilePath $LOG_FILE -Append
    Write-Host "✅ Malaria completo" -ForegroundColor Green
}

Write-Host ""

# ============================================================================
# EXEMPLO 4: COVID (Real)
# ============================================================================
Write-Host "$(Get-Date) ========== EXEMPLO 4: COVID-19 (Real) ==========" -ForegroundColor Yellow
Write-Host "Target: COVID-19 (10,427 moléculas)"
Write-Host "Configuração: 14 qubits máx, 250 trials"
Write-Host "Tempo estimado: 40 min (GPU), 110 min (CPU)"
Write-Host ""
$response = Read-Host "Executar COVID? (s/n)"
if ($response -eq "s" -or $response -eq "S") {
    Write-Host "Iniciando COVID..." -ForegroundColor Cyan
    python vqc_drug_tuner.py `
        --target COVID `
        --max-qubits 14 `
        --trials 250 `
        --seed 42 `
        --out-dir results_vqc_drug `
        | Tee-Object -FilePath $LOG_FILE -Append
    Write-Host "✅ COVID completo" -ForegroundColor Green
}

Write-Host ""

# ============================================================================
# RESUMO FINAL
# ============================================================================
Write-Host "$(Get-Date) ========== RESUMO FINAL ==========" -ForegroundColor Yellow
Write-Host ""
Write-Host "📊 Resultados salvos em: results_vqc_drug/" -ForegroundColor Cyan
Write-Host "  - JSON reports (best params + métricas)"
Write-Host "  - Markdown reports (human-readable)"
Write-Host "  - Plotly HTML (gráficos interativos)"
Write-Host ""
Write-Host "📋 Cache QSAR em: qsar_cache/" -ForegroundColor Cyan
Write-Host "  - Download automático apenas na primeira execução"
Write-Host ""
Write-Host "📝 Log completo: $LOG_FILE" -ForegroundColor Cyan
Write-Host ""
Write-Host "✅ Pipeline completo!" -ForegroundColor Green
Write-Host ""

# List generated files
Write-Host "Arquivos de relatório gerados:" -ForegroundColor Cyan
Get-ChildItem results_vqc_drug -Filter "*.json", "*.md", "*.html" -ErrorAction SilentlyContinue | 
    ForEach-Object { Write-Host "  - $($_.Name) ($('{0:N0}' -f $_.Length) bytes)" }

Write-Host ""
Write-Host "Para visualizar gráficos interativos, abra:" -ForegroundColor Cyan
Write-Host "  results_vqc_drug/optuna_history.html"
Write-Host ""
