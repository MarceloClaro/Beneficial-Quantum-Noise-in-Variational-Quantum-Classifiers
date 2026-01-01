#!/bin/bash
# ============================================================================
# run_vqc_drug_examples.sh
# Execute multiple VQC-Molecular experiments with pre-configured parameters
# ============================================================================

set -e  # Exit on error

echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║     VQC-Molecular v8.0 - Drug Screening Pipeline Examples          ║"
echo "║  Quantum-Enhanced QSAR with Automatic Hyper-parameter Tuning       ║"
echo "╚════════════════════════════════════════════════════════════════════╝"

# Create results directory
mkdir -p results_vqc_drug
mkdir -p qsar_cache
mkdir -p logs

TIMESTAMP=$(date +%Y-%m-%d_%H-%M-%S)
LOG_FILE="logs/vqc_drug_${TIMESTAMP}.log"

{

echo ""
echo "[$(date)] ========== Verificando Dependências =========="
python3 -c "import pennylane as qml; print(f'✓ PennyLane {qml.__version__}')"
python3 -c "import optuna; print(f'✓ Optuna {optuna.__version__}')"
python3 -c "import deepchem as dc; print(f'✓ DeepChem {dc.__version__}')"
python3 -c "import rdkit; print(f'✓ RDKit OK')"
echo ""

# ============================================================================
# EXEMPLO 1: EGFR (Piloto - 20 qubits, 300 trials)
# ============================================================================
echo "[$(date)] ========== EXEMPLO 1: EGFR Kinase (Piloto) =========="
echo "Target: EGFR (6,847 moléculas)"
echo "Configuração: 20 qubits máx, 300 trials Optuna"
echo "Tempo estimado: 45 min (GPU), 120 min (CPU)"
echo ""
read -p "Executar EGFR? (s/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Ss]$ ]]; then
    python3 vqc_drug_tuner.py \
        --target EGFR \
        --max-qubits 20 \
        --trials 300 \
        --seed 42 \
        --out-dir results_vqc_drug
    echo "[$(date)] ✅ EGFR completo"
fi

echo ""

# ============================================================================
# EXEMPLO 2: HIV (Produção - 16 qubits, 200 trials)
# ============================================================================
echo "[$(date)] ========== EXEMPLO 2: HIV (Produção) =========="
echo "Target: HIV (41,913 moléculas - grande!)"
echo "Configuração: 16 qubits máx, 200 trials"
echo "Tempo estimado: 90 min (GPU), 240 min (CPU)"
echo ""
read -p "Executar HIV? (s/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Ss]$ ]]; then
    python3 vqc_drug_tuner.py \
        --target HIV \
        --max-qubits 16 \
        --trials 200 \
        --seed 42 \
        --out-dir results_vqc_drug
    echo "[$(date)] ✅ HIV completo"
fi

echo ""

# ============================================================================
# EXEMPLO 3: Malaria (Rápido - 12 qubits, 150 trials)
# ============================================================================
echo "[$(date)] ========== EXEMPLO 3: Malaria (Rápido) =========="
echo "Target: Malaria (13,281 moléculas)"
echo "Configuração: 12 qubits máx, 150 trials"
echo "Tempo estimado: 30 min (GPU), 80 min (CPU)"
echo ""
read -p "Executar Malaria? (s/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Ss]$ ]]; then
    python3 vqc_drug_tuner.py \
        --target Malaria \
        --max-qubits 12 \
        --trials 150 \
        --seed 42 \
        --out-dir results_vqc_drug
    echo "[$(date)] ✅ Malaria completo"
fi

echo ""

# ============================================================================
# EXEMPLO 4: COVID (Real - 14 qubits, 250 trials)
# ============================================================================
echo "[$(date)] ========== EXEMPLO 4: COVID-19 (Real) =========="
echo "Target: COVID-19 (10,427 moléculas)"
echo "Configuração: 14 qubits máx, 250 trials"
echo "Tempo estimado: 40 min (GPU), 110 min (CPU)"
echo ""
read -p "Executar COVID? (s/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Ss]$ ]]; then
    python3 vqc_drug_tuner.py \
        --target COVID \
        --max-qubits 14 \
        --trials 250 \
        --seed 42 \
        --out-dir results_vqc_drug
    echo "[$(date)] ✅ COVID completo"
fi

echo ""

# ============================================================================
# RESUMO FINAL
# ============================================================================
echo "[$(date)] ========== RESUMO FINAL =========="
echo ""
echo "📊 Resultados salvos em: results_vqc_drug/"
echo "  - JSON reports (best params + métricas)"
echo "  - Markdown reports (human-readable)"
echo "  - Plotly HTML (gráficos interativos)"
echo ""
echo "📋 Cache QSAR em: qsar_cache/"
echo "  - Download automático apenas na primeira execução"
echo ""
echo "📝 Log completo: $LOG_FILE"
echo ""
echo "✅ Pipeline completo!"
echo ""

} | tee "$LOG_FILE"

# Print summary
echo ""
echo "Arquivos de relatório gerados:"
ls -lh results_vqc_drug/*.json results_vqc_drug/*.md results_vqc_drug/*.html 2>/dev/null || echo "  (Nenhum ainda)"
echo ""
echo "Para visualizar gráficos interativos, abra:"
echo "  results_vqc_drug/optuna_history.html"
echo ""
