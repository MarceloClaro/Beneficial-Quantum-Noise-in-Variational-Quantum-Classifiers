#!/bin/bash
# Multiframework Validation Run - CPU Optimized
# Estimated runtime: 2-3 hours
# Configurations: 432 (4 datasets × 3 noise × 3 ansätze × 2 schedules × 3 frameworks × 2 seeds)

set -e  # Exit on error

echo "============================================================================="
echo "MULTIFRAMEWORK VALIDATION RUN - CPU OPTIMIZED"
echo "============================================================================="
echo ""
echo "Configuration: configs/experiment_cpu_optimized.yaml"
echo "Estimated time: 2-3 hours"
echo "Output: results/*/20251227_001/"
echo ""

# Check dependencies
echo "Checking dependencies..."
python -c "import pennylane, qiskit, cirq; print('✓ All frameworks installed')" || {
    echo "✗ Missing dependencies. Run: pip install -r requirements.txt"
    exit 1
}

# Create output directories
RUN_ID="20251227_001"
echo "Creating output directories for run ID: $RUN_ID"
mkdir -p results/pennylane/$RUN_ID
mkdir -p results/qiskit/$RUN_ID
mkdir -p results/cirq/$RUN_ID
mkdir -p results/comparisons/$RUN_ID
mkdir -p logs/pennylane/$RUN_ID
mkdir -p logs/qiskit/$RUN_ID
mkdir -p logs/cirq/$RUN_ID
mkdir -p manifests/pennylane/$RUN_ID
mkdir -p manifests/qiskit/$RUN_ID
mkdir -p manifests/cirq/$RUN_ID

# CPU Optimization
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export NUMEXPR_NUM_THREADS=4

echo "✓ Environment configured for CPU execution (4 threads)"
echo ""

# Phase 5: PennyLane Baseline
echo "============================================================================="
echo "PHASE 5: PennyLane Baseline Execution"
echo "============================================================================="
echo "Estimated time: 40-60 minutes"
echo ""

START_TIME=$(date +%s)
python framework_investigativo_completo.py \
    --config configs/experiment_cpu_optimized.yaml \
    --output results/pennylane/$RUN_ID \
    --run-id $RUN_ID \
    2>&1 | tee logs/pennylane/$RUN_ID/stdout.log

PENNYLANE_TIME=$(($(date +%s) - START_TIME))
echo "✓ PennyLane completed in $PENNYLANE_TIME seconds"
echo ""

# Phase 6: Qiskit Execution
echo "============================================================================="
echo "PHASE 6: Qiskit Execution"
echo "============================================================================="
echo "Estimated time: 60-90 minutes"
echo ""

START_TIME=$(date +%s)
python executar_framework_qiskit.py \
    --config configs/experiment_cpu_optimized.yaml \
    --output results/qiskit/$RUN_ID \
    --run-id $RUN_ID \
    2>&1 | tee logs/qiskit/$RUN_ID/stdout.log

QISKIT_TIME=$(($(date +%s) - START_TIME))
echo "✓ Qiskit completed in $QISKIT_TIME seconds"
echo ""

# Phase 7: Cirq Execution
echo "============================================================================="
echo "PHASE 7: Cirq Execution"
echo "============================================================================="
echo "Estimated time: 50-70 minutes"
echo ""

START_TIME=$(date +%s)
python executar_framework_cirq.py \
    --config configs/experiment_cpu_optimized.yaml \
    --output results/cirq/$RUN_ID \
    --run-id $RUN_ID \
    2>&1 | tee logs/cirq/$RUN_ID/stdout.log

CIRQ_TIME=$(($(date +%s) - START_TIME))
echo "✓ Cirq completed in $CIRQ_TIME seconds"
echo ""

# Phase 8: Comparative Analysis
echo "============================================================================="
echo "PHASE 8: Comparative Analysis"
echo "============================================================================="
echo "Estimated time: 10-15 minutes"
echo ""

START_TIME=$(date +%s)
python generate_comparative_results.py \
    --run-id $RUN_ID \
    --mode validation \
    2>&1 | tee results/comparisons/$RUN_ID/analysis.log

ANALYSIS_TIME=$(($(date +%s) - START_TIME))
echo "✓ Analysis completed in $ANALYSIS_TIME seconds"
echo ""

# Phase 9: Figure Generation
echo "============================================================================="
echo "PHASE 9: Figure Generation"
echo "============================================================================="
echo "Estimated time: 5-10 minutes"
echo ""

if [ -f "tools/generate_all_figures.py" ]; then
    START_TIME=$(date +%s)
    python tools/generate_all_figures.py \
        --run-id $RUN_ID \
        2>&1 | tee figures/$RUN_ID/generation.log
    
    FIGURES_TIME=$(($(date +%s) - START_TIME))
    echo "✓ Figures completed in $FIGURES_TIME seconds"
else
    echo "⚠ Figure generation script not found, skipping..."
    FIGURES_TIME=0
fi
echo ""

# Summary
TOTAL_TIME=$((PENNYLANE_TIME + QISKIT_TIME + CIRQ_TIME + ANALYSIS_TIME + FIGURES_TIME))
echo "============================================================================="
echo "VALIDATION RUN COMPLETE"
echo "============================================================================="
echo ""
echo "Timing Summary:"
echo "  PennyLane:  ${PENNYLANE_TIME}s ($(($PENNYLANE_TIME / 60))m)"
echo "  Qiskit:     ${QISKIT_TIME}s ($(($QISKIT_TIME / 60))m)"
echo "  Cirq:       ${CIRQ_TIME}s ($(($CIRQ_TIME / 60))m)"
echo "  Analysis:   ${ANALYSIS_TIME}s ($(($ANALYSIS_TIME / 60))m)"
echo "  Figures:    ${FIGURES_TIME}s ($(($FIGURES_TIME / 60))m)"
echo "  TOTAL:      ${TOTAL_TIME}s ($(($TOTAL_TIME / 60))m)"
echo ""
echo "Results saved to:"
echo "  - results/pennylane/$RUN_ID/"
echo "  - results/qiskit/$RUN_ID/"
echo "  - results/cirq/$RUN_ID/"
echo "  - results/comparisons/$RUN_ID/"
echo ""
echo "Next steps:"
echo "  1. Verify results: ls -lh results/*/20251227_001/"
echo "  2. Check logs: less logs/*/20251227_001/stdout.log"
echo "  3. Generate article updates: python atualizar_artigos_com_resultados.py --run-id $RUN_ID"
echo ""
