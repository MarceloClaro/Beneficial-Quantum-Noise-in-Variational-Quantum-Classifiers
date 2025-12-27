#!/bin/bash
# Multiframework Minimal Test - CPU Optimized
# Estimated runtime: 10-15 minutes
# Quick test to verify infrastructure before full validation run

set -e

echo "============================================================================="
echo "MULTIFRAMEWORK MINIMAL TEST - CPU OPTIMIZED"
echo "============================================================================="
echo ""
echo "This is a quick test (10-15 min) to verify the infrastructure."
echo "Use run_multiframework_validation.sh for the full validation run (2-3h)."
echo ""

# Check dependencies
echo "Checking dependencies..."
python -c "import pennylane, qiskit, cirq; print('✓ All frameworks installed')" || {
    echo "✗ Missing dependencies. Run: pip install -r requirements.txt"
    exit 1
}

# Run smoke tests
echo ""
echo "Running smoke tests..."
python tests/smoke_multiframework.py || {
    echo "✗ Smoke tests failed. Check the output above."
    exit 1
}

echo ""
echo "✓ Infrastructure validated successfully!"
echo ""
echo "Next steps:"
echo "  1. Run validation: bash scripts/run_multiframework_validation.sh"
echo "  2. Or run full production: bash scripts/run_multiframework_full.sh"
echo ""
