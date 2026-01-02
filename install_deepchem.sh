#!/bin/bash
# =============================================================================
# DeepChem Installation Script
# =============================================================================
# Install DeepChem and dependencies for molecular dataset integration
#
# Usage:
#   bash install_deepchem.sh [python_version] [device_type]
#
# Arguments:
#   python_version: Python version (default: 3.10)
#   device_type: 'cpu' or 'gpu' (default: cpu)
#
# Example:
#   bash install_deepchem.sh 3.10 cpu
#
# Documentation: https://deepchem.readthedocs.io/en/latest/get_started/installation.html

set -e  # Exit on error

PYTHON_VERSION=${1:-3.10}
DEVICE_TYPE=${2:-cpu}

echo "============================================"
echo "DeepChem Installation Script"
echo "============================================"
echo "Python Version: $PYTHON_VERSION"
echo "Device Type: $DEVICE_TYPE"
echo ""

# Check if conda is available
if command -v conda &> /dev/null; then
    echo "✓ Conda found"
    USE_CONDA=true
else
    echo "✗ Conda not found, using pip"
    USE_CONDA=false
fi

# Method 1: Using Conda (Recommended)
if [ "$USE_CONDA" = true ]; then
    echo ""
    echo "Installing DeepChem via Conda..."
    echo ""
    
    # Create conda environment
    ENV_NAME="deepchem_env"
    echo "Creating conda environment: $ENV_NAME"
    conda create -n $ENV_NAME python=$PYTHON_VERSION -y
    
    # Activate environment
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate $ENV_NAME
    
    # Install DeepChem
    if [ "$DEVICE_TYPE" = "gpu" ]; then
        echo "Installing DeepChem with GPU support..."
        conda install -c conda-forge deepchem-gpu -y
    else
        echo "Installing DeepChem with CPU support..."
        conda install -c conda-forge deepchem -y
    fi
    
    # Install RDKit (required for molecular featurization)
    echo "Installing RDKit..."
    conda install -c conda-forge rdkit -y
    
    # Install additional dependencies
    echo "Installing additional dependencies..."
    conda install -c conda-forge scikit-learn pandas numpy scipy matplotlib -y
    
    echo ""
    echo "✓ DeepChem installed via Conda"
    echo ""
    echo "To activate the environment, run:"
    echo "  conda activate $ENV_NAME"
    
else
    # Method 2: Using pip
    echo ""
    echo "Installing DeepChem via pip..."
    echo ""
    
    # Upgrade pip
    python -m pip install --upgrade pip
    
    # Install DeepChem
    if [ "$DEVICE_TYPE" = "gpu" ]; then
        echo "Installing DeepChem with GPU support..."
        pip install deepchem[torch,tensorflow]
    else
        echo "Installing DeepChem..."
        pip install deepchem
    fi
    
    # Note about RDKit
    echo ""
    echo "⚠️  Note: RDKit installation via pip can be problematic."
    echo "    If you encounter issues, consider using conda:"
    echo "    conda install -c conda-forge rdkit"
    echo ""
    
    # Try to install RDKit via pip (may not work on all systems)
    echo "Attempting to install RDKit via pip..."
    pip install rdkit-pypi || echo "⚠️  RDKit pip install failed. Use conda if needed."
    
    echo ""
    echo "✓ DeepChem installed via pip"
fi

# Verify installation
echo ""
echo "============================================"
echo "Verifying DeepChem Installation"
echo "============================================"
echo ""

python << 'VERIFY_EOF'
try:
    import deepchem as dc
    print(f"✓ DeepChem version: {dc.__version__}")
    
    import rdkit
    print(f"✓ RDKit version: {rdkit.__version__}")
    
    # Test loading a simple dataset
    print("\nTesting dataset loading...")
    tasks, datasets, transformers = dc.molnet.load_bace_classification(
        featurizer='ECFP',
        splitter='random',
        reload=False
    )
    print(f"✓ BACE dataset loaded: {len(datasets[0])} training samples")
    
    print("\n✓ DeepChem installation verified successfully!")
    print("\nAvailable MoleculeNet datasets:")
    print("  - BACE: β-secretase inhibition")
    print("  - HIV: Anti-HIV activity")
    print("  - Tox21: Toxicity prediction")
    print("  - BBBP: Blood-brain barrier penetration")
    print("  - and many more...")
    
except ImportError as e:
    print(f"✗ Import error: {e}")
    print("\nTroubleshooting:")
    print("  1. Try installing with conda if you used pip")
    print("  2. Check Python version compatibility")
    print("  3. See: https://deepchem.readthedocs.io/en/latest/get_started/installation.html")
    exit(1)
except Exception as e:
    print(f"✗ Verification error: {e}")
    exit(1)
VERIFY_EOF

echo ""
echo "============================================"
echo "Installation Complete!"
echo "============================================"
echo ""
echo "Next steps:"
echo "  1. Run the advanced framework:"
echo "     python framework_quantum_advanced_v8.py"
echo ""
echo "  2. Check documentation:"
echo "     https://deepchem.readthedocs.io/en/latest/"
echo ""
echo "  3. Explore examples:"
echo "     https://deepchem.readthedocs.io/en/latest/get_started/examples.html"
echo ""
