"""
Smoke tests for Beneficial Quantum Noise in Variational Quantum Classifiers

This test suite validates basic functionality and ensures the repository
is in a working state for Qualis A1 publication standards.
"""

import os
import sys
from pathlib import Path

import pytest

try:
    import torch
    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False

# Add the parent directory to the path so we can import the main module
sys.path.insert(0, str(Path(__file__).parent.parent))


def test_imports():
    """Test that all required imports work correctly."""
    try:
        import numpy as np
        import pandas as pd
        import pennylane as qml
        from sklearn.datasets import make_moons
        import plotly.graph_objects as go
        assert True
    except ImportError as e:
        pytest.fail(f"Import failed: {e}")


def test_repository_structure():
    """Test that the repository has the expected structure."""
    repo_root = Path(__file__).parent.parent
    required_files = [
        "README.md",
        "requirements.txt",
        "framework_investigativo_completo.py",
        "ANALISE_QUALIS_A1.md",
        "INSTALL.md",
        "LICENSE",
    ]

    for file in required_files:
        assert (repo_root / file).exists(), f"Missing required file: {file}"


def test_required_directories():
    """Test that all required directories exist."""
    repo_root = Path(__file__).parent.parent
    required_dirs = ["docs", "examples", "tests", "tools"]

    for dir_name in required_dirs:
        assert (repo_root / dir_name).exists(), f"Missing required directory: {dir_name}"


def test_documentation_files():
    """Test that all documentation files exist and are not empty."""
    repo_root = Path(__file__).parent.parent
    doc_files = [
        "README.md",
        "ANALISE_QUALIS_A1.md",
        "INSTALL.md",
        "STRUCTURE.md",
        "ORGANIZATION_CHECKLIST.md",
        "OBJETIVOS_PROJETO.md",
        "TEMPO_EXPERIMENTO.md",
    ]

    for doc_file in doc_files:
        file_path = repo_root / doc_file
        assert file_path.exists(), f"Missing documentation file: {doc_file}"
        assert file_path.stat().st_size > 0, f"Documentation file is empty: {doc_file}"


def test_requirements_file():
    """Test that requirements.txt exists and contains required packages."""
    repo_root = Path(__file__).parent.parent
    requirements_file = repo_root / "requirements.txt"

    assert requirements_file.exists(), "requirements.txt not found"

    with open(requirements_file) as f:
        requirements = f.read()

    required_packages = [
        "pennylane",
        "numpy",
        "pandas",
        "scipy",
        "scikit-learn",
        "plotly",
        "matplotlib",
        "statsmodels",
    ]

    for package in required_packages:
        assert package in requirements.lower(), f"Missing required package: {package}"


def test_framework_script_syntax():
    """Test that the main framework script has valid Python syntax."""
    repo_root = Path(__file__).parent.parent
    framework_script = repo_root / "framework_investigativo_completo.py"

    assert framework_script.exists(), "Main framework script not found"

    # Try to compile the script to check for syntax errors
    with open(framework_script) as f:
        code = f.read()

    try:
        compile(code, str(framework_script), "exec")
    except SyntaxError as e:
        pytest.fail(f"Syntax error in framework script: {e}")


def test_example_scripts():
    """Test that example scripts exist and have valid syntax."""
    repo_root = Path(__file__).parent.parent
    examples_dir = repo_root / "examples"

    example_files = list(examples_dir.glob("*.py"))
    assert len(example_files) > 0, "No example scripts found"

    for example_file in example_files:
        with open(example_file) as f:
            code = f.read()

        try:
            compile(code, str(example_file), "exec")
        except SyntaxError as e:
            pytest.fail(f"Syntax error in {example_file.name}: {e}")


def test_tool_scripts():
    """Test that tool scripts exist and have valid syntax."""
    repo_root = Path(__file__).parent.parent
    tools_dir = repo_root / "tools"

    tool_files = list(tools_dir.glob("*.py"))
    assert len(tool_files) > 0, "No tool scripts found"

    for tool_file in tool_files:
        with open(tool_file) as f:
            code = f.read()

        try:
            compile(code, str(tool_file), "exec")
        except SyntaxError as e:
            pytest.fail(f"Syntax error in {tool_file.name}: {e}")


def test_ruff_configuration():
    """Test that ruff configuration file is valid."""
    repo_root = Path(__file__).parent.parent
    ruff_config = repo_root / ".ruff.toml"

    assert ruff_config.exists(), ".ruff.toml not found"


def test_pennylane_basic_functionality():
    """Test basic PennyLane functionality."""
    try:
        import pennylane as qml
        import numpy as np

        # Create a simple quantum circuit
        dev = qml.device("default.qubit", wires=2)

        @qml.qnode(dev)
        def circuit(params):
            qml.RY(params[0], wires=0)
            qml.CNOT(wires=[0, 1])
            return qml.expval(qml.PauliZ(0))

        params = np.array([0.5])
        result = circuit(params)

        # PennyLane can return torch.Tensor if torch is available
        valid_types = (float, np.floating)
        if HAS_TORCH:
            valid_types = (float, np.floating, torch.Tensor)
        assert isinstance(result, valid_types)
    except Exception as e:
        pytest.fail(f"PennyLane basic functionality test failed: {e}")


def test_dataset_loading():
    """Test that scikit-learn datasets can be loaded."""
    try:
        from sklearn.datasets import make_moons, make_circles

        # Test make_moons
        X, y = make_moons(n_samples=100, noise=0.1, random_state=42)
        assert X.shape == (100, 2)
        assert y.shape == (100,)

        # Test make_circles
        X, y = make_circles(n_samples=100, noise=0.1, random_state=42)
        assert X.shape == (100, 2)
        assert y.shape == (100,)
    except Exception as e:
        pytest.fail(f"Dataset loading test failed: {e}")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
