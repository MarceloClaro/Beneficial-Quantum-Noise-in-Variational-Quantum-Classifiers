"""
Pytest configuration file for quantum framework tests.

This file configures pytest behavior and defines shared fixtures
for all tests in the test suite.
"""

import sys
import os
from pathlib import Path

# Add parent directory to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
import warnings

# Filter out specific deprecation warnings from quantum libraries
# Only suppress known warnings from quantum frameworks
warnings.filterwarnings("ignore", category=DeprecationWarning, module="pennylane.*")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="qiskit.*")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="cirq.*")
warnings.filterwarnings("ignore", category=PendingDeprecationWarning, module="qiskit.*")


@pytest.fixture(scope="session", autouse=True)
def setup_environment():
    """Set up test environment before running any tests."""
    # Set environment variables for reproducibility
    os.environ["PYTHONHASHSEED"] = "0"
    
    # Disable GPU usage for consistent CPU-only testing
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
    
    yield
    
    # Cleanup after all tests


@pytest.fixture(scope="function")
def numpy_seed():
    """Fixture to set a consistent numpy seed for tests."""
    import numpy as np
    np.random.seed(42)
    return 42


def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "gpu: marks tests that require GPU (deselect with '-m \"not gpu\"')"
    )
