"""
Unit tests for ConstantesFundamentais class.

Tests verify that fundamental constants have correct numerical values
and that initialization strategies work correctly.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from framework_investigativo_completo import ConstantesFundamentais


class TestConstantesFundamentais:
    """Test suite for fundamental constants."""

    def test_pi_value(self):
        """Test that π constant has correct value."""
        assert abs(ConstantesFundamentais.PI - 3.14159265358979) < 1e-10

    def test_e_value(self):
        """Test that Euler's number has correct value."""
        assert abs(ConstantesFundamentais.E - 2.71828182845905) < 1e-10

    def test_golden_ratio_value(self):
        """Test that golden ratio φ has correct value."""
        expected = (1 + np.sqrt(5)) / 2
        assert abs(ConstantesFundamentais.PHI - expected) < 1e-10

    def test_hbar_value(self):
        """Test that reduced Planck constant ℏ has correct value."""
        # Check order of magnitude (should be ~1e-34)
        assert 1e-35 < ConstantesFundamentais.HBAR < 1e-33

    def test_fine_structure_value(self):
        """Test that fine structure constant α has correct value."""
        # Check it's approximately 1/137
        assert 0.007 < ConstantesFundamentais.ALPHA < 0.008

    def test_rydberg_value(self):
        """Test that Rydberg constant has correct value."""
        # Check order of magnitude (should be ~1e7)
        assert 1e7 < ConstantesFundamentais.RYDBERG < 1e8

    def test_inicializar_aleatorio(self):
        """Test random initialization strategy."""
        n_params = 10
        weights = ConstantesFundamentais.inicializar(
            n_params, estrategia='aleatoria', seed=42
        )
        assert weights.shape == (n_params,)
        # Random initialization can be in any range, just check finite
        assert np.all(np.isfinite(weights))

    def test_inicializar_matematico(self):
        """Test mathematical constants initialization."""
        n_params = 12
        weights = ConstantesFundamentais.inicializar(
            n_params, estrategia='matematico', seed=42
        )
        assert weights.shape == (n_params,)
        # Check that some values are close to mathematical constants
        assert any(abs(w - ConstantesFundamentais.PI) < 0.5 for w in weights)

    def test_inicializar_quantico(self):
        """Test quantum constants initialization."""
        n_params = 12
        weights = ConstantesFundamentais.inicializar(
            n_params, estrategia='quantico', seed=42
        )
        assert weights.shape == (n_params,)
        # Quantum constants are scaled appropriately
        assert np.all(np.isfinite(weights))

    def test_inicializar_fibonacci_spiral(self):
        """Test Fibonacci spiral initialization."""
        n_params = 10
        weights = ConstantesFundamentais.inicializar(
            n_params, estrategia='fibonacci_spiral', seed=42
        )
        assert weights.shape == (n_params,)
        # Check values are finite
        assert np.all(np.isfinite(weights))

    def test_inicializar_xavier_quantico(self):
        """Test Xavier quantum initialization."""
        n_params = 10
        weights = ConstantesFundamentais.inicializar(
            n_params, estrategia='xavier_quantico', seed=42
        )
        assert weights.shape == (n_params,)
        # Xavier should produce normally distributed values
        assert np.all(np.isfinite(weights))

    def test_inicializar_reproducibility(self):
        """Test that initialization is reproducible with same seed."""
        n_params = 10
        weights1 = ConstantesFundamentais.inicializar(
            n_params, estrategia='aleatoria', seed=42
        )
        weights2 = ConstantesFundamentais.inicializar(
            n_params, estrategia='aleatoria', seed=42
        )
        np.testing.assert_array_equal(weights1, weights2)

    def test_inicializar_different_seeds(self):
        """Test that different seeds produce different results."""
        n_params = 10
        weights1 = ConstantesFundamentais.inicializar(
            n_params, estrategia='aleatoria', seed=42
        )
        weights2 = ConstantesFundamentais.inicializar(
            n_params, estrategia='aleatoria', seed=43
        )
        assert not np.array_equal(weights1, weights2)

    def test_inicializar_invalid_strategy(self):
        """Test that invalid strategy falls back to random."""
        n_params = 10
        weights = ConstantesFundamentais.inicializar(
            n_params, estrategia='invalid_strategy', seed=42
        )
        assert weights.shape == (n_params,)
        # Should still produce valid weights
        assert np.all(np.isfinite(weights))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
