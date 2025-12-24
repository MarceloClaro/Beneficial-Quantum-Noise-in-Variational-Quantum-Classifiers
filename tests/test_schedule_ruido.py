"""
Unit tests for ScheduleRuido class.

Tests verify that noise scheduling (annealing) curves work correctly
with different schedules: linear, exponential, cosine, and adaptive.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from framework_investigativo_completo import ScheduleRuido


class TestScheduleRuido:
    """Test suite for noise scheduling."""

    def test_linear_schedule(self):
        """Test linear noise schedule."""
        schedule = ScheduleRuido(tipo='linear', nivel_inicial=0.05, nivel_final=0.01)
        
        # Test at different epochs
        nivel_0 = schedule.calcular_nivel(epoca=0, n_epocas=100)
        nivel_50 = schedule.calcular_nivel(epoca=50, n_epocas=100)
        nivel_100 = schedule.calcular_nivel(epoca=99, n_epocas=100)
        
        # Check that noise decreases linearly
        assert abs(nivel_0 - 0.05) < 1e-6
        assert abs(nivel_100 - 0.01) < 1e-6
        assert 0.01 < nivel_50 < 0.05
        
        # Check linearity
        expected_50 = 0.05 - (0.05 - 0.01) * 50 / 99
        assert abs(nivel_50 - expected_50) < 1e-6

    def test_exponential_schedule(self):
        """Test exponential noise schedule."""
        schedule = ScheduleRuido(tipo='exponencial', nivel_inicial=0.05, nivel_final=0.01)
        
        nivel_0 = schedule.calcular_nivel(epoca=0, n_epocas=100)
        nivel_50 = schedule.calcular_nivel(epoca=50, n_epocas=100)
        nivel_100 = schedule.calcular_nivel(epoca=99, n_epocas=100)
        
        # Check that noise decreases exponentially
        assert abs(nivel_0 - 0.05) < 1e-6
        # Final value should be close to nivel_final
        assert nivel_100 < nivel_50 < nivel_0

    def test_cosine_schedule(self):
        """Test cosine annealing schedule."""
        schedule = ScheduleRuido(tipo='cosine', nivel_inicial=0.05, nivel_final=0.01)
        
        nivel_0 = schedule.calcular_nivel(epoca=0, n_epocas=100)
        nivel_50 = schedule.calcular_nivel(epoca=50, n_epocas=100)
        nivel_100 = schedule.calcular_nivel(epoca=99, n_epocas=100)
        
        # Check cosine annealing properties
        assert abs(nivel_0 - 0.05) < 1e-6
        assert abs(nivel_100 - 0.01) < 1e-6
        assert 0.01 < nivel_50 < 0.05

    def test_constant_schedule(self):
        """Test constant noise level (no scheduling)."""
        schedule = ScheduleRuido(tipo='constante', nivel_inicial=0.03, nivel_final=0.03)
        
        nivel_0 = schedule.calcular_nivel(epoca=0, n_epocas=100)
        nivel_50 = schedule.calcular_nivel(epoca=50, n_epocas=100)
        nivel_100 = schedule.calcular_nivel(epoca=99, n_epocas=100)
        
        # All should be the same
        assert abs(nivel_0 - 0.03) < 1e-6
        assert abs(nivel_50 - 0.03) < 1e-6
        assert abs(nivel_100 - 0.03) < 1e-6

    def test_schedule_decreasing(self):
        """Test that all schedules decrease monotonically."""
        schedules = ['linear', 'exponencial', 'cosine']
        
        for tipo in schedules:
            schedule = ScheduleRuido(tipo=tipo, nivel_inicial=0.1, nivel_final=0.01)
            
            niveis = [schedule.calcular_nivel(epoca=e, n_epocas=100) 
                     for e in range(100)]
            
            # Check monotonicity (allowing for small numerical errors)
            for i in range(len(niveis) - 1):
                assert niveis[i] >= niveis[i+1] - 1e-10, \
                    f"Schedule {tipo} not monotonically decreasing at epoch {i}"

    def test_schedule_bounds(self):
        """Test that schedule values stay within bounds."""
        schedule = ScheduleRuido(tipo='linear', nivel_inicial=0.1, nivel_final=0.01)
        
        for epoca in range(150):
            nivel = schedule.calcular_nivel(epoca=epoca, n_epocas=100)
            assert 0.0 <= nivel <= 0.15, \
                f"Noise level {nivel} out of bounds at epoch {epoca}"

    def test_adaptive_schedule(self):
        """Test adaptive noise schedule based on gradient information."""
        schedule = ScheduleRuido(tipo='adaptativo', nivel_inicial=0.05, nivel_final=0.01)
        
        # Simulate gradient variance (high variance should increase noise)
        nivel_1 = schedule.calcular_nivel_adaptativo(
            epoca=10, n_epocas=100, variancia_gradiente=1e-5
        )
        nivel_2 = schedule.calcular_nivel_adaptativo(
            epoca=10, n_epocas=100, variancia_gradiente=1e-8
        )
        
        # Higher gradient variance should lead to higher noise
        # (to help escape barren plateaus)
        assert nivel_1 > nivel_2 or abs(nivel_1 - nivel_2) < 1e-6

    def test_schedule_edge_cases(self):
        """Test edge cases for scheduling."""
        schedule = ScheduleRuido(tipo='linear', nivel_inicial=0.05, nivel_final=0.01)
        
        # Test with single epoch
        nivel = schedule.calcular_nivel(epoca=0, n_epocas=1)
        assert abs(nivel - 0.05) < 1e-6
        
        # Test with zero epochs (edge case)
        nivel = schedule.calcular_nivel(epoca=0, n_epocas=0)
        assert nivel >= 0.0

    def test_invalid_schedule_type(self):
        """Test behavior with invalid schedule type."""
        schedule = ScheduleRuido(tipo='invalid_type', nivel_inicial=0.05, nivel_final=0.01)
        
        # Should fall back to constant or handle gracefully
        nivel = schedule.calcular_nivel(epoca=50, n_epocas=100)
        assert 0.0 <= nivel <= 0.1

    def test_schedule_reproducibility(self):
        """Test that schedules produce reproducible results."""
        schedule = ScheduleRuido(tipo='linear', nivel_inicial=0.05, nivel_final=0.01)
        
        nivel_1 = schedule.calcular_nivel(epoca=25, n_epocas=100)
        nivel_2 = schedule.calcular_nivel(epoca=25, n_epocas=100)
        
        assert nivel_1 == nivel_2

    def test_schedule_smooth_transition(self):
        """Test that schedules provide smooth transitions."""
        schedule = ScheduleRuido(tipo='cosine', nivel_inicial=0.1, nivel_final=0.01)
        
        niveis = [schedule.calcular_nivel(epoca=e, n_epocas=100) 
                 for e in range(100)]
        
        # Check that no two consecutive values differ by more than expected
        max_diff = max(abs(niveis[i] - niveis[i+1]) for i in range(len(niveis)-1))
        assert max_diff < 0.01, "Schedule has abrupt transitions"

    def test_schedule_with_warmup(self):
        """Test that schedule can start from low noise and increase."""
        # Reverse schedule: start low, end high (warmup strategy)
        schedule = ScheduleRuido(tipo='linear', nivel_inicial=0.01, nivel_final=0.05)
        
        nivel_0 = schedule.calcular_nivel(epoca=0, n_epocas=100)
        nivel_100 = schedule.calcular_nivel(epoca=99, n_epocas=100)
        
        # Should increase
        assert nivel_0 < nivel_100


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
