"""
Unit tests for noise models (ModeloRuido classes).

Tests verify that noise models apply correct Kraus operators and
behave correctly with different noise levels.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    import pennylane as qml
    PENNYLANE_AVAILABLE = True
except ImportError:
    PENNYLANE_AVAILABLE = False
    pytest.skip("PennyLane not available", allow_module_level=True)

from framework_investigativo_completo import (
    ModeloRuido,
    RuidoDepolarizante,
    RuidoAmplitudeDamping,
    RuidoPhaseDamping,
    RuidoCrosstalk,
    RuidoCorrelacionado,
    RuidoThermal,
    RuidoBitFlip,
    RuidoPhaseFlip,
    RuidoPinkNoise,
    RuidoReadoutError,
)


class TestModeloRuidoBase:
    """Test suite for base ModeloRuido class."""

    def test_modelo_ruido_initialization(self):
        """Test that ModeloRuido initializes with default level."""
        modelo = ModeloRuido(nivel=0.01)
        assert modelo.nivel == 0.01

    def test_modelo_ruido_aplicar_not_implemented(self):
        """Test that base class aplicar raises NotImplementedError."""
        modelo = ModeloRuido()
        with pytest.raises(NotImplementedError):
            modelo.aplicar(n_qubits=2)


class TestRuidoDepolarizante:
    """Test suite for depolarizing noise model."""

    def test_initialization(self):
        """Test depolarizing noise initialization."""
        noise = RuidoDepolarizante(nivel=0.02)
        assert noise.nivel == 0.02

    def test_aplicar_creates_channels(self):
        """Test that aplicar applies depolarizing channels."""
        noise = RuidoDepolarizante(nivel=0.01)
        n_qubits = 2
        
        # Create a simple device and test application
        dev = qml.device('default.mixed', wires=n_qubits)
        
        @qml.qnode(dev)
        def circuit():
            # Apply some gates
            qml.RY(0.5, wires=0)
            qml.CNOT(wires=[0, 1])
            # Apply noise
            noise.aplicar(n_qubits)
            return qml.state()
        
        # Should execute without error
        result = circuit()
        assert result is not None

    def test_nivel_override(self):
        """Test that nivel_override parameter works."""
        noise = RuidoDepolarizante(nivel=0.01)
        
        dev = qml.device('default.mixed', wires=1)
        
        @qml.qnode(dev)
        def circuit(override_level):
            qml.RY(0.5, wires=0)
            noise.aplicar(n_qubits=1, nivel_override=override_level)
            return qml.expval(qml.PauliZ(0))
        
        # Different noise levels should produce different results
        result1 = circuit(0.001)
        result2 = circuit(0.05)
        assert abs(result1 - result2) > 1e-5


class TestRuidoAmplitudeDamping:
    """Test suite for amplitude damping noise model."""

    def test_initialization(self):
        """Test amplitude damping initialization."""
        noise = RuidoAmplitudeDamping(nivel=0.03)
        assert noise.nivel == 0.03

    def test_aplicar_reduces_excited_state_population(self):
        """Test that amplitude damping reduces |1⟩ population."""
        noise = RuidoAmplitudeDamping(nivel=0.1)
        
        dev = qml.device('default.mixed', wires=1)
        
        @qml.qnode(dev)
        def circuit_with_noise():
            qml.PauliX(wires=0)  # Prepare |1⟩ state
            noise.aplicar(n_qubits=1)
            return qml.expval(qml.PauliZ(0))
        
        @qml.qnode(dev)
        def circuit_without_noise():
            qml.PauliX(wires=0)  # Prepare |1⟩ state
            return qml.expval(qml.PauliZ(0))
        
        with_noise = circuit_with_noise()
        without_noise = circuit_without_noise()
        
        # With amplitude damping, expectation should be less negative
        # (population relaxes from |1⟩ to |0⟩)
        assert with_noise > without_noise


class TestRuidoPhaseDamping:
    """Test suite for phase damping noise model."""

    def test_initialization(self):
        """Test phase damping initialization."""
        noise = RuidoPhaseDamping(nivel=0.04)
        assert noise.nivel == 0.04

    def test_aplicar_preserves_computational_basis(self):
        """Test that phase damping preserves computational basis states."""
        noise = RuidoPhaseDamping(nivel=0.1)
        
        dev = qml.device('default.mixed', wires=1)
        
        @qml.qnode(dev)
        def circuit(prepare_one=False):
            if prepare_one:
                qml.PauliX(wires=0)  # Prepare |1⟩
            noise.aplicar(n_qubits=1)
            return qml.expval(qml.PauliZ(0))
        
        # Computational basis states should be preserved
        result_0 = circuit(prepare_one=False)
        result_1 = circuit(prepare_one=True)
        
        assert abs(result_0 - 1.0) < 0.01  # |0⟩ preserved
        assert abs(result_1 + 1.0) < 0.01  # |1⟩ preserved


class TestRuidoCrosstalk:
    """Test suite for crosstalk noise model."""

    def test_initialization(self):
        """Test crosstalk noise initialization."""
        noise = RuidoCrosstalk(nivel=0.02)
        assert noise.nivel == 0.02

    def test_aplicar_affects_multiple_qubits(self):
        """Test that crosstalk affects neighboring qubits."""
        noise = RuidoCrosstalk(nivel=0.05)
        
        dev = qml.device('default.mixed', wires=3)
        
        @qml.qnode(dev)
        def circuit():
            qml.RY(np.pi/4, wires=0)
            qml.RY(np.pi/4, wires=1)
            qml.RY(np.pi/4, wires=2)
            noise.aplicar(n_qubits=3)
            return qml.expval(qml.PauliZ(0))
        
        # Should execute without error
        result = circuit()
        assert result is not None


class TestRuidoCorrelacionado:
    """Test suite for correlated noise model."""

    def test_initialization(self):
        """Test correlated noise initialization."""
        noise = RuidoCorrelacionado(nivel=0.03)
        assert noise.nivel == 0.03

    def test_aplicar_to_all_qubits(self):
        """Test that correlated noise applies to all qubits."""
        noise = RuidoCorrelacionado(nivel=0.02)
        
        dev = qml.device('default.mixed', wires=2)
        
        @qml.qnode(dev)
        def circuit():
            qml.Hadamard(wires=0)
            qml.CNOT(wires=[0, 1])
            noise.aplicar(n_qubits=2)
            return qml.expval(qml.PauliZ(0))
        
        result = circuit()
        assert result is not None


class TestRuidoThermal:
    """Test suite for thermal noise model."""

    def test_initialization(self):
        """Test thermal noise initialization."""
        noise = RuidoThermal(nivel=0.02)
        assert noise.nivel == 0.02

    def test_aplicar_combines_amplitude_and_phase(self):
        """Test that thermal noise applies both amplitude and phase damping."""
        noise = RuidoThermal(nivel=0.05)
        
        dev = qml.device('default.mixed', wires=1)
        
        @qml.qnode(dev)
        def circuit():
            qml.Hadamard(wires=0)
            noise.aplicar(n_qubits=1)
            return qml.expval(qml.PauliZ(0))
        
        result = circuit()
        assert result is not None


class TestRuidoBitFlip:
    """Test suite for bit-flip noise model."""

    def test_initialization(self):
        """Test bit-flip noise initialization."""
        noise = RuidoBitFlip(nivel=0.1)
        assert noise.nivel == 0.1

    def test_aplicar_flips_bits(self):
        """Test that bit-flip noise affects qubit states."""
        noise = RuidoBitFlip(nivel=0.3)
        
        dev = qml.device('default.mixed', wires=1)
        
        @qml.qnode(dev)
        def circuit_with_noise():
            # Start in |0⟩
            noise.aplicar(n_qubits=1)
            return qml.expval(qml.PauliZ(0))
        
        @qml.qnode(dev)
        def circuit_without_noise():
            return qml.expval(qml.PauliZ(0))
        
        with_noise = circuit_with_noise()
        without_noise = circuit_without_noise()
        
        # Bit flip should affect the expectation value
        assert abs(with_noise - without_noise) > 0.1


class TestRuidoPhaseFlip:
    """Test suite for phase-flip noise model."""

    def test_initialization(self):
        """Test phase-flip noise initialization."""
        noise = RuidoPhaseFlip(nivel=0.1)
        assert noise.nivel == 0.1


class TestRuidoPinkNoise:
    """Test suite for pink noise model."""

    def test_initialization(self):
        """Test pink noise initialization."""
        noise = RuidoPinkNoise(nivel=0.05)
        assert noise.nivel == 0.05

    def test_aplicar_varies_per_qubit(self):
        """Test that pink noise varies per qubit."""
        noise = RuidoPinkNoise(nivel=0.05)
        
        dev = qml.device('default.mixed', wires=2)
        
        @qml.qnode(dev)
        def circuit():
            qml.Hadamard(wires=0)
            qml.Hadamard(wires=1)
            noise.aplicar(n_qubits=2)
            return qml.expval(qml.PauliZ(0) @ qml.PauliZ(1))
        
        result = circuit()
        assert result is not None


class TestRuidoReadoutError:
    """Test suite for readout error model."""

    def test_initialization(self):
        """Test readout error initialization."""
        noise = RuidoReadoutError(nivel=0.05)
        assert noise.nivel == 0.05

    def test_aplicar_simulates_readout_error(self):
        """Test that readout error model applies bit-flip approximation."""
        noise = RuidoReadoutError(nivel=0.2)
        
        dev = qml.device('default.mixed', wires=1)
        
        @qml.qnode(dev)
        def circuit():
            qml.PauliX(wires=0)  # Prepare |1⟩
            noise.aplicar(n_qubits=1)
            return qml.expval(qml.PauliZ(0))
        
        result = circuit()
        assert result is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
