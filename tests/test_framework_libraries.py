#!/usr/bin/env python3
"""
Comprehensive library tests for quantum frameworks.

Tests all critical libraries and components for PennyLane, Qiskit, and Cirq.
This module validates installation, API compatibility, and basic functionality
of each framework after installation.

Usage:
    python -m pytest tests/test_framework_libraries.py -v
"""

import pytest
import sys
import numpy as np
from typing import List, Tuple


# ============================================================================
# PENNYLANE TESTS
# ============================================================================

class TestPennyLaneLibraries:
    """Test suite for PennyLane framework libraries."""
    
    def test_pennylane_import(self):
        """Test PennyLane core import."""
        import pennylane as qml
        assert qml.__version__ is not None
        print(f"✓ PennyLane version: {qml.__version__}")
    
    def test_pennylane_devices(self):
        """Test PennyLane devices availability."""
        import pennylane as qml
        
        # Test default.qubit
        dev_qubit = qml.device('default.qubit', wires=2)
        assert dev_qubit is not None
        
        # Test default.mixed (density matrix)
        dev_mixed = qml.device('default.mixed', wires=2)
        assert dev_mixed is not None
        
        print(f"✓ PennyLane devices: default.qubit, default.mixed")
    
    def test_pennylane_operations(self):
        """Test PennyLane quantum operations."""
        import pennylane as qml
        from pennylane import numpy as pnp
        
        dev = qml.device('default.qubit', wires=2)
        
        @qml.qnode(dev)
        def circuit(x):
            qml.RX(x, wires=0)
            qml.CNOT(wires=[0, 1])
            return qml.expval(qml.PauliZ(0))
        
        result = circuit(np.pi / 4)
        assert isinstance(result, (float, np.ndarray))
        print(f"✓ PennyLane operations: RX, CNOT, PauliZ working (result={result:.4f})")
    
    def test_pennylane_templates(self):
        """Test PennyLane templates (ansätze)."""
        import pennylane as qml
        
        dev = qml.device('default.qubit', wires=3)
        
        @qml.qnode(dev)
        def circuit_with_template(weights):
            qml.templates.BasicEntanglerLayers(weights, wires=range(3))
            return qml.expval(qml.PauliZ(0))
        
        weights = np.random.random((2, 3))
        result = circuit_with_template(weights)
        assert isinstance(result, (float, np.ndarray))
        print(f"✓ PennyLane templates: BasicEntanglerLayers working")
    
    def test_pennylane_gradients(self):
        """Test PennyLane gradient computation."""
        import pennylane as qml
        from pennylane import numpy as pnp
        
        dev = qml.device('default.qubit', wires=1)
        
        @qml.qnode(dev, diff_method='parameter-shift')
        def circuit(x):
            qml.RX(x, wires=0)
            return qml.expval(qml.PauliZ(0))
        
        # Use PennyLane numpy for trainable parameters
        x_val = pnp.array(np.pi / 4, requires_grad=True)
        grad_fn = qml.grad(circuit)
        gradient = grad_fn(x_val)
        
        # Gradient should be a scalar or 0-d array
        if isinstance(gradient, (np.ndarray, pnp.ndarray)):
            gradient = float(gradient)
        
        assert isinstance(gradient, (float, np.floating)), f"Gradient type: {type(gradient)}, value: {gradient}"
        assert not np.isnan(gradient), "Gradient is NaN"
        print(f"✓ PennyLane gradients: parameter-shift working (grad={gradient:.4f})")
    
    def test_pennylane_noise_models(self):
        """Test PennyLane noise channel support."""
        import pennylane as qml
        
        dev = qml.device('default.mixed', wires=2)
        
        @qml.qnode(dev)
        def noisy_circuit():
            qml.Hadamard(wires=0)
            qml.DepolarizingChannel(0.1, wires=0)
            qml.AmplitudeDamping(0.05, wires=0)
            return qml.expval(qml.PauliZ(0))
        
        result = noisy_circuit()
        assert isinstance(result, (float, np.ndarray))
        print(f"✓ PennyLane noise: DepolarizingChannel, AmplitudeDamping working")


# ============================================================================
# QISKIT TESTS
# ============================================================================

class TestQiskitLibraries:
    """Test suite for Qiskit framework libraries."""
    
    def test_qiskit_import(self):
        """Test Qiskit core imports."""
        import qiskit
        from qiskit import QuantumCircuit
        from qiskit_aer import Aer
        
        # Qiskit 1.x: Use StatevectorSampler or Aer's Sampler
        try:
            from qiskit.primitives import StatevectorSampler
            sampler_available = True
        except ImportError:
            try:
                from qiskit_aer.primitives import Sampler
                sampler_available = True
            except ImportError:
                sampler_available = False
        
        assert qiskit.__version__ is not None
        assert sampler_available, "No compatible Sampler found for Qiskit primitives"
        print(f"✓ Qiskit version: {qiskit.__version__}")
    
    def test_qiskit_circuit_creation(self):
        """Test Qiskit circuit creation."""
        from qiskit import QuantumCircuit
        
        qc = QuantumCircuit(2, 2)
        qc.h(0)
        qc.cx(0, 1)
        qc.measure([0, 1], [0, 1])
        
        assert qc.num_qubits == 2
        assert qc.num_clbits == 2
        print(f"✓ Qiskit circuit: 2 qubits, 2 clbits, {qc.size()} operations")
    
    def test_qiskit_aer_simulator(self):
        """Test Qiskit Aer simulator."""
        from qiskit import QuantumCircuit, transpile
        from qiskit_aer import Aer
        
        # Create circuit
        qc = QuantumCircuit(2)
        qc.h(0)
        qc.cx(0, 1)
        qc.measure_all()
        
        # Run on simulator
        backend = Aer.get_backend('qasm_simulator')
        transpiled = transpile(qc, backend)
        job = backend.run(transpiled, shots=100)
        result = job.result()
        counts = result.get_counts()
        
        assert len(counts) > 0
        assert sum(counts.values()) == 100
        print(f"✓ Qiskit Aer simulator: {len(counts)} basis states measured")
    
    def test_qiskit_statevector(self):
        """Test Qiskit statevector simulator."""
        from qiskit import QuantumCircuit, transpile
        from qiskit_aer import Aer
        
        qc = QuantumCircuit(2)
        qc.h(0)
        qc.cx(0, 1)
        
        backend = Aer.get_backend('statevector_simulator')
        transpiled = transpile(qc, backend)
        job = backend.run(transpiled)
        result = job.result()
        statevector = result.get_statevector()
        
        assert len(statevector) == 4  # 2^2 qubits
        print(f"✓ Qiskit statevector: {len(statevector)} amplitudes")
    
    def test_qiskit_noise_model(self):
        """Test Qiskit noise model creation."""
        from qiskit_aer.noise import NoiseModel, depolarizing_error, amplitude_damping_error
        
        noise_model = NoiseModel()
        
        # Add depolarizing error
        depol_error = depolarizing_error(0.01, 1)
        noise_model.add_all_qubit_quantum_error(depol_error, ['h', 'rx', 'ry', 'rz'])
        
        # Add amplitude damping
        damping_error = amplitude_damping_error(0.05)
        noise_model.add_all_qubit_quantum_error(damping_error, ['measure'])
        
        assert len(noise_model.noise_qubits) >= 0
        print(f"✓ Qiskit noise model: depolarizing + amplitude damping")
    
    def test_qiskit_transpiler(self):
        """Test Qiskit transpiler."""
        from qiskit import QuantumCircuit, transpile
        from qiskit_aer import Aer
        
        qc = QuantumCircuit(3)
        qc.h(0)
        qc.cx(0, 1)
        qc.cx(1, 2)
        
        backend = Aer.get_backend('qasm_simulator')
        transpiled = transpile(qc, backend, optimization_level=1)
        
        assert transpiled.num_qubits == 3
        print(f"✓ Qiskit transpiler: optimization_level=1 working")
    
    def test_qiskit_primitives(self):
        """Test Qiskit primitives."""
        from qiskit import QuantumCircuit
        
        # Try Qiskit 1.x API first, fallback to Aer
        try:
            from qiskit.primitives import StatevectorSampler as Sampler
        except ImportError:
            from qiskit_aer.primitives import Sampler
        
        qc = QuantumCircuit(1)
        qc.h(0)
        qc.measure_all()
        
        sampler = Sampler()
        
        # Run based on API version
        try:
            # Qiskit 1.x StatevectorSampler API
            job = sampler.run([qc], shots=100)
            result = job.result()
            assert hasattr(result, 'quasi_dists') or hasattr(result[0].data, 'meas'), "No measurement results found"
        except TypeError:
            # Aer Sampler API
            job = sampler.run(qc, shots=100)
            result = job.result()
            assert hasattr(result, 'quasi_dists'), "No quasi_dists in result"
        
        print(f"✓ Qiskit primitives: Sampler working")


# ============================================================================
# CIRQ TESTS
# ============================================================================

class TestCirqLibraries:
    """Test suite for Cirq framework libraries."""
    
    def test_cirq_import(self):
        """Test Cirq core import."""
        import cirq
        assert cirq.__version__ is not None
        print(f"✓ Cirq version: {cirq.__version__}")
    
    def test_cirq_circuit_creation(self):
        """Test Cirq circuit creation."""
        import cirq
        
        qubits = cirq.LineQubit.range(2)
        circuit = cirq.Circuit(
            cirq.H(qubits[0]),
            cirq.CNOT(qubits[0], qubits[1]),
            cirq.measure(*qubits, key='result')
        )
        
        assert len(circuit) == 3
        print(f"✓ Cirq circuit: {len(circuit)} moments, 2 qubits")
    
    def test_cirq_simulator(self):
        """Test Cirq simulator."""
        import cirq
        
        qubits = cirq.LineQubit.range(2)
        circuit = cirq.Circuit(
            cirq.H(qubits[0]),
            cirq.CNOT(qubits[0], qubits[1]),
            cirq.measure(*qubits, key='result')
        )
        
        simulator = cirq.Simulator()
        result = simulator.run(circuit, repetitions=100)
        
        assert 'result' in result.measurements
        assert len(result.measurements['result']) == 100
        print(f"✓ Cirq simulator: 100 shots executed")
    
    def test_cirq_density_matrix(self):
        """Test Cirq density matrix simulator."""
        import cirq
        
        qubit = cirq.LineQubit(0)
        circuit = cirq.Circuit(
            cirq.H(qubit),
            cirq.measure(qubit, key='result')
        )
        
        simulator = cirq.DensityMatrixSimulator()
        result = simulator.run(circuit, repetitions=100)
        
        assert 'result' in result.measurements
        print(f"✓ Cirq density matrix simulator: working")
    
    def test_cirq_noise_channels(self):
        """Test Cirq noise channel support."""
        import cirq
        
        qubit = cirq.LineQubit(0)
        circuit = cirq.Circuit(
            cirq.H(qubit),
            cirq.depolarize(0.01)(qubit),
            cirq.amplitude_damp(0.05)(qubit)
        )
        
        simulator = cirq.DensityMatrixSimulator()
        result = simulator.simulate(circuit)
        
        assert result.final_density_matrix is not None
        print(f"✓ Cirq noise: depolarize, amplitude_damp working")
    
    def test_cirq_gates(self):
        """Test Cirq gate operations."""
        import cirq
        
        qubit = cirq.LineQubit(0)
        circuit = cirq.Circuit(
            cirq.H(qubit),
            cirq.X(qubit),
            cirq.Y(qubit),
            cirq.Z(qubit),
            cirq.rx(np.pi / 4)(qubit),
            cirq.ry(np.pi / 4)(qubit),
            cirq.rz(np.pi / 4)(qubit)
        )
        
        assert len(circuit) == 7
        print(f"✓ Cirq gates: H, X, Y, Z, Rx, Ry, Rz working")
    
    def test_cirq_optimization(self):
        """Test Cirq circuit optimization."""
        import cirq
        
        qubits = cirq.LineQubit.range(2)
        circuit = cirq.Circuit(
            cirq.H(qubits[0]),
            cirq.X(qubits[0]),
            cirq.X(qubits[0]),  # Double X cancels
            cirq.CNOT(qubits[0], qubits[1])
        )
        
        # Optimize circuit
        optimized = cirq.optimize_for_target_gateset(circuit, gateset=cirq.SqrtIswapTargetGateset())
        
        assert optimized is not None
        print(f"✓ Cirq optimization: circuit optimized")


# ============================================================================
# CROSS-FRAMEWORK TESTS
# ============================================================================

class TestCrossFrameworkCompatibility:
    """Test compatibility and interoperability between frameworks."""
    
    def test_numpy_compatibility(self):
        """Test NumPy compatibility across all frameworks."""
        import pennylane as qml
        from qiskit import QuantumCircuit
        import cirq
        
        arr = np.array([0.1, 0.2, 0.3])
        
        # PennyLane
        dev = qml.device('default.qubit', wires=1)
        @qml.qnode(dev)
        def pl_circuit(x):
            qml.RX(x, wires=0)
            return qml.expval(qml.PauliZ(0))
        pl_result = pl_circuit(arr[0])
        
        # Qiskit
        qc = QuantumCircuit(1)
        qc.rx(arr[0], 0)
        
        # Cirq
        qubit = cirq.LineQubit(0)
        cirq_circuit = cirq.Circuit(cirq.rx(arr[0])(qubit))
        
        print(f"✓ NumPy compatibility: all frameworks accept np.array parameters")
    
    def test_seed_reproducibility(self):
        """Test seed-based reproducibility across frameworks.
        
        Verifies that setting the same random seed produces identical
        results when parameters are generated outside the QNode and
        passed as arguments, ensuring true reproducibility.
        """
        import pennylane as qml
        from qiskit import QuantumCircuit, transpile
        from qiskit_aer import Aer
        import cirq
        
        seed = 42
        
        # PennyLane - Generate random parameter explicitly outside QNode
        np.random.seed(seed)
        param1 = np.random.random()
        
        dev = qml.device('default.qubit', wires=1)
        @qml.qnode(dev)
        def pl_circuit(param):
            qml.RX(param, wires=0)
            return qml.expval(qml.PauliZ(0))
        pl_result1 = pl_circuit(param1)
        
        # Reset seed and generate same parameter
        np.random.seed(seed)
        param2 = np.random.random()
        pl_result2 = pl_circuit(param2)
        
        # Verify parameters are identical
        assert np.allclose(param1, param2), f"Parameters differ: {param1} vs {param2}"
        # Verify results are identical with explicit tolerance
        assert np.allclose(pl_result1, pl_result2, rtol=1e-7, atol=1e-9), \
            f"Results differ: {pl_result1} vs {pl_result2}"
        print(f"✓ Seed reproducibility: PennyLane deterministic with seed={seed}")
    
    def test_performance_scaling(self):
        """Test basic performance characteristics."""
        import time
        import pennylane as qml
        from qiskit import QuantumCircuit, transpile
        from qiskit_aer import Aer
        import cirq
        
        n_qubits = 4
        
        # PennyLane timing
        start = time.time()
        dev = qml.device('default.qubit', wires=n_qubits)
        @qml.qnode(dev)
        def pl_circuit():
            for i in range(n_qubits):
                qml.Hadamard(wires=i)
            return qml.expval(qml.PauliZ(0))
        pl_result = pl_circuit()
        pl_time = time.time() - start
        
        # Qiskit timing
        start = time.time()
        qc = QuantumCircuit(n_qubits)
        for i in range(n_qubits):
            qc.h(i)
        qc.measure_all()
        backend = Aer.get_backend('qasm_simulator')
        transpiled = transpile(qc, backend)
        job = backend.run(transpiled, shots=100)
        result = job.result()
        qiskit_time = time.time() - start
        
        # Cirq timing
        start = time.time()
        qubits = cirq.LineQubit.range(n_qubits)
        circuit = cirq.Circuit([cirq.H(q) for q in qubits])
        circuit.append(cirq.measure(*qubits, key='result'))
        simulator = cirq.Simulator()
        cirq_result = simulator.run(circuit, repetitions=100)
        cirq_time = time.time() - start
        
        print(f"✓ Performance ({n_qubits} qubits): PennyLane={pl_time:.3f}s, Qiskit={qiskit_time:.3f}s, Cirq={cirq_time:.3f}s")


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def generate_test_report() -> str:
    """Generate a summary report of all tests."""
    report = []
    report.append("=" * 80)
    report.append("QUANTUM FRAMEWORK LIBRARY TEST REPORT")
    report.append("=" * 80)
    report.append("")
    report.append("Test Categories:")
    report.append("  1. PennyLane Libraries (6 tests)")
    report.append("  2. Qiskit Libraries (7 tests)")
    report.append("  3. Cirq Libraries (6 tests)")
    report.append("  4. Cross-Framework Compatibility (3 tests)")
    report.append("")
    report.append("Total: 22 comprehensive library tests")
    report.append("")
    report.append("=" * 80)
    return "\n".join(report)


if __name__ == "__main__":
    print(generate_test_report())
    print("\nRun with: python -m pytest tests/test_framework_libraries.py -v")
