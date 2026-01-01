#!/usr/bin/env python3
"""
Smoke Tests for Multiframework Execution

Purpose: Quick validation that all frameworks are importable and can execute
         a minimal configuration without errors.

Usage: python tests/smoke_multiframework.py
"""

import sys
import os
import yaml
import json
import numpy as np
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

def test_imports():
    """Test that all frameworks can be imported."""
    print("=" * 60)
    print("TEST 1: Framework Imports")
    print("=" * 60)
    
    results = {}
    
    # Test PennyLane
    try:
        import pennylane as qml
        print(f"✅ PennyLane imported successfully (v{qml.__version__})")
        results['pennylane'] = True
    except ImportError as e:
        print(f"❌ PennyLane import failed: {e}")
        results['pennylane'] = False
    
    # Test Qiskit
    try:
        import qiskit
        print(f"✅ Qiskit imported successfully (v{qiskit.__version__})")
        results['qiskit'] = True
    except ImportError as e:
        print(f"❌ Qiskit import failed: {e}")
        results['qiskit'] = False
    
    # Test Cirq
    try:
        import cirq
        print(f"✅ Cirq imported successfully (v{cirq.__version__})")
        results['cirq'] = True
    except ImportError as e:
        print(f"❌ Cirq import failed: {e}")
        results['cirq'] = False
    
    print()
    return all(results.values()), results


def test_config_loading():
    """Test that unified config can be loaded."""
    print("=" * 60)
    print("TEST 2: Configuration Loading")
    print("=" * 60)
    
    config_path = Path("configs/experiment_unified.yaml")
    
    if not config_path.exists():
        print(f"❌ Config file not found: {config_path}")
        return False, None
    
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        print(f"✅ Config loaded successfully")
        print(f"   - Version: {config['metadata']['version']}")
        print(f"   - Run ID: {config['metadata']['run_id']}")
        print(f"   - Global seed: {config['seeds']['global']}")
        return True, config
    except Exception as e:
        print(f"❌ Config loading failed: {e}")
        return False, None


def test_pennylane_minimal():
    """Run minimal PennyLane circuit."""
    print("=" * 60)
    print("TEST 3: PennyLane Minimal Execution")
    print("=" * 60)
    
    try:
        import pennylane as qml
        
        # Set seed
        np.random.seed(42)
        
        # Create device
        dev = qml.device('default.qubit', wires=2, shots=100)
        
        # Define circuit
        @qml.qnode(dev)
        def circuit(params):
            qml.RY(params[0], wires=0)
            qml.RY(params[1], wires=1)
            qml.CNOT(wires=[0, 1])
            return qml.expval(qml.PauliZ(0))
        
        # Execute
        params = np.array([0.1, 0.2])
        result = circuit(params)
        
        print(f"✅ PennyLane circuit executed")
        print(f"   - Result: {result:.4f}")
        print(f"   - Type: {type(result)}")
        
        return True, result
    except Exception as e:
        print(f"❌ PennyLane execution failed: {e}")
        import traceback
        traceback.print_exc()
        return False, None


def test_qiskit_minimal():
    """Run minimal Qiskit circuit."""
    print("=" * 60)
    print("TEST 4: Qiskit Minimal Execution")
    print("=" * 60)
    
    try:
        from qiskit import QuantumCircuit
        from qiskit_aer import Aer
        import numpy as np
        
        # Set seed
        np.random.seed(42)
        
        # Create circuit
        qc = QuantumCircuit(2, 1)
        qc.ry(0.1, 0)
        qc.ry(0.2, 1)
        qc.cx(0, 1)
        qc.measure(0, 0)
        
        # Execute
        backend = Aer.get_backend('qasm_simulator')
        job = backend.run(qc, shots=100, seed_simulator=42)
        result = job.result()
        counts = result.get_counts()
        
        print(f"✅ Qiskit circuit executed")
        print(f"   - Counts: {counts}")
        print(f"   - Total shots: {sum(counts.values())}")
        
        return True, counts
    except Exception as e:
        print(f"❌ Qiskit execution failed: {e}")
        import traceback
        traceback.print_exc()
        return False, None


def test_cirq_minimal():
    """Run minimal Cirq circuit."""
    print("=" * 60)
    print("TEST 5: Cirq Minimal Execution")
    print("=" * 60)
    
    try:
        import cirq
        import numpy as np
        
        # Set seed
        np.random.seed(42)
        
        # Create qubits
        q0, q1 = cirq.GridQubit.rect(1, 2)
        
        # Create circuit
        circuit = cirq.Circuit(
            cirq.ry(0.1)(q0),
            cirq.ry(0.2)(q1),
            cirq.CNOT(q0, q1),
            cirq.measure(q0, key='result')
        )
        
        # Execute
        simulator = cirq.Simulator(seed=42)
        result = simulator.run(circuit, repetitions=100)
        counts = result.histogram(key='result')
        
        print(f"✅ Cirq circuit executed")
        print(f"   - Histogram: {counts}")
        print(f"   - Total samples: {sum(counts.values())}")
        
        return True, counts
    except Exception as e:
        print(f"❌ Cirq execution failed: {e}")
        import traceback
        traceback.print_exc()
        return False, None


def test_output_schema():
    """Test that outputs can be written in expected schema."""
    print("=" * 60)
    print("TEST 6: Output Schema Validation")
    print("=" * 60)
    
    try:
        # Create test output directory
        test_dir = Path("results/test_smoke")
        test_dir.mkdir(parents=True, exist_ok=True)
        
        # Test metrics.csv schema
        metrics_path = test_dir / "metrics.csv"
        with open(metrics_path, 'w') as f:
            f.write("framework,run_id,dataset,seed,config_id,accuracy,loss,time\n")
            f.write("pennylane,test_001,iris,42,cfg_001,0.95,0.12,10.5\n")
        
        print(f"✅ Metrics CSV written: {metrics_path}")
        
        # Test manifest.json schema
        manifest_path = test_dir / "manifest.json"
        manifest = {
            "manifest_version": "1.0",
            "run_id": "test_001",
            "framework": "pennylane",
            "timestamp": "2025-12-27T10:00:00Z",
            "seeds": {"global": 42},
            "status": "success"
        }
        with open(manifest_path, 'w') as f:
            json.dump(manifest, f, indent=2)
        
        print(f"✅ Manifest JSON written: {manifest_path}")
        
        # Cleanup
        metrics_path.unlink()
        manifest_path.unlink()
        test_dir.rmdir()
        
        return True
    except Exception as e:
        print(f"❌ Output schema test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_determinism():
    """Test that same seed produces same results."""
    print("=" * 60)
    print("TEST 7: Determinism (Reproducibility)")
    print("=" * 60)
    
    try:
        import pennylane as qml
        
        def run_circuit(seed):
            np.random.seed(seed)
            dev = qml.device('default.qubit', wires=1)
            
            @qml.qnode(dev)
            def circuit():
                qml.RY(np.random.rand(), wires=0)
                return qml.expval(qml.PauliZ(0))
            
            return circuit()
        
        # Run twice with same seed
        result1 = run_circuit(42)
        result2 = run_circuit(42)
        
        if np.allclose(result1, result2):
            print(f"✅ Determinism verified")
            print(f"   - Run 1: {result1:.6f}")
            print(f"   - Run 2: {result2:.6f}")
            print(f"   - Diff: {abs(result1 - result2):.2e}")
            return True
        else:
            print(f"❌ Non-deterministic results!")
            print(f"   - Run 1: {result1:.6f}")
            print(f"   - Run 2: {result2:.6f}")
            print(f"   - Diff: {abs(result1 - result2):.6f}")
            return False
            
    except Exception as e:
        print(f"❌ Determinism test failed: {e}")
        return False


def main():
    """Run all smoke tests."""
    print("\n" + "=" * 60)
    print("MULTIFRAMEWORK SMOKE TESTS")
    print("=" * 60)
    print()
    
    results = {}
    
    # Test 1: Imports
    passed, data = test_imports()
    results['imports'] = passed
    print()
    
    # Test 2: Config loading
    passed, config = test_config_loading()
    results['config'] = passed
    print()
    
    # Test 3: PennyLane
    passed, _ = test_pennylane_minimal()
    results['pennylane'] = passed
    print()
    
    # Test 4: Qiskit
    passed, _ = test_qiskit_minimal()
    results['qiskit'] = passed
    print()
    
    # Test 5: Cirq
    passed, _ = test_cirq_minimal()
    results['cirq'] = passed
    print()
    
    # Test 6: Output schema
    passed = test_output_schema()
    results['output_schema'] = passed
    print()
    
    # Test 7: Determinism
    passed = test_determinism()
    results['determinism'] = passed
    print()
    
    # Summary
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    
    total = len(results)
    passed = sum(results.values())
    
    for test_name, test_passed in results.items():
        status = "✅ PASS" if test_passed else "❌ FAIL"
        print(f"{status} - {test_name}")
    
    print()
    print(f"Total: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 All smoke tests passed! System is ready for execution.")
        return 0
    else:
        print(f"\n⚠️  {total - passed} test(s) failed. Please fix before proceeding.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
