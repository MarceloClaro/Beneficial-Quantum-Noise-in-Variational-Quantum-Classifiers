# Quantum Framework Library Testing Guide

## Overview

Comprehensive testing infrastructure for validating quantum framework installations and functionality across PennyLane, Qiskit, and Cirq.

## Test Files

### 1. `tests/test_framework_libraries.py`

Comprehensive test suite with 22 tests across 4 categories:

#### Test Categories

**PennyLane Tests (6 tests):**
- `test_pennylane_import` - Core library import
- `test_pennylane_devices` - Device availability (default.qubit, default.mixed)
- `test_pennylane_operations` - Quantum operations (RX, CNOT, PauliZ)
- `test_pennylane_templates` - Ansätze templates (BasicEntanglerLayers)
- `test_pennylane_gradients` - Gradient computation (parameter-shift)
- `test_pennylane_noise_models` - Noise channels (depolarizing, amplitude damping)

**Qiskit Tests (7 tests):**
- `test_qiskit_import` - Core library import (Qiskit 2.x API)
- `test_qiskit_circuit_creation` - Circuit creation
- `test_qiskit_aer_simulator` - Aer simulator (qasm_simulator)
- `test_qiskit_statevector` - Statevector simulator
- `test_qiskit_noise_model` - Noise model creation
- `test_qiskit_transpiler` - Circuit transpilation
- `test_qiskit_primitives` - Primitives API (Sampler)

**Cirq Tests (6 tests):**
- `test_cirq_import` - Core library import
- `test_cirq_circuit_creation` - Circuit creation
- `test_cirq_simulator` - Simulator execution
- `test_cirq_density_matrix` - Density matrix simulator
- `test_cirq_noise_channels` - Noise channels (depolarize, amplitude_damp)
- `test_cirq_gates` - Gate operations (H, X, Y, Z, Rx, Ry, Rz)
- `test_cirq_optimization` - Circuit optimization

**Cross-Framework Tests (3 tests):**
- `test_numpy_compatibility` - NumPy array compatibility
- `test_seed_reproducibility` - Deterministic execution with seeds
- `test_performance_scaling` - Performance characteristics

## GitHub Actions Workflows

### 1. `quantum-frameworks-test.yml`

Automated testing workflow triggered on:
- Push to main/develop/copilot branches
- Pull requests
- Manual dispatch
- Weekly schedule (Sunday midnight)

**Jobs:**

1. **test-quantum-frameworks**
   - Tests all 3 Python versions (3.9, 3.10, 3.11)
   - Runs 22 comprehensive tests
   - Generates test summaries
   - Uploads test artifacts

2. **test-compatibility-matrix**
   - Tests framework interoperability
   - Generates compatibility report
   - Validates cross-framework features

3. **notification**
   - Reports final test status
   - Creates workflow summary

## Usage

### Local Testing

```bash
# Run all library tests
python -m pytest tests/test_framework_libraries.py -v

# Run specific framework tests
python -m pytest tests/test_framework_libraries.py::TestPennyLaneLibraries -v
python -m pytest tests/test_framework_libraries.py::TestQiskitLibraries -v
python -m pytest tests/test_framework_libraries.py::TestCirqLibraries -v

# Run cross-framework tests
python -m pytest tests/test_framework_libraries.py::TestCrossFrameworkCompatibility -v

# With detailed output
python -m pytest tests/test_framework_libraries.py -v --tb=short --durations=10
```

### GitHub Actions

The workflow runs automatically on push/PR. To manually trigger:

1. Go to repository → Actions
2. Select "Quantum Frameworks Library Tests"
3. Click "Run workflow"
4. Select branch and run

### Viewing Results

**GitHub Actions:**
- Check "Actions" tab for workflow runs
- View detailed logs for each job
- Download test artifacts (30-day retention)
- Check workflow summary for quick overview

**Local:**
```bash
# Generate test report
python tests/test_framework_libraries.py

# Run with coverage
python -m pytest tests/test_framework_libraries.py --cov=. --cov-report=html
```

## Expected Results

### Successful Run

```
tests/test_framework_libraries.py::TestPennyLaneLibraries::test_pennylane_import PASSED
tests/test_framework_libraries.py::TestPennyLaneLibraries::test_pennylane_devices PASSED
tests/test_framework_libraries.py::TestPennyLaneLibraries::test_pennylane_operations PASSED
tests/test_framework_libraries.py::TestPennyLaneLibraries::test_pennylane_templates PASSED
tests/test_framework_libraries.py::TestPennyLaneLibraries::test_pennylane_gradients PASSED
tests/test_framework_libraries.py::TestPennyLaneLibraries::test_pennylane_noise_models PASSED
tests/test_framework_libraries.py::TestQiskitLibraries::test_qiskit_import PASSED
tests/test_framework_libraries.py::TestQiskitLibraries::test_qiskit_circuit_creation PASSED
tests/test_framework_libraries.py::TestQiskitLibraries::test_qiskit_aer_simulator PASSED
tests/test_framework_libraries.py::TestQiskitLibraries::test_qiskit_statevector PASSED
tests/test_framework_libraries.py::TestQiskitLibraries::test_qiskit_noise_model PASSED
tests/test_framework_libraries.py::TestQiskitLibraries::test_qiskit_transpiler PASSED
tests/test_framework_libraries.py::TestQiskitLibraries::test_qiskit_primitives PASSED
tests/test_framework_libraries.py::TestCirqLibraries::test_cirq_import PASSED
tests/test_framework_libraries.py::TestCirqLibraries::test_cirq_circuit_creation PASSED
tests/test_framework_libraries.py::TestCirqLibraries::test_cirq_simulator PASSED
tests/test_framework_libraries.py::TestCirqLibraries::test_cirq_density_matrix PASSED
tests/test_framework_libraries.py::TestCirqLibraries::test_cirq_noise_channels PASSED
tests/test_framework_libraries.py::TestCirqLibraries::test_cirq_gates PASSED
tests/test_framework_libraries.py::TestCirqLibraries::test_cirq_optimization PASSED
tests/test_framework_libraries.py::TestCrossFrameworkCompatibility::test_numpy_compatibility PASSED
tests/test_framework_libraries.py::TestCrossFrameworkCompatibility::test_seed_reproducibility PASSED
tests/test_framework_libraries.py::TestCrossFrameworkCompatibility::test_performance_scaling PASSED

============================== 22 passed in 12.34s ===============================
```

### Test Summary Output

```
================================================================================
QUANTUM FRAMEWORK LIBRARY TEST REPORT
================================================================================

Test Categories:
  1. PennyLane Libraries (6 tests)
  2. Qiskit Libraries (7 tests)
  3. Cirq Libraries (6 tests)
  4. Cross-Framework Compatibility (3 tests)

Total: 22 comprehensive library tests

================================================================================
```

## Troubleshooting

### Common Issues

**1. Import Errors**
```bash
# Solution: Reinstall packages
pip install --upgrade -r requirements.txt
```

**2. Qiskit API Compatibility**
```bash
# Ensure Qiskit 2.x API is used
pip install qiskit>=2.0 qiskit-aer>=0.17
```

**3. Test Timeout**
```bash
# Increase timeout for slow systems
python -m pytest tests/test_framework_libraries.py -v --timeout=300
```

**4. Memory Issues**
```bash
# Run tests individually
python -m pytest tests/test_framework_libraries.py::TestPennyLaneLibraries -v
python -m pytest tests/test_framework_libraries.py::TestQiskitLibraries -v
python -m pytest tests/test_framework_libraries.py::TestCirqLibraries -v
```

### Debugging

```bash
# Verbose output with full tracebacks
python -m pytest tests/test_framework_libraries.py -vv --tb=long

# Show print statements
python -m pytest tests/test_framework_libraries.py -v -s

# Run specific test with debugging
python -m pytest tests/test_framework_libraries.py::TestQiskitLibraries::test_qiskit_aer_simulator -vv -s
```

## Integration with Existing Tests

The library tests complement existing test files:

- `tests/smoke_multiframework.py` - Quick infrastructure validation (7 tests)
- `tests/test_framework_libraries.py` - Comprehensive library tests (22 tests)
- `tests/test_*.py` - Domain-specific tests

**Recommended Test Sequence:**
1. Smoke tests (fast validation)
2. Library tests (comprehensive validation)
3. Domain tests (specific functionality)

## Performance Benchmarks

Expected execution times (approximate):

| Test Category | Time (local) | Time (CI) |
|---------------|-------------|-----------|
| PennyLane (6 tests) | 2-4 sec | 3-6 sec |
| Qiskit (7 tests) | 3-6 sec | 5-10 sec |
| Cirq (6 tests) | 2-4 sec | 3-6 sec |
| Cross-Framework (3 tests) | 1-3 sec | 2-5 sec |
| **Total (22 tests)** | **8-17 sec** | **13-27 sec** |

## Continuous Integration

### Matrix Strategy

Tests run on:
- Python 3.9, 3.10, 3.11
- Ubuntu latest
- 30-minute timeout per job

### Caching

Pip packages are cached to speed up workflow:
- Cache key: `${{ runner.os }}-pip-quantum-${{ matrix.python-version }}`
- Restores from previous runs with same dependencies

### Artifacts

Test results retained for 30 days:
- pytest reports
- Cache files
- Test logs

## API Compatibility Notes

### Qiskit 2.x Changes

The tests use Qiskit 2.x API:
```python
# Old API (pre-2.0)
from qiskit import Aer, execute

# New API (2.x)
from qiskit_aer import Aer
backend.run(circuit)  # instead of execute(circuit, backend)
```

### PennyLane Devices

Supported devices:
- `default.qubit` - Fast statevector simulator
- `default.mixed` - Density matrix simulator (for noise)
- `lightning.qubit` - Optimized C++ simulator

### Cirq Simulators

Available simulators:
- `cirq.Simulator()` - Statevector simulator
- `cirq.DensityMatrixSimulator()` - Density matrix simulator
- `cirq.Simulator(noise=noise_model)` - Noisy simulation

## Contributing

To add new tests:

1. Add test method to appropriate class
2. Follow naming convention: `test_<feature>_<aspect>`
3. Include assertions and print statements
4. Update documentation
5. Run locally before committing
6. Verify CI passes

## References

- **PennyLane Documentation:** https://docs.pennylane.ai/
- **Qiskit Documentation:** https://docs.quantum.ibm.com/
- **Cirq Documentation:** https://quantumai.google/cirq
- **Pytest Documentation:** https://docs.pytest.org/

## Support

For issues with tests:
1. Check GitHub Actions logs
2. Run tests locally with verbose output
3. Review troubleshooting section
4. Check framework-specific documentation
5. Open issue with test output and error messages
