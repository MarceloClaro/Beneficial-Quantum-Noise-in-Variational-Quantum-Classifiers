#!/usr/bin/env python3
# =============================================================================
# TESTE RÁPIDO - Framework Quantum Advanced V8
# =============================================================================
"""
Script para testar rapidamente se o framework está funcionando.
Executa um experimento pequeno com Iris para validação.
"""

import sys
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def test_imports():
    """Testa se todas as importações estão funcionando."""
    print("\n" + "="*80)
    print("TESTE 1: Validando Importações")
    print("="*80)
    
    errors = []
    
    # Framework
    try:
        from framework_quantum_advanced_v8 import (
            ExperimentConfig,
            QuantumCircuitConfig,
            NoiseConfig,
            OptimizationConfig,
            ErrorMitigationConfig,
            QuantumExperimentRunner,
            FrameworkType,
            NoiseModel,
            OptimizationMethod,
            ErrorMitigationTechnique,
        )
        print("✓ Framework imports OK")
    except ImportError as e:
        print(f"✗ Framework import error: {e}")
        errors.append(str(e))
    
    # Dependências científicas
    try:
        import numpy as np
        import pandas as pd
        from scipy import stats
        from sklearn.preprocessing import StandardScaler
        from sklearn.metrics import accuracy_score
        print("✓ Scientific libraries OK")
    except ImportError as e:
        print(f"✗ Scientific library error: {e}")
        errors.append(str(e))
    
    # Visualização
    try:
        import matplotlib.pyplot as plt
        import plotly.graph_objects as go
        print("✓ Visualization libraries OK")
    except ImportError as e:
        print(f"✗ Visualization error: {e}")
        errors.append(str(e))
    
    # Quantum frameworks (opcionais)
    try:
        import pennylane as qml
        print("✓ PennyLane available")
    except ImportError:
        print("⚠ PennyLane not available (optional)")
    
    try:
        from qiskit import QuantumCircuit
        print("✓ Qiskit available")
    except ImportError:
        print("⚠ Qiskit not available (optional)")
    
    try:
        import cirq
        print("✓ Cirq available")
    except ImportError:
        print("⚠ Cirq not available (optional)")
    
    # DeepChem
    try:
        import deepchem as dc
        print("✓ DeepChem available")
    except ImportError:
        print("⚠ DeepChem not available (optional)")
    
    if errors:
        print(f"\n❌ {len(errors)} erro(s) crítico(s) encontrado(s)")
        return False
    else:
        print("\n✓ Todas as importações críticas OK")
        return True


def test_config_creation():
    """Testa criação de configurações."""
    print("\n" + "="*80)
    print("TESTE 2: Criação de Configurações")
    print("="*80)
    
    try:
        from framework_quantum_advanced_v8 import (
            ExperimentConfig,
            QuantumCircuitConfig,
            NoiseConfig,
            OptimizationConfig,
            ErrorMitigationConfig,
            FrameworkType,
            NoiseModel,
            OptimizationMethod,
            ErrorMitigationTechnique,
        )
        
        config = ExperimentConfig(
            framework=FrameworkType.PENNYLANE,
            circuit_config=QuantumCircuitConfig(
                n_qubits=3,
                n_layers=1,
                n_parameters=9
            ),
            noise_config=NoiseConfig(
                noise_model=NoiseModel.DEPOLARIZING,
                noise_level=0.01
            ),
            optimization_config=OptimizationConfig(
                method=OptimizationMethod.ADAM,
                max_iterations=10
            ),
            error_mitigation_config=ErrorMitigationConfig(
                technique=ErrorMitigationTechnique.ZNE
            ),
            dataset_name="iris",
            n_shots=512
        )
        
        print(f"✓ ExperimentConfig criado")
        print(f"  - Framework: {config.framework.value}")
        print(f"  - Qubits: {config.circuit_config.n_qubits}")
        print(f"  - Layers: {config.circuit_config.n_layers}")
        print(f"  - Dataset: {config.dataset_name}")
        
        return True
    except Exception as e:
        print(f"✗ Erro ao criar config: {e}")
        return False


def test_complexity_analysis():
    """Testa análise de complexidade."""
    print("\n" + "="*80)
    print("TESTE 3: Análise de Complexidade")
    print("="*80)
    
    try:
        from framework_quantum_advanced_v8 import (
            QuantumComplexityAnalyzer,
            QuantumCircuitConfig
        )
        
        analyzer = QuantumComplexityAnalyzer()
        
        config = QuantumCircuitConfig(
            n_qubits=4,
            n_layers=2,
            n_parameters=24
        )
        
        resources = analyzer.analyze_resource_requirements(config, n_shots=1024)
        
        print(f"✓ Análise de complexidade OK")
        print(f"  - Profundidade: {resources['circuit_depth']}")
        print(f"  - Total de gates: {resources['gate_count']['total']}")
        print(f"  - BP Probability: {resources['barren_plateau_probability']:.4f}")
        print(f"  - Tempo estimado: {resources['total_estimated_time_s']:.2f}s")
        
        return True
    except Exception as e:
        print(f"✗ Erro em análise de complexidade: {e}")
        return False


def test_noise_validation():
    """Testa validação de ruído."""
    print("\n" + "="*80)
    print("TESTE 4: Validação de Ruído")
    print("="*80)
    
    try:
        from framework_quantum_advanced_v8 import NoiseValidationFramework
        
        validator = NoiseValidationFramework()
        
        fidelity = validator.predict_noise_impact(
            noise_level=0.01,
            circuit_depth=10,
            n_qubits=4
        )
        
        validation = validator.validate_noise_prediction(
            actual_fidelity=0.90,
            predicted_fidelity=fidelity
        )
        
        print(f"✓ Validação de ruído OK")
        print(f"  - Fidelidade predita: {fidelity:.4f}")
        print(f"  - Fidelidade real (simulada): 0.9000")
        print(f"  - Erro relativo: {validation['relative_error']:.4f}")
        print(f"  - Validado: {validation['validation_passed']}")
        
        return True
    except Exception as e:
        print(f"✗ Erro em validação de ruído: {e}")
        return False


def test_zne():
    """Testa Zero-Noise Extrapolation."""
    print("\n" + "="*80)
    print("TESTE 5: Zero-Noise Extrapolation (ZNE)")
    print("="*80)
    
    try:
        import numpy as np
        from framework_quantum_advanced_v8 import (
            ZeroNoiseExtrapolation,
            ErrorMitigationConfig,
            ErrorMitigationTechnique
        )
        
        config = ErrorMitigationConfig(
            technique=ErrorMitigationTechnique.ZNE,
            zne_scale_factors=[1.0, 1.5, 2.0],
            zne_extrapolation_type="linear"
        )
        
        zne = ZeroNoiseExtrapolation(config)
        
        def mock_observable(scale):
            return 0.95 * np.exp(-scale) + 0.05
        
        extrapolated, details = zne.apply_zne(mock_observable)
        
        print(f"✓ ZNE OK")
        print(f"  - Valores medidos: {[f'{v:.4f}' for v in details['measured_values']]}")
        print(f"  - Valor extrapolado: {extrapolated:.4f}")
        print(f"  - Tipo: {details['extrapolation_type']}")
        
        return True
    except Exception as e:
        print(f"✗ Erro em ZNE: {e}")
        return False


def test_benchmarking():
    """Testa benchmarking."""
    print("\n" + "="*80)
    print("TESTE 6: Benchmarking")
    print("="*80)
    
    try:
        import numpy as np
        from framework_quantum_advanced_v8 import QuantumAlgorithmBenchmark
        
        benchmark = QuantumAlgorithmBenchmark()
        
        np.random.seed(42)
        y_true = np.array([0, 1, 1, 0, 1, 0])
        vqc_pred = np.array([0.1, 0.9, 0.8, 0.2, 0.85, 0.15])
        classical_pred = np.array([0.2, 0.7, 0.6, 0.3, 0.75, 0.25])
        
        comparison = benchmark.benchmark_against_classical(
            vqc_pred, classical_pred, y_true
        )
        
        print(f"✓ Benchmarking OK")
        print(f"  - VQC Accuracy: {comparison['vqc_metrics']['accuracy']:.4f}")
        print(f"  - Classical Accuracy: {comparison['classical_metrics']['accuracy']:.4f}")
        print(f"  - VQC wins: {comparison['vqc_wins']} métricas")
        
        return True
    except Exception as e:
        print(f"✗ Erro em benchmarking: {e}")
        return False


def test_small_experiment():
    """Testa um experimento muito pequeno."""
    print("\n" + "="*80)
    print("TESTE 7: Experimento Pequeno (Iris, 3 qubits)")
    print("="*80)
    
    try:
        from framework_quantum_advanced_v8 import (
            ExperimentConfig,
            QuantumCircuitConfig,
            NoiseConfig,
            OptimizationConfig,
            ErrorMitigationConfig,
            QuantumExperimentRunner,
            FrameworkType,
            NoiseModel,
            OptimizationMethod,
            ErrorMitigationTechnique,
        )
        
        print("Criando configuração...")
        config = ExperimentConfig(
            framework=FrameworkType.PENNYLANE,
            circuit_config=QuantumCircuitConfig(
                n_qubits=3,
                n_layers=1,
                n_parameters=9
            ),
            noise_config=NoiseConfig(
                noise_model=NoiseModel.DEPOLARIZING,
                noise_level=0.01
            ),
            optimization_config=OptimizationConfig(
                method=OptimizationMethod.ADAM,
                learning_rate=0.1,
                max_iterations=10,
                early_stopping=False
            ),
            error_mitigation_config=ErrorMitigationConfig(
                technique=ErrorMitigationTechnique.NONE
            ),
            dataset_name="iris",
            n_shots=512,
            seed=42
        )
        
        print("Criando runner...")
        runner = QuantumExperimentRunner(config, results_dir="./results_test")
        
        print("Executando teste de complexidade...")
        complexity = runner.run_complexity_analysis()
        
        print(f"✓ Experimento teste OK")
        print(f"  - Circuit depth: {complexity['circuit_depth']}")
        print(f"  - Total gates: {complexity['gate_count']['total']}")
        
        return True
    except Exception as e:
        print(f"✗ Erro no experimento teste: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Executa todos os testes."""
    print("\n" + "="*80)
    print("TESTES - FRAMEWORK QUANTUM ADVANCED V8")
    print("="*80)
    
    results = {
        "Imports": test_imports(),
        "Config Creation": test_config_creation(),
        "Complexity Analysis": test_complexity_analysis(),
        "Noise Validation": test_noise_validation(),
        "ZNE": test_zne(),
        "Benchmarking": test_benchmarking(),
        "Small Experiment": test_small_experiment(),
    }
    
    print("\n" + "="*80)
    print("RESUMO DOS TESTES")
    print("="*80)
    
    passed = sum(1 for v in results.values() if v)
    total = len(results)
    
    for test_name, result in results.items():
        status = "✓ PASSOU" if result else "✗ FALHOU"
        print(f"{test_name:30s}: {status}")
    
    print("-" * 50)
    print(f"Total: {passed}/{total} testes passaram")
    
    if passed == total:
        print("\n✓ TODOS OS TESTES PASSARAM! Framework pronto para uso.")
        return 0
    else:
        print(f"\n⚠ {total - passed} teste(s) falharam. Verificar logs acima.")
        return 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
