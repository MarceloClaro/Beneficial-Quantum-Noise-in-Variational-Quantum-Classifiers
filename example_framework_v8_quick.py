#!/usr/bin/env python3
"""
Quick Example: Advanced Quantum Framework V8.0

Demonstrates basic usage of the framework with synthetic data.
For real usage, install DeepChem: bash install_deepchem.sh

Usage:
    python example_framework_v8_quick.py
"""

import numpy as np
from framework_quantum_advanced_v8 import (
    AdvancedConfig,
    AdvancedVQC,
    DeepChemDatasetLoader,
    QuantumComplexityAnalyzer,
    ZeroNoiseExtrapolation
)

def main():
    print("="*70)
    print("ADVANCED QUANTUM FRAMEWORK V8.0 - QUICK EXAMPLE")
    print("="*70)
    
    # 1. Load Dataset
    print("\n1. Loading Dataset...")
    loader = DeepChemDatasetLoader(verbose=True)
    dataset = loader.load_dataset("BACE", max_samples=100)
    
    print(f"   Train samples: {len(dataset['X_train'])}")
    print(f"   Test samples: {len(dataset['X_test'])}")
    print(f"   Features: {dataset['n_features']}")
    
    # 2. Configure Framework
    print("\n2. Configuring Framework...")
    config = AdvancedConfig(
        framework="pennylane",
        n_qubits=4,
        n_layers=2,
        n_epochs=20,
        learning_rate=0.01,
        noise_level=0.01,
        error_mitigation="zne",
        results_dir="resultados_example_v8",
        verbose=True
    )
    
    print(f"   Framework: {config.framework}")
    print(f"   Qubits: {config.n_qubits}, Layers: {config.n_layers}")
    print(f"   Error Mitigation: {config.error_mitigation}")
    
    # 3. Analyze Quantum Complexity
    print("\n3. Analyzing Quantum Complexity...")
    analyzer = QuantumComplexityAnalyzer(verbose=True)
    metrics = analyzer.analyze(
        n_qubits=config.n_qubits,
        n_layers=config.n_layers,
        entanglement="linear"
    )
    
    # 4. Train VQC
    print("\n4. Training VQC...")
    vqc = AdvancedVQC(config)
    vqc.fit(dataset['X_train'], dataset['y_train'])
    
    # 5. Evaluate
    print("\n5. Evaluating...")
    train_acc = vqc.score(dataset['X_train'], dataset['y_train'])
    test_acc = vqc.score(dataset['X_test'], dataset['y_test'])
    
    print(f"   Training Accuracy: {train_acc:.4f}")
    print(f"   Test Accuracy: {test_acc:.4f}")
    print(f"   Generalization Gap: {(train_acc - test_acc):.4f}")
    
    # 6. Demonstrate ZNE
    print("\n6. Demonstrating Zero-Noise Extrapolation...")
    zne = ZeroNoiseExtrapolation(
        scale_factors=[1.0, 1.5, 2.0, 2.5],
        method="linear"
    )
    
    # Simulated expectation values at different noise scales
    expectation_values = [0.85, 0.80, 0.75, 0.70]
    mitigated_value = zne.extrapolate(expectation_values)
    
    print(f"   Noise scales: {zne.scale_factors}")
    print(f"   Expectation values: {expectation_values}")
    print(f"   Extrapolated (zero-noise): {mitigated_value:.4f}")
    
    # 7. Summary
    print("\n" + "="*70)
    print("EXAMPLE COMPLETED SUCCESSFULLY ✓")
    print("="*70)
    print("\nKey Results:")
    print(f"  • Dataset: {dataset['name']}")
    print(f"  • Test Accuracy: {test_acc:.4f}")
    print(f"  • Circuit Depth: {metrics['circuit_depth']}")
    print(f"  • Total Gates: {metrics['total_gates']}")
    print(f"  • Barren Plateau Risk: {metrics['barren_plateau_risk']}")
    print(f"\nResults saved to: {config.results_dir}/")
    
    print("\nNext Steps:")
    print("  1. Install DeepChem for real molecular datasets:")
    print("     bash install_deepchem.sh")
    print("  2. Run full framework:")
    print("     python -m framework_quantum_advanced_v8")
    print("  3. See documentation:")
    print("     cat README_FRAMEWORK_ADVANCED_V8.md")


if __name__ == "__main__":
    main()
