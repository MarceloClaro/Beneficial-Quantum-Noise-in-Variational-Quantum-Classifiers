#!/usr/bin/env python3
"""
Tests for Advanced Quantum Framework V8.0

Test suite covering:
- Multi-framework support
- ZNE error mitigation
- DeepChem integration
- Complexity analysis
- End-to-end training

Usage:
    python -m pytest tests/test_framework_advanced_v8.py -v
    python tests/test_framework_advanced_v8.py  # Direct execution
"""

import pytest
import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from framework_quantum_advanced_v8 import (
    AdvancedConfig,
    AdvancedVQC,
    ZeroNoiseExtrapolation,
    QuantumComplexityAnalyzer,
    DeepChemDatasetLoader,
    QuantumFramework,
    ErrorMitigationType
)


class TestZeroNoiseExtrapolation:
    """Test Zero-Noise Extrapolation functionality."""
    
    def test_zne_initialization(self):
        """Test ZNE initializes correctly."""
        zne = ZeroNoiseExtrapolation(
            scale_factors=[1.0, 1.5, 2.0],
            method="linear"
        )
        assert zne.scale_factors == [1.0, 1.5, 2.0]
        assert zne.method == "linear"
    
    def test_zne_linear_extrapolation(self):
        """Test linear extrapolation."""
        zne = ZeroNoiseExtrapolation(
            scale_factors=[1.0, 2.0, 3.0],
            method="linear"
        )
        
        # Simulated expectation values with linear noise dependence
        # E(λ) = 0.8 - 0.1*λ
        expectation_values = [0.7, 0.6, 0.5]
        
        result = zne.extrapolate(expectation_values)
        
        # At λ=0, should extrapolate to ~0.8
        assert 0.75 <= result <= 0.85
        assert zne.fit_params is not None
    
    def test_zne_exponential_extrapolation(self):
        """Test exponential extrapolation."""
        zne = ZeroNoiseExtrapolation(
            scale_factors=[1.0, 2.0, 3.0],
            method="exponential"
        )
        
        # Simulated exponential decay
        expectation_values = [0.7, 0.55, 0.45]
        
        result = zne.extrapolate(expectation_values)
        
        # Should give reasonable result
        assert 0.6 <= result <= 1.0
    
    def test_zne_invalid_length(self):
        """Test error handling for mismatched lengths."""
        zne = ZeroNoiseExtrapolation(scale_factors=[1.0, 2.0, 3.0])
        
        with pytest.raises(ValueError):
            zne.extrapolate([0.5, 0.6])  # Wrong length


class TestQuantumComplexityAnalyzer:
    """Test Quantum Complexity Analyzer."""
    
    def test_complexity_analysis(self):
        """Test basic complexity analysis."""
        analyzer = QuantumComplexityAnalyzer(verbose=False)
        
        metrics = analyzer.analyze(
            n_qubits=4,
            n_layers=2,
            entanglement="linear"
        )
        
        assert metrics["n_qubits"] == 4
        assert metrics["n_layers"] == 2
        assert metrics["circuit_depth"] == 2
        assert "total_gates" in metrics
        assert "barren_plateau_risk" in metrics
    
    def test_full_entanglement_complexity(self):
        """Test complexity with full entanglement."""
        analyzer = QuantumComplexityAnalyzer(verbose=False)
        
        metrics = analyzer.analyze(
            n_qubits=4,
            n_layers=2,
            entanglement="full"
        )
        
        # Full entanglement should have more 2-qubit gates
        assert metrics["two_qubit_gates"] > 0
        assert metrics["time_complexity"] == "O(n²·L)"
    
    def test_barren_plateau_risk_assessment(self):
        """Test barren plateau risk assessment."""
        analyzer = QuantumComplexityAnalyzer(verbose=False)
        
        # Small circuit - low risk
        metrics_small = analyzer.analyze(n_qubits=2, n_layers=1)
        assert metrics_small["barren_plateau_risk"] == "LOW"
        
        # Large circuit - higher risk
        metrics_large = analyzer.analyze(n_qubits=10, n_layers=5)
        assert metrics_large["barren_plateau_risk"] in ["MEDIUM", "HIGH"]
    
    def test_export_report(self, tmp_path):
        """Test report export."""
        analyzer = QuantumComplexityAnalyzer(verbose=False)
        analyzer.analyze(n_qubits=4, n_layers=2)
        
        report_path = tmp_path / "complexity_report.md"
        analyzer.export_report(str(report_path))
        
        assert report_path.exists()
        content = report_path.read_text()
        assert "Quantum Circuit Complexity Analysis" in content


class TestDeepChemDatasetLoader:
    """Test DeepChem dataset loader."""
    
    def test_loader_initialization(self):
        """Test loader initializes correctly."""
        loader = DeepChemDatasetLoader(verbose=False)
        assert loader.cache_dir is not None
        assert os.path.exists(loader.cache_dir)
    
    def test_synthetic_dataset_loading(self):
        """Test synthetic dataset generation (when DeepChem not available)."""
        loader = DeepChemDatasetLoader(verbose=False)
        
        dataset = loader.load_dataset("TEST_SYNTHETIC", max_samples=100)
        
        assert "X_train" in dataset
        assert "X_test" in dataset
        assert "y_train" in dataset
        assert "y_test" in dataset
        assert len(dataset["X_train"]) == 100
        assert dataset["X_train"].shape[1] == 16
    
    def test_dataset_shapes(self):
        """Test dataset has correct shapes."""
        loader = DeepChemDatasetLoader(verbose=False)
        dataset = loader.load_dataset("BACE", max_samples=200)
        
        # Training set
        assert dataset["X_train"].shape[0] <= 200
        assert dataset["X_train"].shape[1] > 0
        
        # Test set
        assert dataset["X_test"].shape[0] <= 40  # 200 / 5
        
        # Labels are binary
        assert set(np.unique(dataset["y_train"])).issubset({0, 1})


class TestAdvancedConfig:
    """Test configuration dataclass."""
    
    def test_default_config(self):
        """Test default configuration."""
        config = AdvancedConfig()
        
        assert config.framework == "pennylane"
        assert config.n_qubits == 4
        assert config.n_layers == 2
        assert config.n_epochs == 50
        assert config.error_mitigation == "zne"
    
    def test_custom_config(self):
        """Test custom configuration."""
        config = AdvancedConfig(
            framework="qiskit",
            n_qubits=6,
            n_layers=3,
            error_mitigation="combined"
        )
        
        assert config.framework == "qiskit"
        assert config.n_qubits == 6
        assert config.n_layers == 3
        assert config.error_mitigation == "combined"


class TestAdvancedVQC:
    """Test Advanced VQC classifier."""
    
    def test_vqc_initialization(self):
        """Test VQC initializes correctly."""
        config = AdvancedConfig(n_epochs=10, verbose=False)
        vqc = AdvancedVQC(config)
        
        assert vqc.config == config
        assert vqc.params is None  # Not initialized until fit
        assert hasattr(vqc, 'complexity_analyzer')
    
    def test_vqc_fit_predict(self):
        """Test VQC training and prediction."""
        # Create simple synthetic dataset
        np.random.seed(42)
        X = np.random.randn(50, 16)
        y = (X[:, 0] + X[:, 1] > 0).astype(int)
        
        # Train VQC
        config = AdvancedConfig(n_epochs=5, verbose=False)
        vqc = AdvancedVQC(config)
        vqc.fit(X, y)
        
        # Check parameters initialized
        assert vqc.params is not None
        assert len(vqc.params) > 0
        
        # Check predictions
        y_pred = vqc.predict(X)
        assert len(y_pred) == len(y)
        assert set(np.unique(y_pred)).issubset({0, 1})
    
    def test_vqc_score(self):
        """Test VQC scoring."""
        np.random.seed(42)
        X = np.random.randn(50, 16)
        y = (X[:, 0] + X[:, 1] > 0).astype(int)
        
        config = AdvancedConfig(n_epochs=5, verbose=False)
        vqc = AdvancedVQC(config)
        vqc.fit(X, y)
        
        score = vqc.score(X, y)
        
        assert 0.0 <= score <= 1.0
        # With simple dataset, should get reasonable accuracy
        assert score > 0.4
    
    def test_vqc_predict_proba(self):
        """Test probability predictions."""
        np.random.seed(42)
        X = np.random.randn(20, 16)
        y = (X[:, 0] > 0).astype(int)
        
        config = AdvancedConfig(n_epochs=3, verbose=False)
        vqc = AdvancedVQC(config)
        vqc.fit(X, y)
        
        probas = vqc.predict_proba(X)
        
        assert probas.shape == (len(X), 2)
        assert np.allclose(probas.sum(axis=1), 1.0, atol=0.01)
        assert np.all((probas >= 0) & (probas <= 1))
    
    def test_vqc_with_zne(self):
        """Test VQC with ZNE error mitigation."""
        np.random.seed(42)
        X = np.random.randn(30, 16)
        y = (X[:, 0] > 0).astype(int)
        
        config = AdvancedConfig(
            n_epochs=3,
            error_mitigation="zne",
            noise_level=0.01,
            verbose=False
        )
        vqc = AdvancedVQC(config)
        
        assert vqc.zne is not None
        
        vqc.fit(X, y)
        score = vqc.score(X, y)
        
        assert 0.0 <= score <= 1.0
    
    def test_vqc_history_tracking(self):
        """Test training history is tracked."""
        np.random.seed(42)
        X = np.random.randn(40, 16)
        y = (X[:, 0] > 0).astype(int)
        
        config = AdvancedConfig(n_epochs=10, verbose=False)
        vqc = AdvancedVQC(config)
        vqc.fit(X, y)
        
        assert len(vqc.history['loss']) == 10
        assert len(vqc.history['accuracy']) == 10


class TestIntegration:
    """Integration tests for complete workflow."""
    
    def test_end_to_end_workflow(self, tmp_path):
        """Test complete workflow from data loading to results."""
        # Load dataset
        loader = DeepChemDatasetLoader(verbose=False)
        dataset = loader.load_dataset("BACE", max_samples=50)
        
        # Configure VQC
        config = AdvancedConfig(
            n_epochs=5,
            n_qubits=4,
            n_layers=1,
            results_dir=str(tmp_path / "results"),
            verbose=False
        )
        
        # Train and evaluate
        vqc = AdvancedVQC(config)
        vqc.fit(dataset['X_train'], dataset['y_train'])
        
        train_acc = vqc.score(dataset['X_train'], dataset['y_train'])
        test_acc = vqc.score(dataset['X_test'], dataset['y_test'])
        
        # Verify results
        assert 0.0 <= train_acc <= 1.0
        assert 0.0 <= test_acc <= 1.0
        
        # Export complexity report
        report_path = tmp_path / "complexity.md"
        vqc.complexity_analyzer.export_report(str(report_path))
        assert report_path.exists()
    
    def test_multiple_datasets(self):
        """Test training on multiple datasets."""
        loader = DeepChemDatasetLoader(verbose=False)
        
        results = []
        for dataset_name in ["BACE", "HIV"]:
            dataset = loader.load_dataset(dataset_name, max_samples=30)
            
            config = AdvancedConfig(n_epochs=3, verbose=False)
            vqc = AdvancedVQC(config)
            vqc.fit(dataset['X_train'], dataset['y_train'])
            
            acc = vqc.score(dataset['X_test'], dataset['y_test'])
            results.append({"dataset": dataset_name, "accuracy": acc})
        
        assert len(results) == 2
        assert all(0.0 <= r["accuracy"] <= 1.0 for r in results)


def run_tests():
    """Run all tests directly."""
    print("="*70)
    print("Running Advanced Quantum Framework V8.0 Tests")
    print("="*70)
    
    # Run pytest
    import pytest
    exit_code = pytest.main([__file__, "-v", "--tb=short"])
    
    return exit_code


if __name__ == "__main__":
    exit(run_tests())
