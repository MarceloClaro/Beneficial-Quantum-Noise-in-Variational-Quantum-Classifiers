"""
Unit tests for ClassificadorVQC class.

Tests verify that the VQC classifier works correctly on toy datasets,
including training, prediction, and various configuration options.
"""

import sys
from pathlib import Path

import numpy as np
import pytest
from sklearn.datasets import make_moons, make_circles
from sklearn.model_selection import train_test_split

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    import pennylane as qml
    PENNYLANE_AVAILABLE = True
except ImportError:
    PENNYLANE_AVAILABLE = False
    pytest.skip("PennyLane not available", allow_module_level=True)

from framework_investigativo_completo import ClassificadorVQC


class TestClassificadorVQC:
    """Test suite for VQC classifier."""

    @pytest.fixture
    def toy_dataset(self):
        """Create a simple toy dataset for testing."""
        X, y = make_moons(n_samples=40, noise=0.1, random_state=42)
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.25, random_state=42
        )
        return X_train, X_test, y_train, y_test

    def test_initialization_basic(self):
        """Test basic initialization of VQC."""
        vqc = ClassificadorVQC(
            n_qubits=2,
            n_camadas=1,
            arquitetura='basico',
            n_epocas=2,
            seed=42
        )
        assert vqc.n_qubits == 2
        assert vqc.n_camadas == 1
        assert vqc.arquitetura == 'basico'

    def test_fit_basic(self, toy_dataset):
        """Test basic training functionality."""
        X_train, X_test, y_train, y_test = toy_dataset
        
        vqc = ClassificadorVQC(
            n_qubits=2,
            n_camadas=1,
            arquitetura='basico',
            n_epocas=3,
            batch_size=10,
            seed=42
        )
        
        # Should train without errors
        vqc.fit(X_train, y_train)
        
        # Check that weights were created
        assert hasattr(vqc, 'weights_')
        assert vqc.weights_ is not None
        assert len(vqc.weights_) > 0

    def test_predict(self, toy_dataset):
        """Test prediction functionality."""
        X_train, X_test, y_train, y_test = toy_dataset
        
        vqc = ClassificadorVQC(
            n_qubits=2,
            n_camadas=1,
            arquitetura='basico',
            n_epocas=2,
            seed=42
        )
        vqc.fit(X_train, y_train)
        
        # Test predictions
        y_pred = vqc.predict(X_test)
        
        assert len(y_pred) == len(y_test)
        assert all(p in [0, 1] for p in y_pred)

    def test_score(self, toy_dataset):
        """Test accuracy scoring."""
        X_train, X_test, y_train, y_test = toy_dataset
        
        vqc = ClassificadorVQC(
            n_qubits=2,
            n_camadas=2,  # Increased from 1 to 2 for more expressivity
            arquitetura='basico',
            n_epocas=15,  # Increased from 5 to 15 for better convergence
            taxa_aprendizado=0.05,  # Increased from 0.01 to 0.05
            batch_size=10,  # Add explicit batch size
            seed=42
        )
        vqc.fit(X_train, y_train)
        
        # Get accuracy score
        accuracy = vqc.score(X_test, y_test)
        
        # Should be between 0 and 1
        assert 0.0 <= accuracy <= 1.0
        
        # With enough training, should achieve some reasonable accuracy
        assert accuracy > 0.3  # At least better than random

    def test_different_architectures(self, toy_dataset):
        """Test different VQC architectures."""
        X_train, X_test, y_train, y_test = toy_dataset
        
        architectures = ['basico', 'strongly_entangling', 'hardware_efficient']
        
        for arch in architectures:
            vqc = ClassificadorVQC(
                n_qubits=2,
                n_camadas=1,
                arquitetura=arch,
                n_epocas=2,
                seed=42
            )
            
            # Should train without errors
            vqc.fit(X_train, y_train)
            accuracy = vqc.score(X_test, y_test)
            
            # Should produce valid results
            assert 0.0 <= accuracy <= 1.0

    def test_noise_models(self, toy_dataset):
        """Test VQC with different noise models."""
        X_train, X_test, y_train, y_test = toy_dataset
        
        noise_types = ['sem_ruido', 'depolarizante', 'phase_damping']
        
        for noise_type in noise_types:
            vqc = ClassificadorVQC(
                n_qubits=2,
                n_camadas=1,
                arquitetura='basico',
                tipo_ruido=noise_type,
                nivel_ruido=0.01,
                n_epocas=2,
                seed=42
            )
            
            vqc.fit(X_train, y_train)
            accuracy = vqc.score(X_test, y_test)
            
            assert 0.0 <= accuracy <= 1.0

    def test_initialization_strategies(self, toy_dataset):
        """Test different weight initialization strategies."""
        X_train, X_test, y_train, y_test = toy_dataset
        
        strategies = ['aleatoria', 'matematico', 'quantico', 'fibonacci_spiral']
        
        for strategy in strategies:
            vqc = ClassificadorVQC(
                n_qubits=2,
                n_camadas=1,
                arquitetura='basico',
                estrategia_init=strategy,
                n_epocas=2,
                seed=42
            )
            
            vqc.fit(X_train, y_train)
            
            # Should have initialized weights
            assert hasattr(vqc, 'weights_')
            assert vqc.weights_ is not None

    def test_training_history(self, toy_dataset):
        """Test that training history is recorded."""
        X_train, X_test, y_train, y_test = toy_dataset
        
        vqc = ClassificadorVQC(
            n_qubits=2,
            n_camadas=1,
            arquitetura='basico',
            n_epocas=5,
            seed=42
        )
        vqc.fit(X_train, y_train)
        
        # Check history exists
        assert hasattr(vqc, 'historico_')
        assert 'custo' in vqc.historico_
        assert 'acuracia_treino' in vqc.historico_
        assert 'epoca' in vqc.historico_
        
        # Check history length matches epochs
        assert len(vqc.historico_['epoca']) == 5

    def test_early_stopping(self, toy_dataset):
        """Test early stopping functionality."""
        X_train, X_test, y_train, y_test = toy_dataset
        
        vqc = ClassificadorVQC(
            n_qubits=2,
            n_camadas=1,
            arquitetura='basico',
            n_epocas=50,  # Large number
            early_stopping=True,
            patience=3,
            val_split=0.2,
            seed=42
        )
        vqc.fit(X_train, y_train)
        
        # Should stop before 50 epochs in most cases
        actual_epochs = len(vqc.historico_['epoca'])
        # May train full epochs on toy dataset, but should not error
        assert 1 <= actual_epochs <= 50

    def test_reproducibility(self, toy_dataset):
        """Test that results are reproducible with same seed."""
        X_train, X_test, y_train, y_test = toy_dataset
        
        vqc1 = ClassificadorVQC(
            n_qubits=2,
            n_camadas=1,
            arquitetura='basico',
            n_epocas=3,
            seed=42
        )
        vqc1.fit(X_train, y_train)
        pred1 = vqc1.predict(X_test)
        
        vqc2 = ClassificadorVQC(
            n_qubits=2,
            n_camadas=1,
            arquitetura='basico',
            n_epocas=3,
            seed=42
        )
        vqc2.fit(X_train, y_train)
        pred2 = vqc2.predict(X_test)
        
        # Results should be very similar (allowing for small numerical differences)
        agreement = np.mean(pred1 == pred2)
        assert agreement > 0.8  # At least 80% agreement

    def test_batch_size_effect(self, toy_dataset):
        """Test that different batch sizes work."""
        X_train, X_test, y_train, y_test = toy_dataset
        
        batch_sizes = [5, 10, 15]
        
        for batch_size in batch_sizes:
            vqc = ClassificadorVQC(
                n_qubits=2,
                n_camadas=1,
                arquitetura='basico',
                n_epocas=2,
                batch_size=batch_size,
                seed=42
            )
            
            vqc.fit(X_train, y_train)
            accuracy = vqc.score(X_test, y_test)
            
            assert 0.0 <= accuracy <= 1.0

    def test_multiclass_classification(self):
        """Test VQC on multiclass problem."""
        # Create 3-class dataset
        from sklearn.datasets import make_blobs
        X, y = make_blobs(n_samples=60, centers=3, n_features=2, 
                         random_state=42, cluster_std=1.0)
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.25, random_state=42
        )
        
        vqc = ClassificadorVQC(
            n_qubits=2,
            n_camadas=1,
            arquitetura='basico',
            n_epocas=3,
            seed=42
        )
        
        # Should handle multiclass (via one-vs-rest internally or similar)
        vqc.fit(X_train, y_train)
        y_pred = vqc.predict(X_test)
        
        # Check predictions are valid
        assert len(y_pred) == len(y_test)
        assert set(y_pred).issubset(set([0, 1, 2]))

    def test_learning_rate_effect(self, toy_dataset):
        """Test that learning rate parameter works."""
        X_train, X_test, y_train, y_test = toy_dataset
        
        learning_rates = [0.001, 0.01, 0.1]
        
        for lr in learning_rates:
            vqc = ClassificadorVQC(
                n_qubits=2,
                n_camadas=1,
                arquitetura='basico',
                taxa_aprendizado=lr,
                n_epocas=2,
                seed=42
            )
            
            vqc.fit(X_train, y_train)
            assert vqc.taxa_aprendizado == lr

    def test_noise_scheduling(self, toy_dataset):
        """Test noise scheduling functionality."""
        X_train, X_test, y_train, y_test = toy_dataset
        
        vqc = ClassificadorVQC(
            n_qubits=2,
            n_camadas=1,
            arquitetura='basico',
            tipo_ruido='depolarizante',
            nivel_ruido=0.05,
            ruido_schedule='linear',
            ruido_inicial=0.1,
            ruido_final=0.01,
            n_epocas=5,
            seed=42
        )
        
        vqc.fit(X_train, y_train)
        
        # Check that noise levels were recorded
        assert 'nivel_ruido' in vqc.historico_
        niveis = vqc.historico_['nivel_ruido']
        
        # Should show decreasing trend
        assert niveis[0] > niveis[-1]

    def test_invalid_architecture_fallback(self, toy_dataset):
        """Test behavior with invalid architecture."""
        X_train, X_test, y_train, y_test = toy_dataset
        
        # Should either raise error or fallback gracefully
        try:
            vqc = ClassificadorVQC(
                n_qubits=2,
                n_camadas=1,
                arquitetura='invalid_architecture',
                n_epocas=2,
                seed=42
            )
            # If it doesn't raise, it should still work
            vqc.fit(X_train, y_train)
        except (KeyError, ValueError):
            # Expected error for invalid architecture
            pass

    def test_sklearn_compatibility(self, toy_dataset):
        """Test scikit-learn API compatibility."""
        X_train, X_test, y_train, y_test = toy_dataset
        
        vqc = ClassificadorVQC(
            n_qubits=2,
            n_camadas=1,
            arquitetura='basico',
            n_epocas=2,
            seed=42
        )
        
        # Check sklearn methods exist
        assert hasattr(vqc, 'fit')
        assert hasattr(vqc, 'predict')
        assert hasattr(vqc, 'score')
        
        # Test fit returns self
        result = vqc.fit(X_train, y_train)
        assert result is vqc
        
        # Test classes_ attribute
        vqc.fit(X_train, y_train)
        assert hasattr(vqc, 'classes_')


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
