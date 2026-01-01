"""Pytest suite – 30 s smoke test."""
import pytest
import sys
from pathlib import Path

def test_imports():
    """Test all module imports."""
    from src import data, models, tune, audit, plots, cli
    assert data is not None
    assert models is not None
    assert tune is not None

def test_data_loading():
    """Test data loading pipeline (small sample)."""
    from src.data import MorganFeaturizer, PCAReducer
    import numpy as np
    
    feat = MorganFeaturizer(radius=2, n_bits=1024)
    fp = feat.featurize("CCO")  # ethanol
    assert fp.shape == (1024,)
    assert fp.dtype == np.float32

def test_model_creation():
    """Test VQC model instantiation."""
    from src.models import VQCAudit
    
    model = VQCAudit(
        n_qubits=4, n_layers=2, noise_type="none",
        noise_level=0.0, lr=0.01, epochs=1,
        constant_init="pi", arch="tree",
        optimizer="adam", loss="bce",
        batch_size=16, trial_id=0
    )
    
    assert model.n_qubits == 4
    assert model.n_layers == 2

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
