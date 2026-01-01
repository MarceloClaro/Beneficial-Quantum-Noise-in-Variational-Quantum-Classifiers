"""Smoke tests para módulos v10.1-A1."""
import numpy as np
from src.v10_noise_ansatz import extract_backend_noise, ansatz_sha256
from src.v10_fisher import crlb_effect_size, choose_constant_with_max_fisher, PHYSICAL_CONSTANTS
from src.v10_lindblad import LindbladOptimalSchedule, load_backend_schedule
from src.v10_qasr import QASRSampler
from src.v10_meta import meta_warm_start
from src.v10_stats import westfall_young_step_down, qualis_report_from_trials
from src.v10_zne import richardson_extrapolate, should_trigger_zne
from src.v10_psqk import psqk_shared_weights
from src.v10_qng import qng_metric_block_diag, condition_number
from src.v10_qda import qda_augment_shadows, shadow_statistics
import pandas as pd


def test_noise_profile_hash():
    noise = extract_backend_noise("default")
    assert noise
    digest = ansatz_sha256(noise)
    assert len(digest) == 64


def test_fisher_picker():
    g1 = np.random.randn(20)
    g2 = np.random.randn(18)
    res = crlb_effect_size(g1, g2, n_boot=100)  # n_boot reduzido para teste
    assert "d_obs" in res and "fisher_lb" in res
    name, value = choose_constant_with_max_fisher(PHYSICAL_CONSTANTS, samples=16)
    assert name in PHYSICAL_CONSTANTS and np.isfinite(value)


def test_lindblad_schedule():
    sched = load_backend_schedule("default")
    lut = sched.lookup_table(5)
    assert len(lut) == 5
    summary = sched.summary(5)
    assert "mean" in summary


def test_qasr_sampler():
    sampler = QASRSampler()
    assert sampler is not None


def test_meta_warm_start():
    cfg = {"lr": 0.01, "n_layers": 2, "n_qubits": 8}
    out = meta_warm_start("HIV", cfg)
    assert out["meta_used"] is True


def test_stats_report():
    df = pd.DataFrame({"value": [0.8, 0.82, 0.78]})
    report = qualis_report_from_trials(df)
    assert not report.empty
    reject, padj = westfall_young_step_down(np.array([0.01, 0.02, 0.5]))
    assert len(padj) == 3


def test_zne_helpers():
    zero = richardson_extrapolate([0.7, 0.65, 0.6], [1.0, 1.5, 2.0])
    assert np.isfinite(zero)
    assert should_trigger_zne([1, 0.9, 0.85, 0.84, 0.84, 0.84]) is True


def test_psqk_share():
    theta = np.array([1.0, 2.0, 3.0, 4.0])
    shared, mask = psqk_shared_weights(theta, share_rate=0.5)
    assert mask.sum() == 2
    assert np.isfinite(shared).all()


def test_qng_metric():
    jac = np.ones((3, 4))
    metric = qng_metric_block_diag(jac, block_size=2)
    assert metric.shape == (4, 4)
    assert condition_number(metric) >= 0


def test_qda_aug():
    shadows = np.ones((2, 3))
    aug = qda_augment_shadows(shadows, n_aug=1, noise=0.0)
    assert aug.shape[0] == 4
    mean, std = shadow_statistics(aug)
    assert np.isfinite(mean) and np.isfinite(std)


if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-q"]) 
