"""
VQC-Molecular v10.1-A1
Mathematically-optimised quantum drug screening – QUALIS A1 compliant
"""
__version__ = "10.1.0"
__author__ = "Marcelo Claro Laranjeira"
__email__ = "marceloclaro@gmail.com"
__all__ = [
    "data", "models", "tune", "audit", "plots", "cli",
    "VQCAudit", "run_study", "power_curve", "get_seed", "LindbladOptimalSchedule",
    "v10_noise_ansatz", "v10_fisher", "v10_lindblad", "v10_qasr", "v10_meta", "v10_stats",
    "v10_zne", "v10_psqk", "v10_qng", "v10_qda"
]

# Facade imports – user only needs "from vqc_drug_a1 import run_study"
from .data import load_split, MorganFeaturizer, PCAReducer
from .models import VQCAudit
from .tune import run_study, LindbladOptimalSchedule, QASRSampler
from .audit import get_seed, power_curve, checksum_file, pre_register
from .plots import fig_auc_heatmap, fig_forest_effect, latex_boilerplate
from .cli import main  # entry point
from .v10_fisher import choose_constant_with_max_fisher
from .v10_noise_ansatz import extract_backend_noise, noise_aware_ansatz
from .v10_lindblad import load_backend_schedule
