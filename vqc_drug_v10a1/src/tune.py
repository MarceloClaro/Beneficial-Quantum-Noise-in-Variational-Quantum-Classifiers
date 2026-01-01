"""
Mathematically-optimised quantum hyper-parameter search (v10.1-A1).
- QASR bandit (UCB) com penalização de barren plateau
- Lindblad-optimal schedule registrado no trial
- Fisher-CRLB constant picker
- Meta-learning warm-start
"""
import numpy as np
import optuna
import torch
from optuna.pruners import MedianPruner
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from .models import VQCAudit
from .audit import get_seed
from .v10_qasr import QASRSampler
from .v10_lindblad import LindbladOptimalSchedule
from .v10_fisher import choose_constant_with_max_fisher, PHYSICAL_CONSTANTS
from .v10_meta import meta_warm_start
from .v10_stats import qualis_report_from_trials

# Hyper-parameter space
HYPER_SPACE = {
    "n_qubits": [4, 6, 8, 10, 12, 14, 16, 18, 20],
    "n_layers": [1, 2, 3, 4, 5, 6],
    "noise_type": ["none", "depolarizing", "amplitude_damping"],
    "noise_level": [0.0, 0.001, 0.005, 0.01, 0.015, 0.02],
    "constant_init": ["random", "pi", "e", "phi", "hbar", "alpha", "fisher"],
    "arch": ["tree", "star", "brickwork"],
    "optimizer": ["adam", "sgd"],
    "loss": ["mse", "bce", "hinge"],
    "lr": [1e-3, 5e-3, 1e-2, 2e-2],
    "batch_size": [16, 32, 64],
}

def objective_ultra(trial, X, y, folds, target: str):
    """Single trial with cross-validation e rastreamento de métricas A1."""
    # Sample hyperparameters
    cfg = {k: trial.suggest_categorical(k, v) for k, v in HYPER_SPACE.items()}
    cfg["epochs"] = 40

    if cfg.get("constant_init") == "fisher":
        best_const, _ = choose_constant_with_max_fisher(PHYSICAL_CONSTANTS)
        cfg["constant_init"] = best_const

    cfg = meta_warm_start(target, cfg)

    # Cross-validation
    aucs = []
    schedule = LindbladOptimalSchedule()
    for fold, (tr_idx, val_idx) in enumerate(folds):
        X_tr, X_val = X[tr_idx], X[val_idx]
        y_tr, y_val = y[tr_idx], y[val_idx]
        X_tr = torch.tensor(X_tr, dtype=torch.float32)
        y_tr = torch.tensor(y_tr, dtype=torch.float32)
        X_val = torch.tensor(X_val, dtype=torch.float32)
        y_val = torch.tensor(y_val, dtype=torch.float32)

        model = VQCAudit(**cfg, trial_id=trial.number)
        model.fit(X_tr, y_tr, X_val, y_val)

        proba = model.predict_proba(X_val)[:, 1]
        aucs.append(roc_auc_score(y_val.cpu().numpy(), proba))

    mean_auc = np.mean(aucs)
    trial.set_user_attr("std_auc", np.std(aucs))
    trial.set_user_attr("aucs", aucs)
    trial.set_user_attr("lindblad", schedule.summary(cfg["epochs"]))
    trial.set_user_attr("meta_used", cfg.get("meta_used", False))
    trial.set_user_attr("constant_init", cfg.get("constant_init"))
    return mean_auc

def run_study(X, y, target: str, n_trials=500, max_qubits=20):
    """Full pipeline v10."""
    folds = list(StratifiedKFold(n_splits=5, shuffle=True, 
                                  random_state=get_seed()).split(X, y))
    
    sampler = QASRSampler(seed=get_seed())
    
    study = optuna.create_study(
        direction="maximize",
        sampler=sampler,
        pruner=MedianPruner(n_warmup_steps=10),
        study_name=f"{target}_v10"
    )
    
    study.optimize(lambda trial: objective_ultra(trial, X, y, folds, target),
                   n_trials=n_trials, show_progress_bar=True)

    study.set_user_attr("qualis_report", qualis_report_from_trials(study.trials_dataframe()).to_dict())
    return study
