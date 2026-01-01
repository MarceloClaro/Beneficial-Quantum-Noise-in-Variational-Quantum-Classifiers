# ============================================================================
# preregister.py
# Pré-registro automatizado com SHA-256 para conformidade Qualis A1
# ============================================================================
import json
import hashlib
import datetime
import os
from typing import Tuple, Dict, Any

def pre_register(target: str, n_trials: int, alpha: float = 0.05,
                 power: float = 0.8, primary_endpoint: str = "delta_AUC") -> Tuple[str, str]:
    """
    Cria pré-registro com hash SHA-256 imutável.
    
    Args:
        target: Nome do alvo (EGFR, HIV, Malaria, COVID)
        n_trials: Número de trials Optuna planejados
        alpha: Nível de significância (default 0.05)
        power: Potência estatística desejada (default 0.8)
        primary_endpoint: Endpoint primário (default delta_AUC)
    
    Returns:
        Tupla (caminho_arquivo, hash_SHA256)
    """
    
    reg = {
        "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
        "version": "VQC-Molecular v8.1-A1",
        "target": target,
        "primary_endpoint": primary_endpoint,
        "secondary_endpoints": [
            "VQC_ROC_AUC",
            "Baseline_ROC_AUC",
            "Cohen_d_effect_size",
            "Bootstrap_CI_95_percent"
        ],
        "study_design": {
            "type": "Optimization + Baseline Comparison",
            "n_trials_planned": n_trials,
            "cv_folds": 5,
            "cv_strategy": "StratifiedKFold",
            "random_seed": 42
        },
        "statistical_assumptions": {
            "alpha": alpha,
            "power": power,
            "alternative_hypothesis": "VQC_AUC > Baseline_AUC + 0.02",
            "effect_size_threshold": 0.3  # Cohen d
        },
        "hypotheses": [
            f"H0: delta_{primary_endpoint} ≤ 0.00 (non-inferiority)",
            f"H1: delta_{primary_endpoint} > 0.02 (superiority, clinical relevance)"
        ],
        "multiple_comparison_method": "Bonferroni-Holm",
        "effect_size_measures": [
            "Cohen_d (standardized mean difference)",
            "Hedges_g (unbiased Cohen d)",
            "Glass_delta (uses baseline SD)"
        ],
        "stopping_rules": {
            "futility": "upper_95_CI < 0.005",
            "efficacy": "lower_95_CI > 0.015",
            "early_stopping_enabled": True,
            "early_stopping_rounds": 50
        },
        "datasets": {
            "EGFR": {
                "n_molecules": 6847,
                "pct_active": 8,
                "source": "ChEMBL EGFR kinase inhibitors"
            },
            "HIV": {
                "n_molecules": 41913,
                "pct_active": 4,
                "source": "MoleculeNet HIV reverse transcriptase"
            },
            "Malaria": {
                "n_molecules": 13281,
                "pct_active": 6,
                "source": "MoleculeNet Plasmodium falciparum"
            },
            "COVID": {
                "n_molecules": 10427,
                "pct_active": 5,
                "source": "COVID-19 Moonshot SARS-CoV-2 protease"
            }
        },
        "vqc_hyperparameter_search_space": {
            "n_qubits": {"min": 4, "max": 20, "step": 2},
            "n_layers": {"min": 1, "max": 8},
            "noise_type": ["depolarizing", "amplitude_damping", "phase_damping", "none"],
            "noise_level": {"min": 0.0005, "max": 0.02, "scale": "log"},
            "learning_rate": {"min": 0.0005, "max": 0.2, "scale": "log"},
            "epochs": {"min": 10, "max": 50},
            "batch_size": [16, 32, 64]
        },
        "baseline_model": {
            "name": "DeepChem GraphConvModel",
            "version": "4.6.0",
            "task": "Binary classification",
            "epochs": 40,
            "batch_size": 50
        },
        "reproducibility": {
            "global_seed": 42,
            "np_random_seed": 42,
            "tf_random_seed": 42,
            "python_version": "3.10+",
            "major_dependencies": [
                "pennylane >= 0.32",
                "optuna >= 3.4",
                "deepchem >= 4.6",
                "scikit-learn >= 1.3",
                "numpy >= 1.24",
                "pandas >= 2.0"
            ]
        },
        "data_availability": {
            "raw_data": "Public (QSAR datasets auto-downloaded)",
            "cache_directory": "qsar_cache/",
            "data_version_control": "SHA-256 checksums maintained"
        },
        "ethics": {
            "irb_approval": "Not required (public computational data, no human subjects)",
            "data_privacy": "All data from public repositories",
            "license": "MIT (code), CC-BY (data)"
        },
        "preregistration_timestamp": datetime.datetime.utcnow().isoformat() + "Z",
        "preregistration_status": "LOCKED"
    }
    
    # Serialize com ordem determinística para hash consistente
    reg_json = json.dumps(reg, sort_keys=True, separators=(',', ':'), indent=2)
    
    # Calcular SHA-256
    sha256_hash = hashlib.sha256(reg_json.encode('utf-8')).hexdigest()
    
    # Criar diretório de resultados com timestamp
    timestamp = datetime.datetime.utcnow().strftime("%Y-%m-%d_%H-%M-%S")
    results_dir = f"results_{timestamp}"
    os.makedirs(results_dir, exist_ok=True)
    
    # Salvar protocolo pré-registrado
    proto_filename = os.path.join(results_dir, 
                                  f"01_protocolo_pre_registrado_{sha256_hash[:8]}.json")
    
    with open(proto_filename, 'w', encoding='utf-8') as f:
        f.write(reg_json)
    
    print(f"\n✅ PRÉ-REGISTRO CRIADO:")
    print(f"   Arquivo: {proto_filename}")
    print(f"   SHA-256: {sha256_hash}")
    print(f"   Status:  BLOQUEADO (imutável)\n")
    
    return proto_filename, sha256_hash


def validate_preregistration(proto_file: str) -> bool:
    """Valida integridade do pré-registro comparando hash."""
    with open(proto_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    current_hash = hashlib.sha256(content.encode('utf-8')).hexdigest()
    stored_hash = os.path.basename(proto_file).split('_')[-1].split('.')[0]
    
    is_valid = current_hash.startswith(stored_hash)
    status = "✅ VÁLIDO" if is_valid else "❌ CORROMPIDO"
    print(f"{status}: {os.path.basename(proto_file)}")
    
    return is_valid


if __name__ == "__main__":
    # Exemplo
    proto_file, proto_hash = pre_register(
        target="EGFR",
        n_trials=500,
        alpha=0.05,
        power=0.8,
        primary_endpoint="delta_AUC_VQC_vs_GraphConv"
    )
    validate_preregistration(proto_file)
