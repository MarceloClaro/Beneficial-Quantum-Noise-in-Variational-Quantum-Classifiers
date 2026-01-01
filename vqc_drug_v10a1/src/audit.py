"""
SHA-256, checksums, pre-registration, seeds, environment snapshot – QUALIS A1 audit trail
"""
import os
import json
import datetime
import hashlib
import logging
import subprocess
import sys
from typing import Dict, Any

GLOBAL_SEEDS = [42, 123, 456, 789, 999]
LOGGER = logging.getLogger("VQC-Audit")

def get_seed(idx: int = None) -> int:
    """Round-robin seeds for cross-runs."""
    if idx is None:
        idx = int(datetime.datetime.utcnow().timestamp()) % len(GLOBAL_SEEDS)
    return GLOBAL_SEEDS[idx]

def checksum_file(path: str, out: str = "checksums.sha256", append: bool = True):
    """Compute and store SHA-256 of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            h.update(chunk)
    digest = h.hexdigest()
    mode = "a" if append else "w"
    with open(out, mode) as f:
        f.write(f"{digest}  {path}\n")
    LOGGER.debug("Checksum %s → %s", path, digest)
    return digest

def checksum_dir(path: str = ".", out: str = "checksums_final.sha256"):
    """Recurse and checksum all files (exclude itself)."""
    with open(out, "w") as f:
        for root, _, files in os.walk(path):
            for file in files:
                if "checksums" in file:
                    continue
                full = os.path.join(root, file)
                if os.path.isfile(full):
                    h = hashlib.sha256()
                    with open(full, "rb") as g:
                        for chunk in iter(lambda: g.read(4096), b""):
                            h.update(chunk)
                    f.write(f"{h.hexdigest()}  {full}\n")
    LOGGER.info("Final checksums saved → %s", out)

def pre_register(target: str, n_trials: int, alpha: float, power: float,
                 primary_endpoint: str) -> str:
    """Gera protocolo pré-registrado com hash SHA-256."""
    reg = {
        "protocol_id": "VQC-A1-2025-001",
        "timestamp": datetime.datetime.utcnow().isoformat(),
        "target": target,
        "n_trials": n_trials,
        "alpha": alpha,
        "power": power,
        "primary_endpoint": primary_endpoint,
        "seeds": GLOBAL_SEEDS,
        "multiple_comparisons": "Westfall-Young step-down",
        "effect_size_threshold": 0.35,
        "stopping_rules": {
            "futility": "upper 95 % CI < 0.01",
            "efficacy": "lower 95 % CI > 0.02"
        },
        "hash": None
    }
    reg_str = json.dumps(reg, sort_keys=True, separators=(',', ':'))
    digest = hashlib.sha256(reg_str.encode()).hexdigest()
    reg["hash"] = digest
    fname = f"01_protocolo_pre_registrado_{digest[:8]}.json"
    with open(fname, "w") as f:
        json.dump(reg, f, indent=2)
    LOGGER.info("Pre-registration SHA-256: %s", digest)
    return digest

def snapshot_environment(out: str = "environment_snapshot.txt"):
    """Captures pip + conda for full reproducibility."""
    pip = subprocess.check_output([sys.executable, "-m", "pip", "freeze"]).decode()
    try:
        conda = subprocess.check_output(["conda", "list", "--export"]).decode()
    except FileNotFoundError:
        conda = "# conda not found\n"
    with open(out, "w") as f:
        f.write("===== PIP =====\n" + pip + "\n===== CONDA =====\n" + conda)
    checksum_file(out)
    LOGGER.info("Environment snapshot → %s", out)

def power_curve(effect=0.35, alpha=0.05, power_min=0.8):
    """Return N per arm for effect size (simplified)."""
    from statsmodels.stats.power import tt_solve_power
    n = int(tt_solve_power(effect_size=effect, alpha=alpha, power=power_min, ratio=1))
    return {effect: n}
