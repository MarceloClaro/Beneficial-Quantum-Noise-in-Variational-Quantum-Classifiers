#!/usr/bin/env python3
"""
Execution Manifest Generator

Purpose: Generate comprehensive execution manifests for reproducibility and auditability.

Usage:
    from tools.manifest_generator import ManifestGenerator
    
    mg = ManifestGenerator()
    mg.start_execution("20251227_001", "pennylane")
    # ... run experiments ...
    mg.end_execution(status="success", outputs={"metrics": "path/to/metrics.csv"})
    mg.save("manifests/pennylane/20251227_001/manifest_execucao.json")
"""

import json
import subprocess
import platform
import sys
import os
from datetime import datetime
from pathlib import Path
import importlib.metadata

class ManifestGenerator:
    """Generate execution manifests with full environment tracking."""
    
    def __init__(self):
        """Initialize manifest generator."""
        self.manifest = {
            "manifest_version": "1.0",
            "generation_timestamp": datetime.now().isoformat()
        }
        self.start_time = None
        self.end_time = None
    
    def start_execution(self, run_id, framework, config_path=None):
        """
        Start tracking an execution.
        
        Args:
            run_id: Unique identifier for this run
            framework: Framework name (pennylane, qiskit, cirq)
            config_path: Path to configuration file
        """
        self.start_time = datetime.now()
        
        self.manifest.update({
            "run_id": run_id,
            "framework": framework,
            "timestamp_start": self.start_time.isoformat(),
            "config_file": str(config_path) if config_path else None
        })
        
        # Capture environment
        self._capture_environment()
        self._capture_git_info()
        self._capture_dependencies()
        
        return self
    
    def end_execution(self, status="success", outputs=None, error_message=None):
        """
        End tracking an execution.
        
        Args:
            status: Execution status (success, failure, partial)
            outputs: Dictionary of output files
            error_message: Error message if status is failure
        """
        self.end_time = datetime.now()
        duration = (self.end_time - self.start_time).total_seconds()
        
        self.manifest.update({
            "timestamp_end": self.end_time.isoformat(),
            "duration_seconds": duration,
            "status": status,
            "outputs": outputs or {},
            "error_message": error_message
        })
        
        return self
    
    def add_seeds(self, seeds_dict):
        """
        Add seed information.
        
        Args:
            seeds_dict: Dictionary of seeds (e.g., {"global": 42, "numpy": 42})
        """
        self.manifest["seeds"] = seeds_dict
        return self
    
    def add_config_params(self, params):
        """
        Add configuration parameters.
        
        Args:
            params: Dictionary of configuration parameters
        """
        self.manifest["config_params"] = params
        return self
    
    def add_command(self, command):
        """
        Add command executed.
        
        Args:
            command: Command string or list
        """
        if isinstance(command, list):
            command = " ".join(command)
        self.manifest["command"] = command
        return self
    
    def _capture_environment(self):
        """Capture system environment information."""
        self.manifest["environment"] = {
            "os": platform.system(),
            "os_version": platform.version(),
            "os_release": platform.release(),
            "python_version": sys.version.split()[0],
            "python_implementation": platform.python_implementation(),
            "architecture": platform.machine(),
            "processor": platform.processor() or "Unknown",
            "hostname": platform.node()
        }
        
        # Try to get CPU and RAM info
        try:
            import psutil
            self.manifest["environment"]["cpu_count"] = psutil.cpu_count()
            self.manifest["environment"]["ram_total_gb"] = round(psutil.virtual_memory().total / (1024**3), 2)
        except ImportError:
            pass
    
    def _capture_git_info(self):
        """Capture git repository information."""
        try:
            # Get git commit hash
            commit = subprocess.check_output(
                ["git", "rev-parse", "HEAD"],
                stderr=subprocess.DEVNULL
            ).decode().strip()
            
            # Get git branch
            branch = subprocess.check_output(
                ["git", "rev-parse", "--abbrev-ref", "HEAD"],
                stderr=subprocess.DEVNULL
            ).decode().strip()
            
            # Check for uncommitted changes
            status = subprocess.check_output(
                ["git", "status", "--porcelain"],
                stderr=subprocess.DEVNULL
            ).decode().strip()
            
            self.manifest["git"] = {
                "commit": commit,
                "commit_short": commit[:7],
                "branch": branch,
                "uncommitted_changes": len(status) > 0,
                "dirty": len(status) > 0
            }
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.manifest["git"] = {
                "available": False,
                "reason": "Git not available or not a git repository"
            }
    
    def _capture_dependencies(self):
        """Capture Python dependencies with versions."""
        dependencies = {}
        
        # Key packages to track
        key_packages = [
            "pennylane",
            "qiskit",
            "cirq",
            "numpy",
            "scipy",
            "torch",
            "tensorflow",
            "jax",
            "matplotlib",
            "pandas",
            "scikit-learn"
        ]
        
        for package in key_packages:
            try:
                version = importlib.metadata.version(package)
                dependencies[package] = version
            except importlib.metadata.PackageNotFoundError:
                dependencies[package] = "not installed"
        
        self.manifest["dependencies"] = dependencies
    
    def save(self, filepath):
        """
        Save manifest to JSON file.
        
        Args:
            filepath: Path to save manifest
        """
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        
        with open(filepath, 'w') as f:
            json.dump(self.manifest, f, indent=2)
        
        print(f"Manifest saved to: {filepath}")
        return filepath
    
    def to_dict(self):
        """Return manifest as dictionary."""
        return self.manifest
    
    def to_json(self):
        """Return manifest as JSON string."""
        return json.dumps(self.manifest, indent=2)


def generate_manifest_for_execution(run_id, framework, config_path, seeds, status="success", outputs=None):
    """
    Convenience function to generate a complete manifest.
    
    Args:
        run_id: Unique run identifier
        framework: Framework name
        config_path: Path to config file
        seeds: Dictionary of seeds
        status: Execution status
        outputs: Dictionary of output files
    
    Returns:
        Path to saved manifest
    """
    mg = ManifestGenerator()
    mg.start_execution(run_id, framework, config_path)
    mg.add_seeds(seeds)
    mg.add_command(f"python framework_{framework}.py --config {config_path}")
    mg.end_execution(status=status, outputs=outputs)
    
    manifest_path = f"manifests/{framework}/{run_id}/manifest_execucao.json"
    return mg.save(manifest_path)


if __name__ == "__main__":
    # Example usage
    print("Manifest Generator - Example Usage\n")
    
    mg = ManifestGenerator()
    mg.start_execution("20251227_example", "pennylane", "configs/experiment_unified.yaml")
    mg.add_seeds({"global": 42, "numpy": 42, "random": 42})
    mg.add_command("python framework_investigativo_completo.py --config configs/experiment_unified.yaml")
    
    # Simulate some execution time
    import time
    time.sleep(1)
    
    mg.end_execution(
        status="success",
        outputs={
            "metrics": "results/pennylane/20251227_example/metrics.csv",
            "summary": "results/pennylane/20251227_example/summary.csv"
        }
    )
    
    # Print manifest
    print(mg.to_json())
    
    # Save manifest
    mg.save("manifests/example/manifest_execucao.json")
