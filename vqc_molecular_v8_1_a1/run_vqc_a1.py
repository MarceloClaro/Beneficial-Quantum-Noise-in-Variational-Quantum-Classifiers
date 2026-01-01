#!/usr/bin/env python
"""
VQC-Molecular v8.1-A1 - Quantum Drug Discovery with Qualis A1 Compliance
==========================================================================

Ponto de entrada simplificado para executar o framework.

Uso:
    python run_vqc_a1.py --target EGFR --trials 300
    python run_vqc_a1.py --target HIV --trials 500 --seed 42
"""

import sys
import os
from pathlib import Path

# Ensure the project directory is in the path
project_dir = Path(__file__).parent.absolute()
sys.path.insert(0, str(project_dir))

# Import the main framework
from vqc_drug_qualis_a1 import run_qualis_a1_study

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="VQC-Molecular v8.1-A1 - Quantum Drug Discovery Framework",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_vqc_a1.py --target EGFR --trials 300
  python run_vqc_a1.py --target HIV --trials 500 --seed 123
  python run_vqc_a1.py --target Malaria --trials 200 --n-qubits 8
        """
    )
    
    parser.add_argument(
        "--target",
        type=str,
        default="EGFR",
        choices=["EGFR", "HIV", "Malaria", "COVID"],
        help="Target molecular (default: EGFR)"
    )
    
    parser.add_argument(
        "--trials",
        type=int,
        default=300,
        help="Número de tentativas Optuna (default: 300)"
    )
    
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Seed para reprodutibilidade (default: 42)"
    )
    
    parser.add_argument(
        "--n-qubits",
        type=int,
        default=None,
        help="Número de qubits (default: auto-detectado)"
    )
    
    parser.add_argument(
        "--no-deepchem",
        action="store_true",
        help="Pular comparação DeepChem (mais rápido)"
    )
    
    args = parser.parse_args()
    
    print("\n" + "="*70)
    print("🔬 VQC-Molecular v8.1-A1 - Quantum Drug Discovery")
    print("="*70)
    print(f"\n📊 Configuração:")
    print(f"   Alvo: {args.target}")
    print(f"   Tentativas Optuna: {args.trials}")
    print(f"   Seed: {args.seed}")
    print(f"   Qubits: {args.n_qubits if args.n_qubits else 'Auto'}")
    print(f"   DeepChem: {'Desativado' if args.no_deepchem else 'Ativado'}")
    print("\n" + "="*70 + "\n")
    
    # Run the main study
    run_qualis_a1_study(
        target=args.target,
        n_trials=args.trials,
        seed=args.seed,
        n_qubits=args.n_qubits,
        include_deepchem=not args.no_deepchem
    )
    
    print("\n" + "="*70)
    print("✅ Estudo concluído com sucesso!")
    print("="*70 + "\n")
