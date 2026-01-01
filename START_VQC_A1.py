#!/usr/bin/env python
"""
Inicializador VQC-Molecular v8.1-A1
Fornece acesso direto ao framework desde a pasta raiz.

Uso:
    python -m vqc_molecular_v8_1_a1.run_vqc_a1 --target EGFR --trials 300
"""

import sys
from pathlib import Path

# Adicionar pasta vqc_molecular_v8_1_a1 ao path
vqc_dir = Path(__file__).parent / "vqc_molecular_v8_1_a1"
sys.path.insert(0, str(vqc_dir))

print("\n" + "="*70)
print("🔬 VQC-Molecular v8.1-A1 - Framework Inicializado")
print("="*70)
print(f"\n📂 Diretório: {vqc_dir}")
print(f"\n💡 Dica: Execute o seguinte comando:")
print(f"   cd vqc_molecular_v8_1_a1/")
print(f"   python run_vqc_a1.py --target EGFR --trials 300")
print("\n" + "="*70 + "\n")
