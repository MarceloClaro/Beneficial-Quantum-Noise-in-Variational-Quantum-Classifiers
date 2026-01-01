"""
Script de instalação de dependências - VQC Framework
Instala dependências via pip e fornece instruções para RDKit (requer conda)
"""
import subprocess
import sys
from pathlib import Path

print("="*70)
print("🔧 INSTALAÇÃO DE DEPENDÊNCIAS - VQC FRAMEWORK")
print("="*70)
print()

# Caminho do requirements
req_file = Path(__file__).parent / "requirements.txt"

# Instalar via pip
print("📦 Instalando dependências via pip...")
print(f"   Arquivo: {req_file}")
print()

try:
    subprocess.run(
        [sys.executable, "-m", "pip", "install", "-r", str(req_file)],
        check=True
    )
    print()
    print("✅ Dependências pip instaladas com sucesso!")
    print()
except subprocess.CalledProcessError as e:
    print()
    print(f"⚠️  Algumas dependências falharam (código {e.returncode})")
    print()

# Instruções para RDKit
print("="*70)
print("⚠️  ATENÇÃO: RDKit requer Conda")
print("="*70)
print()
print("O RDKit não está disponível via pip no Windows/Python 3.13.")
print("Para instalar:")
print()
print("   1. Instale Miniconda/Anaconda (se ainda não tiver)")
print("   2. Execute em um prompt Anaconda:")
print()
print("      conda install -c conda-forge rdkit")
print()
print("   3. Ou crie um ambiente conda com RDKit:")
print()
print("      conda create -n vqc python=3.11 rdkit -c conda-forge")
print("      conda activate vqc")
print("      pip install -r requirements.txt")
print()
print("="*70)
print()

# Verificar instalações
print("🔍 Verificando módulos instalados...")
print()

modulos = [
    ("pennylane", "PennyLane"),
    ("numpy", "NumPy"),
    ("pandas", "Pandas"),
    ("optuna", "Optuna"),
    ("opentelemetry", "OpenTelemetry"),
    ("qiskit", "Qiskit"),
    ("cirq", "Cirq"),
]

sucesso = 0
falhas = []

for modulo, nome in modulos:
    try:
        __import__(modulo)
        print(f"   ✅ {nome}")
        sucesso += 1
    except ImportError:
        print(f"   ❌ {nome} - NÃO INSTALADO")
        falhas.append(nome)

# RDKit separado
try:
    __import__("rdkit")
    print(f"   ✅ RDKit")
    sucesso += 1
except ImportError:
    print(f"   ⚠️  RDKit - Requer conda (veja instruções acima)")
    falhas.append("RDKit (requer conda)")

print()
print("="*70)
print(f"📊 RESUMO: {sucesso}/{len(modulos)+1} módulos instalados")
print("="*70)

if falhas:
    print()
    print("⚠️  Módulos faltantes:")
    for f in falhas:
        print(f"   - {f}")
    print()

print()
print("💡 Para executar os frameworks:")
print()
print("   Sem RDKit (QAOA, VQC básico):")
print("      python pipeline_com_tracing.py")
print()
print("   Com RDKit (VQC molecular, HIV):")
print("      conda activate vqc  # (ambiente com RDKit)")
print("      python execucao_rapida_hiv.py")
print()
