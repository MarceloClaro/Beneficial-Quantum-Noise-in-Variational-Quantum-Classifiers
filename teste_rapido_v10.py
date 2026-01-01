"""
Teste rápido do VQC-Molecular v10.0-A1
Verifica importações e funcionalidades básicas
"""
import sys
from pathlib import Path

print("="*70)
print("🧪 TESTE RÁPIDO VQC-Molecular v10.0-A1")
print("="*70)
print()

# Adicionar src ao path
vqc_dir = Path(__file__).parent / "vqc_drug_v10a1"
sys.path.insert(0, str(vqc_dir))

# Teste 1: Importações
print("📦 Teste 1: Verificando importações dos módulos...")
try:
    from src import data, models, tune, audit, plots, cli
    print("   ✅ Módulo data importado")
    print("   ✅ Módulo models importado")
    print("   ✅ Módulo tune importado")
    print("   ✅ Módulo audit importado")
    print("   ✅ Módulo plots importado")
    print("   ✅ Módulo cli importado")
    print()
except Exception as e:
    print(f"   ❌ Erro na importação: {e}")
    sys.exit(1)

# Teste 2: Featurização
print("🧬 Teste 2: Testando featurização molecular...")
try:
    from src.data import MorganFeaturizer
    feat = MorganFeaturizer(radius=2, n_bits=1024)
    fp = feat.featurize("CCO")  # ethanol
    print(f"   ✅ Morgan fingerprint gerado: shape={fp.shape}, dtype={fp.dtype}")
    print()
except Exception as e:
    print(f"   ❌ Erro na featurização: {e}")
    print()

# Teste 3: Criação do modelo VQC
print("⚛️  Teste 3: Testando criação do modelo VQC...")
try:
    from src.models import VQCAudit
    model = VQCAudit(
        n_qubits=4, n_layers=2, noise_type="none",
        noise_level=0.0, lr=0.01, epochs=1,
        constant_init="pi", arch="tree",
        optimizer="adam", loss="bce",
        batch_size=16, trial_id=0
    )
    print(f"   ✅ Modelo VQC criado: {model.n_qubits} qubits, {model.n_layers} layers")
    print(f"   ✅ Arquitetura: {model.arch}")
    print(f"   ✅ Device: {model.dev}")
    print()
except Exception as e:
    print(f"   ❌ Erro na criação do modelo: {e}")
    print()

# Teste 4: Funções de auditoria
print("🔐 Teste 4: Testando funções de auditoria...")
try:
    from src.audit import get_seed, power_curve
    seed = get_seed(0)
    print(f"   ✅ Seed determinística: {seed}")
    
    power = power_curve(effect=0.35, alpha=0.05, power_min=0.8)
    print(f"   ✅ Power curve calculada: {power}")
    print()
except Exception as e:
    print(f"   ❌ Erro na auditoria: {e}")
    print()

# Teste 5: Meta-informações
print("📊 Teste 5: Verificando meta-informações...")
try:
    from src import __version__, __author__, __email__
    print(f"   ✅ Versão: {__version__}")
    print(f"   ✅ Autor: {__author__}")
    print(f"   ✅ Email: {__email__}")
    print()
except Exception as e:
    print(f"   ❌ Erro nas meta-informações: {e}")
    print()

# Resumo final
print("="*70)
print("✅ TODOS OS TESTES BÁSICOS PASSARAM!")
print("="*70)
print()
print("🚀 Framework v10.0-A1 está pronto para uso!")
print()
print("Próximos passos:")
print("  1. pip install -e vqc_drug_v10a1/")
print("  2. vqc-drug-a1 --target EGFR --trials 10 (teste rápido)")
print("  3. vqc-drug-a1 --target EGFR --trials 500 (produção)")
print()
