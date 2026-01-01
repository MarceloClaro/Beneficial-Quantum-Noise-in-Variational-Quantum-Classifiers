"""
Execução Rápida - VQC v10.0-A1 com HIV Dataset
Teste com 10 trials para validação rápida (~5 minutos)
"""
import sys
import os
from pathlib import Path
import time

# Adicionar vqc_drug_v10a1 ao path
base_dir = Path(__file__).parent
vqc_dir = base_dir / "vqc_drug_v10a1"
sys.path.insert(0, str(vqc_dir))

print("="*70)
print("🚀 EXECUÇÃO RÁPIDA - VQC v10.0-A1 com HIV")
print("="*70)
print()
print("📋 Configuração:")
print("   Dataset: HIV")
print("   Max qubits: 8 (reduzido para velocidade)")
print("   Trials: 10 (teste rápido)")
print("   Tempo estimado: ~3-5 minutos")
print()

# Importar módulos
try:
    from src.data import load_split
    from src.tune import run_study
    from src.audit import get_seed
    print("✅ Módulos importados com sucesso")
    print()
except ImportError as e:
    print(f"❌ Erro ao importar módulos: {e}")
    print()
    print("💡 Instale as dependências:")
    print("   pip install pennylane torch optuna rdkit scikit-learn pandas numpy")
    sys.exit(1)

# Carregar dados
print("📦 Carregando dataset HIV...")
start_time = time.time()

try:
    X_train, X_test, y_train, y_test = load_split("HIV", n_qubits=8)
    load_time = time.time() - start_time
    
    print(f"✅ Dataset carregado em {load_time:.1f}s")
    print(f"   Train: {len(y_train)} amostras")
    print(f"   Test: {len(y_test)} amostras")
    print(f"   Ativos: {y_train.mean():.1%}")
    print(f"   Features: {X_train.shape[1]} (PCA reduzido)")
    print()
except Exception as e:
    print(f"❌ Erro ao carregar dados: {e}")
    print()
    print("💡 Verifique conexão com internet (download automático)")
    sys.exit(1)

# Executar estudo com Optuna
print("🔍 Executando otimização com Optuna (10 trials)...")
print("   (Isso pode levar 3-5 minutos...)")
print()

try:
    study_start = time.time()
    study = run_study(X_train, y_train, target="HIV", n_trials=10, max_qubits=8)
    study_time = time.time() - study_start
    
    print()
    print("="*70)
    print("✅ OTIMIZAÇÃO CONCLUÍDA!")
    print("="*70)
    print()
    
    # Resultados
    print("📊 RESULTADOS:")
    print(f"   Tempo total: {study_time:.1f}s ({study_time/60:.1f} min)")
    print(f"   Melhor AUC: {study.best_value:.4f}")
    print()
    
    print("🎯 Melhores Hiperparâmetros:")
    for key, value in study.best_params.items():
        print(f"   {key}: {value}")
    print()
    
    # Top 3 trials
    df = study.trials_dataframe()
    df_sorted = df.sort_values('value', ascending=False).head(3)
    
    print("🏆 Top 3 Configurações:")
    for idx, (i, row) in enumerate(df_sorted.iterrows(), 1):
        print(f"   #{idx}: AUC={row['value']:.4f}")
        print(f"       qubits={row.get('params_n_qubits', 'N/A')}, "
              f"layers={row.get('params_n_layers', 'N/A')}, "
              f"noise={row.get('params_noise_type', 'N/A')}")
    print()
    
    # Salvar resultados
    output_dir = base_dir / "results_HIV_quick_test"
    output_dir.mkdir(exist_ok=True)
    
    # CSV com todos os trials
    csv_path = output_dir / "trials_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"💾 Resultados salvos em: {csv_path}")
    print()
    
    # Relatório resumido
    report_path = output_dir / "quick_report.txt"
    with open(report_path, "w") as f:
        f.write("="*70 + "\n")
        f.write("VQC-Molecular v10.0-A1 - Quick Test Report\n")
        f.write("="*70 + "\n\n")
        f.write(f"Dataset: HIV\n")
        f.write(f"Trials: 10\n")
        f.write(f"Qubits: 8\n")
        f.write(f"Tempo: {study_time:.1f}s\n\n")
        f.write(f"Melhor AUC: {study.best_value:.4f}\n\n")
        f.write("Melhores Parâmetros:\n")
        for key, value in study.best_params.items():
            f.write(f"  {key}: {value}\n")
    
    print(f"📄 Relatório salvo em: {report_path}")
    print()
    
except Exception as e:
    print(f"❌ Erro durante otimização: {e}")
    print()
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Resumo final
print("="*70)
print("🎉 TESTE RÁPIDO CONCLUÍDO COM SUCESSO!")
print("="*70)
print()
print("📈 Próximos passos:")
print("   1. Aumentar trials (50-500) para melhores resultados")
print("   2. Aumentar qubits (12-20) para maior capacidade")
print("   3. Executar pipeline completo: vqc-drug-a1 --target HIV --trials 500")
print()
print(f"⏱️  Tempo total de execução: {time.time() - start_time:.1f}s")
print()
