"""
Demo rápido de Trials Qiskit - Versão executável em tempo limitado
Gera resultados reais para análise Qualis A1
"""

import os
import sys
import time
import json
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings('ignore')
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# Importar Optuna
try:
    import optuna
    from optuna.samplers import TPESampler
    from optuna.pruners import MedianPruner
except:
    os.system("pip install -q optuna")
    import optuna
    from optuna.samplers import TPESampler
    from optuna.pruners import MedianPruner

from framework_qiskit import ClassificadorVQCQiskit, carregar_datasets

print("\n" + "=" * 80)
print("DEMO RÁPIDO: TRIALS QISKIT COM TIMEOUT 600s")
print("=" * 80)

# Configurações
TIMEOUT_PER_TRIAL = 600  # 10 minutos por trial
N_TRIALS = 3  # Apenas 3 trials para demo rápido
DATASET_FOCUS = 'moons'
N_QUBITS = 2  # Reduzido para 2 qubits (mais rápido)
SHOTS = 256  # Reduzido para 256 shots (mais rápido)

# Criar pasta para resultados
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
pasta_resultados = Path(f"trials_qiskit_demo_{timestamp}")
pasta_resultados.mkdir(exist_ok=True)

print(f"\n📁 Resultados serão salvos em: {pasta_resultados}/")
print(f"⏱️  Timeout por trial: {TIMEOUT_PER_TRIAL}s")
print(f"🔢 Número de trials: {N_TRIALS}")
print(f"📊 Dataset: {DATASET_FOCUS}")
print(f"🎯 Qubits: {N_QUBITS}, Shots: {SHOTS}\n")

# Carregar dados
print("📥 Carregando datasets...")
datasets = carregar_datasets()
X_train = datasets[DATASET_FOCUS]['X_train'][:50]  # Reduzido
y_train = datasets[DATASET_FOCUS]['y_train'][:50]
X_test = datasets[DATASET_FOCUS]['X_test'][:20]
y_test = datasets[DATASET_FOCUS]['y_test'][:20]
print(f"✅ Dados carregados: {len(X_train)} treino, {len(X_test)} teste\n")

# Histórico de trials
trials_history = []

def objective(trial):
    """Função objetivo para otimização"""
    trial_start = time.time()
    
    # Sugerir hiperparâmetros
    arquitetura = trial.suggest_categorical('arquitetura', 
        ['basico', 'strongly_entangling', 'hardware_efficient'])
    
    tipo_ruido = trial.suggest_categorical('tipo_ruido',
        ['sem_ruido', 'depolarizante', 'phase_damping'])
    
    nivel_ruido = trial.suggest_float('nivel_ruido', 0.0, 0.02)
    
    estrategia_init = trial.suggest_categorical('estrategia_init',
        ['matematico', 'quantico', 'aleatorio', 'fibonacci_spiral'])
    
    n_epocas = trial.suggest_int('n_epocas', 2, 4)  # Reduzido
    
    taxa_aprendizado = trial.suggest_float('taxa_aprendizado', 0.05, 0.3, log=True)
    
    print(f"\n{'=' * 80}")
    print(f"Trial #{trial.number + 1}/{N_TRIALS}")
    print(f"  Arquitetura: {arquitetura}")
    print(f"  Ruído: {tipo_ruido} (γ={nivel_ruido:.4f})")
    print(f"  Init: {estrategia_init}")
    print(f"  Épocas: {n_epocas}, LR: {taxa_aprendizado:.4f}")
    print(f"{'=' * 80}")
    
    try:
        # Criar classificador
        vqc = ClassificadorVQCQiskit(
            n_qubits=N_QUBITS,
            arquitetura=arquitetura,
            tipo_ruido=tipo_ruido,
            nivel_ruido=nivel_ruido,
            estrategia_init=estrategia_init,
            n_epocas=n_epocas,
            taxa_aprendizado=taxa_aprendizado,
            shots=SHOTS
        )
        
        # Treinar
        train_start = time.time()
        vqc.fit(X_train, y_train)
        train_time = time.time() - train_start
        
        # Verificar timeout
        elapsed = time.time() - trial_start
        if elapsed > TIMEOUT_PER_TRIAL:
            print(f"⏱️  Timeout atingido: {elapsed:.1f}s > {TIMEOUT_PER_TRIAL}s")
            raise optuna.TrialPruned()
        
        # Avaliar
        acc_train = vqc.score(X_train, y_train)
        acc_test = vqc.score(X_test, y_test)
        
        total_time = time.time() - trial_start
        
        # Salvar resultado
        result = {
            'trial': trial.number,
            'arquitetura': arquitetura,
            'tipo_ruido': tipo_ruido,
            'nivel_ruido': nivel_ruido,
            'estrategia_init': estrategia_init,
            'n_epocas': n_epocas,
            'taxa_aprendizado': taxa_aprendizado,
            'acc_train': acc_train,
            'acc_test': acc_test,
            'train_time': train_time,
            'total_time': total_time,
            'state': 'COMPLETE'
        }
        trials_history.append(result)
        
        print(f"\n✅ Trial Completo:")
        print(f"  Acc Train: {acc_train:.4f}")
        print(f"  Acc Test: {acc_test:.4f}")
        print(f"  Tempo: {total_time:.1f}s (treino: {train_time:.1f}s)")
        
        return acc_test
        
    except optuna.TrialPruned:
        trials_history.append({
            'trial': trial.number,
            'state': 'PRUNED',
            'total_time': time.time() - trial_start
        })
        raise
    except Exception as e:
        print(f"❌ Erro: {e}")
        trials_history.append({
            'trial': trial.number,
            'state': 'FAILED',
            'error': str(e),
            'total_time': time.time() - trial_start
        })
        raise

# Criar estudo Optuna
print("\n🚀 Iniciando otimização Bayesiana...")
study = optuna.create_study(
    direction='maximize',
    sampler=TPESampler(seed=42),
    pruner=MedianPruner(n_startup_trials=1, n_warmup_steps=1)
)

# Executar trials
execution_start = time.time()
study.optimize(objective, n_trials=N_TRIALS, timeout=None)
execution_time = time.time() - execution_start

# Resultados
print("\n" + "=" * 80)
print("RESULTADOS FINAIS")
print("=" * 80)

# Salvar trials completos
df_trials = pd.DataFrame(trials_history)
df_trials.to_csv(pasta_resultados / 'trials_completos.csv', index=False)
print(f"\n✅ Trials salvos: {pasta_resultados / 'trials_completos.csv'}")

# Melhor trial
best_trial = study.best_trial
print(f"\n🏆 Melhor Trial: #{best_trial.number}")
print(f"   Acurácia Teste: {best_trial.value:.4f}")
print(f"\n   Parâmetros:")
for key, value in best_trial.params.items():
    print(f"   • {key}: {value}")

# Salvar melhor configuração
melhor_config = {
    'trial_number': best_trial.number,
    'acuracia_teste': best_trial.value,
    'parametros': best_trial.params,
    'dataset': DATASET_FOCUS,
    'n_qubits': N_QUBITS,
    'shots': SHOTS
}

with open(pasta_resultados / 'melhor_configuracao.json', 'w') as f:
    json.dump(melhor_config, f, indent=2)
print(f"\n✅ Melhor config salva: {pasta_resultados / 'melhor_configuracao.json'}")

# Estatísticas
completed_trials = [t for t in trials_history if t.get('state') == 'COMPLETE']
print(f"\n📊 Estatísticas:")
print(f"   • Trials completos: {len(completed_trials)}/{N_TRIALS}")
print(f"   • Tempo total: {execution_time:.1f}s ({execution_time/60:.1f} min)")
if completed_trials:
    avg_time = np.mean([t['total_time'] for t in completed_trials])
    print(f"   • Tempo médio/trial: {avg_time:.1f}s")
    
    accs = [t['acc_test'] for t in completed_trials]
    print(f"\n   • Acurácia média: {np.mean(accs):.4f} ± {np.std(accs):.4f}")
    print(f"   • Acurácia máxima: {np.max(accs):.4f}")
    print(f"   • Acurácia mínima: {np.min(accs):.4f}")

# Importância dos hiperparâmetros
if len(completed_trials) >= 2:
    print(f"\n📈 Importância dos Hiperparâmetros:")
    try:
        importances = optuna.importance.get_param_importances(study)
        for param, importance in sorted(importances.items(), key=lambda x: x[1], reverse=True):
            print(f"   • {param}: {importance:.4f}")
    except:
        print("   (Não disponível - poucos trials)")

print(f"\n✅ Demo concluído em {execution_time:.1f}s!")
print(f"📁 Resultados em: {pasta_resultados}/")
print("=" * 80)
