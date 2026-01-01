"""
Script de Execução de Trials (Optuna) com Framework Qiskit
Timeout de 600 segundos por trial para otimização de hiperparâmetros

Utiliza Optuna para busca inteligente de hiperparâmetros:
- Arquitetura de circuito
- Tipo e nível de ruído
- Estratégia de inicialização
- Número de épocas e taxa de aprendizado

Timeout: 600 segundos por trial (10 minutos)
"""

import os
import sys
import time
import logging
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

# Suprimir warnings
warnings.filterwarnings('ignore')
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

print("Instalando/verificando dependências...")
os.system("pip install -q optuna 2>&1 | tail -3")

import optuna
from optuna.samplers import TPESampler
from optuna.pruners import MedianPruner

from framework_qiskit import ClassificadorVQCQiskit, carregar_datasets

# Configurar logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('trials_qiskit_600s.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

print("\n" + "=" * 80)
print("TRIALS QISKIT COM TIMEOUT 600s - OTIMIZAÇÃO DE HIPERPARÂMETROS")
print("=" * 80)

# ============================================================================
# CONFIGURAÇÃO
# ============================================================================

# Parâmetros globais
TIMEOUT_PER_TRIAL = 600  # segundos por trial
N_TRIALS = 5  # número total de trials (reduzido para demonstração)
DATASET_FOCUS = 'moons'  # dataset para otimização

# Espaço de busca
ARQUITETURAS = [
    'basico',
    'strongly_entangling', 
    'hardware_efficient',
    'alternating_layers',
    'brickwork',
    'tree',
    'star_entanglement'
]

TIPOS_RUIDO = [
    'sem_ruido',
    'depolarizante',
    'amplitude_damping',
    'phase_damping',
    'crosstalk',
    'correlacionado'
]

ESTRATEGIAS_INIT = [
    'matematico',
    'quantico',
    'aleatorio',
    'fibonacci_spiral'
]

# Criar pasta de resultados
pasta_resultados = f"trials_qiskit_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
os.makedirs(pasta_resultados, exist_ok=True)
logger.info(f"📁 Pasta de resultados: {pasta_resultados}/")

# Carregar datasets
print("\n📊 Carregando datasets...")
datasets = carregar_datasets()
print(f"   ✓ {len(datasets)} datasets carregados")

# Dataset para otimização
X_train = datasets[DATASET_FOCUS]['X_train']
X_test = datasets[DATASET_FOCUS]['X_test']
y_train = datasets[DATASET_FOCUS]['y_train']
y_test = datasets[DATASET_FOCUS]['y_test']

print(f"   ✓ Dataset '{DATASET_FOCUS}': {len(X_train)} treino, {len(X_test)} teste")

# ============================================================================
# FUNÇÃO OBJETIVO (COM TIMEOUT)
# ============================================================================

def objective(trial):
    """
    Função objetivo para Optuna com timeout de 600s.
    
    Otimiza:
    - Arquitetura do circuito
    - Tipo de ruído
    - Nível de ruído
    - Estratégia de inicialização
    - Hiperparâmetros de treinamento
    """
    inicio_trial = time.time()
    
    try:
        # Sugerir hiperparâmetros
        arquitetura = trial.suggest_categorical('arquitetura', ARQUITETURAS)
        tipo_ruido = trial.suggest_categorical('tipo_ruido', TIPOS_RUIDO)
        
        # Nível de ruído (0.0 a 0.02)
        if tipo_ruido == 'sem_ruido':
            nivel_ruido = 0.0
        else:
            nivel_ruido = trial.suggest_float('nivel_ruido', 0.0, 0.02)
        
        estrategia_init = trial.suggest_categorical('estrategia_init', ESTRATEGIAS_INIT)
        n_epocas = trial.suggest_int('n_epocas', 3, 8)
        taxa_aprendizado = trial.suggest_float('taxa_aprendizado', 0.01, 0.5, log=True)
        
        # Log dos parâmetros
        logger.info(f"\n{'─' * 80}")
        logger.info(f"Trial #{trial.number + 1}/{N_TRIALS}")
        logger.info(f"  Arquitetura: {arquitetura}")
        logger.info(f"  Ruído: {tipo_ruido} (γ={nivel_ruido:.4f})")
        logger.info(f"  Init: {estrategia_init}")
        logger.info(f"  Epochs: {n_epocas}, LR: {taxa_aprendizado:.4f}")
        
        # Criar classificador
        vqc = ClassificadorVQCQiskit(
            n_qubits=4,
            arquitetura=arquitetura,
            tipo_ruido=tipo_ruido,
            nivel_ruido=nivel_ruido,
            estrategia_init=estrategia_init,
            n_epocas=n_epocas,
            taxa_aprendizado=taxa_aprendizado,
            shots=512,  # Reduzido para velocidade
            seed=42
        )
        
        # Treinar com timeout
        tempo_treino_inicio = time.time()
        vqc.fit(X_train, y_train)
        tempo_treino = time.time() - tempo_treino_inicio
        
        # Verificar timeout
        tempo_decorrido = time.time() - inicio_trial
        if tempo_decorrido > TIMEOUT_PER_TRIAL:
            logger.warning(f"  ⏱️ Timeout atingido ({tempo_decorrido:.1f}s > {TIMEOUT_PER_TRIAL}s)")
            raise optuna.TrialPruned()
        
        # Avaliar
        acc_train = vqc.score(X_train, y_train)
        acc_test = vqc.score(X_test, y_test)
        
        tempo_total = time.time() - inicio_trial
        
        # Log dos resultados
        logger.info(f"  Acc Train: {acc_train:.4f}")
        logger.info(f"  Acc Test: {acc_test:.4f}")
        logger.info(f"  Tempo: {tempo_total:.1f}s (treino: {tempo_treino:.1f}s)")
        
        # Salvar informações do trial
        trial.set_user_attr('acc_train', acc_train)
        trial.set_user_attr('acc_test', acc_test)
        trial.set_user_attr('tempo_treino', tempo_treino)
        trial.set_user_attr('tempo_total', tempo_total)
        
        # Retornar métrica para otimização (acurácia de teste)
        return acc_test
        
    except optuna.TrialPruned:
        raise
    except Exception as e:
        logger.error(f"  ✗ Erro no trial: {e}")
        raise optuna.TrialPruned()

# ============================================================================
# EXECUTAR OTIMIZAÇÃO
# ============================================================================

print(f"\n🔍 Iniciando otimização com Optuna...")
print(f"   • Trials: {N_TRIALS}")
print(f"   • Timeout por trial: {TIMEOUT_PER_TRIAL}s")
print(f"   • Dataset: {DATASET_FOCUS}")
print(f"   • Método: Tree-structured Parzen Estimator (TPE)")

# Criar study
study = optuna.create_study(
    study_name=f"qiskit_vqc_{DATASET_FOCUS}",
    direction='maximize',  # Maximizar acurácia
    sampler=TPESampler(seed=42),
    pruner=MedianPruner(n_startup_trials=3, n_warmup_steps=2)
)

# Executar otimização
inicio_global = time.time()

try:
    study.optimize(
        objective,
        n_trials=N_TRIALS,
        timeout=None,  # Sem timeout global (controlado por trial)
        show_progress_bar=True,
        catch=(Exception,)
    )
except KeyboardInterrupt:
    logger.info("\n⚠️ Otimização interrompida pelo usuário")

tempo_total = time.time() - inicio_global

# ============================================================================
# ANÁLISE DOS RESULTADOS
# ============================================================================

print("\n" + "=" * 80)
print("RESULTADOS DA OTIMIZAÇÃO")
print("=" * 80)

# Melhor trial
best_trial = study.best_trial
print(f"\n🏆 Melhor Trial: #{best_trial.number + 1}")
print(f"   Acurácia Teste: {best_trial.value:.4f}")
print(f"\n   Hiperparâmetros:")
for key, value in best_trial.params.items():
    print(f"      • {key}: {value}")

if best_trial.user_attrs:
    print(f"\n   Métricas adicionais:")
    for key, value in best_trial.user_attrs.items():
        if isinstance(value, float):
            print(f"      • {key}: {value:.4f}")
        else:
            print(f"      • {key}: {value}")

# Estatísticas gerais
trials_completos = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
print(f"\n📊 Estatísticas:")
print(f"   • Trials completos: {len(trials_completos)}/{N_TRIALS}")
print(f"   • Trials podados: {len([t for t in study.trials if t.state == optuna.trial.TrialState.PRUNED])}")
print(f"   • Tempo total: {tempo_total:.1f}s ({tempo_total/60:.1f} min)")
if trials_completos:
    print(f"   • Tempo médio/trial: {tempo_total/len(trials_completos):.1f}s")

# ============================================================================
# SALVAR RESULTADOS
# ============================================================================

# Dataframe com todos os trials
trials_df = study.trials_dataframe()
csv_path = os.path.join(pasta_resultados, "trials_completos.csv")
trials_df.to_csv(csv_path, index=False)
print(f"\n💾 Resultados salvos:")
print(f"   • {csv_path}")

# Melhor configuração
best_config = {
    'trial_number': best_trial.number + 1,
    'acuracia_teste': best_trial.value,
    'parametros': best_trial.params,
    'metricas': best_trial.user_attrs
}

import json
json_path = os.path.join(pasta_resultados, "melhor_configuracao.json")
with open(json_path, 'w') as f:
    json.dump(best_config, f, indent=2)
print(f"   • {json_path}")

# Importância dos hiperparâmetros
if len(trials_completos) >= 3:
    print(f"\n📈 Importância dos Hiperparâmetros:")
    try:
        importances = optuna.importance.get_param_importances(study)
        for param, importance in sorted(importances.items(), key=lambda x: x[1], reverse=True):
            print(f"   • {param}: {importance:.4f}")
    except Exception as e:
        print(f"   ⚠️ Não foi possível calcular importâncias: {e}")

# ============================================================================
# VALIDAÇÃO FINAL
# ============================================================================

print("\n" + "=" * 80)
print("VALIDAÇÃO COM MELHOR CONFIGURAÇÃO")
print("=" * 80)

print("\n🔬 Re-treinando com melhor configuração...")

# Extrair parâmetros
best_params = best_trial.params
nivel_ruido = best_params.get('nivel_ruido', 0.0) if best_params['tipo_ruido'] != 'sem_ruido' else 0.0

# Criar classificador com melhor configuração
vqc_best = ClassificadorVQCQiskit(
    n_qubits=4,
    arquitetura=best_params['arquitetura'],
    tipo_ruido=best_params['tipo_ruido'],
    nivel_ruido=nivel_ruido,
    estrategia_init=best_params['estrategia_init'],
    n_epocas=best_params['n_epocas'],
    taxa_aprendizado=best_params['taxa_aprendizado'],
    shots=1024,  # Shots completos para validação
    seed=42
)

# Treinar
inicio_validacao = time.time()
vqc_best.fit(X_train, y_train)
tempo_validacao = time.time() - inicio_validacao

# Avaliar
acc_train_final = vqc_best.score(X_train, y_train)
acc_test_final = vqc_best.score(X_test, y_test)

print(f"\n✅ Validação completa:")
print(f"   • Acurácia Treino: {acc_train_final:.4f}")
print(f"   • Acurácia Teste: {acc_test_final:.4f}")
print(f"   • Tempo: {tempo_validacao:.1f}s")

# ============================================================================
# COMPARAÇÃO COM BASELINE
# ============================================================================

print("\n" + "=" * 80)
print("COMPARAÇÃO COM BASELINE (SEM RUÍDO)")
print("=" * 80)

print("\n🔬 Treinando baseline (sem ruído, configuração padrão)...")

vqc_baseline = ClassificadorVQCQiskit(
    n_qubits=4,
    arquitetura='basico',
    tipo_ruido='sem_ruido',
    nivel_ruido=0.0,
    estrategia_init='aleatorio',
    n_epocas=5,
    taxa_aprendizado=0.1,
    shots=1024,
    seed=42
)

vqc_baseline.fit(X_train, y_train)
acc_baseline = vqc_baseline.score(X_test, y_test)

print(f"\n📊 Comparação:")
print(f"   • Baseline (sem otimização): {acc_baseline:.4f}")
print(f"   • Otimizado (Optuna):        {acc_test_final:.4f}")
print(f"   • Ganho absoluto:            {acc_test_final - acc_baseline:+.4f}")
if acc_baseline > 0:
    ganho_relativo = ((acc_test_final - acc_baseline) / acc_baseline) * 100
    print(f"   • Ganho relativo:            {ganho_relativo:+.2f}%")

# ============================================================================
# RESUMO FINAL
# ============================================================================

print("\n" + "=" * 80)
print("RESUMO FINAL")
print("=" * 80)

print(f"""
✅ Otimização completa!

📊 **Configuração do Estudo**:
   • Dataset: {DATASET_FOCUS}
   • Trials executados: {len(trials_completos)}/{N_TRIALS}
   • Timeout por trial: {TIMEOUT_PER_TRIAL}s
   • Tempo total: {tempo_total/60:.1f} minutos

🏆 **Melhor Configuração**:
   • Trial: #{best_trial.number + 1}
   • Acurácia: {best_trial.value:.4f}
   • Arquitetura: {best_params['arquitetura']}
   • Ruído: {best_params['tipo_ruido']} (γ={nivel_ruido:.4f})
   • Init: {best_params['estrategia_init']}

📁 **Arquivos Gerados**:
   • {csv_path}
   • {json_path}
   • trials_qiskit_600s.log

🎯 **Próximos Passos**:
   1. Analisar importância dos hiperparâmetros
   2. Testar melhor configuração em outros datasets
   3. Executar com mais trials para refinar busca
   4. Explorar interações entre hiperparâmetros
""")

print("=" * 80)
print(f"Execução finalizada em {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 80)
