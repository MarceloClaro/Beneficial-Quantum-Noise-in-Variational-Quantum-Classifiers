"""
Tarefa 2: Adaptação do Optuna para QAOA
Otimização bayesiana de hiperparâmetros QAOA para maximizar approximation_ratio.

Referências:
- Optuna Documentation: https://optuna.org/
- Akiba et al. (2019). "Optuna: A Next-generation Hyperparameter Optimization Framework"
"""

from typing import Dict, Any, Optional, Callable
import numpy as np

try:
    import optuna
    from optuna.samplers import TPESampler
    from optuna.pruners import MedianPruner
    OPTUNA_AVAILABLE = True
except ImportError:
    OPTUNA_AVAILABLE = False
    print("⚠️  Optuna não disponível. Instale com: pip install optuna")


def objetivo_optuna_qaoa(
    trial: 'optuna.Trial',
    problem,
    evaluate_fn: Callable,
    config: Dict[str, Any]
) -> float:
    """
    Função objetivo para otimização Optuna no contexto QAOA.
    
    Otimiza:
    - p_layers: Profundidade do circuito (1-10)
    - noise_level: Nível de ruído (0.0-0.05)
    - noise_schedule: Tipo de schedule (constant, linear, exponential)
    - optimizer: Otimizador clássico (COBYLA, SLSQP, Powell)
    
    Args:
        trial: Trial do Optuna
        problem: Instância do problema MaxCut
        evaluate_fn: Função que executa QAOA e retorna approximation_ratio
        config: Configuração base do experimento
    
    Returns:
        approximation_ratio: Métrica a maximizar (0.0-1.0)
    """
    # Sugerir hiperparâmetros
    p_layers = trial.suggest_int('p_layers', 1, 10)
    
    # Tipo de ruído
    noise_model = trial.suggest_categorical(
        'noise_model',
        ['sem_ruido', 'depolarizing', 'amplitude_damping', 'phase_damping', 'thermal']
    )
    
    # Nível de ruído (log-scale para melhor exploração)
    if noise_model != 'sem_ruido':
        noise_level = trial.suggest_float('noise_level', 1e-4, 5e-2, log=True)
    else:
        noise_level = 0.0
    
    # Schedule de ruído
    noise_schedule = trial.suggest_categorical(
        'noise_schedule',
        ['constant', 'linear', 'exponential']
    )
    
    # Otimizador clássico
    optimizer = trial.suggest_categorical(
        'optimizer',
        ['COBYLA', 'SLSQP', 'Powell', 'Nelder-Mead']
    )
    
    # Estratégia de inicialização
    init_strategy = trial.suggest_categorical(
        'init_strategy',
        ['random', 'heuristic', 'zeros']
    )
    
    # Shots (para problemas grandes, menos shots = mais rápido)
    if problem.n_nodes > 50:
        shots = trial.suggest_categorical('shots', [512, 1024])
    else:
        shots = trial.suggest_categorical('shots', [1024, 2048, 4096])
    
    # Criar configuração para este trial
    trial_config = {
        'model': {
            'p_layers': p_layers
        },
        'noise': {
            'enabled': noise_level > 0,
            'model': noise_model,
            'schedule': noise_schedule,
            'params': {'p': [noise_level]}
        },
        'optimization': {
            'optimizer': optimizer,
            'maxiter': config.get('optimization', {}).get('maxiter', 100),
            'initial_params': init_strategy
        },
        'frameworks': {
            'qiskit': {
                'shots': shots,
                'method': config.get('frameworks', {}).get('qiskit', {}).get('method', 'automatic')
            }
        }
    }
    
    # Executar QAOA com esta configuração
    try:
        approximation_ratio = evaluate_fn(problem, trial_config)
        
        # Reportar métricas intermediárias para pruning
        trial.report(approximation_ratio, step=0)
        
        # Verificar se deve fazer pruning
        if trial.should_prune():
            raise optuna.TrialPruned()
        
        return approximation_ratio
    
    except Exception as e:
        # Em caso de erro, retornar valor baixo
        print(f"Erro no trial {trial.number}: {e}")
        return 0.0


def otimizar_hiperparametros_qaoa(
    problem,
    evaluate_fn: Callable,
    config: Dict[str, Any],
    n_trials: int = 100,
    timeout: Optional[int] = None,
    n_jobs: int = 1
) -> Dict[str, Any]:
    """
    Executa otimização bayesiana de hiperparâmetros QAOA.
    
    Args:
        problem: Instância do problema MaxCut
        evaluate_fn: Função que executa QAOA e retorna approximation_ratio
        config: Configuração base do experimento
        n_trials: Número de trials
        timeout: Tempo limite em segundos (opcional)
        n_jobs: Número de jobs paralelos (1 = sequencial)
    
    Returns:
        Dicionário com resultados da otimização
    """
    if not OPTUNA_AVAILABLE:
        raise ImportError("Optuna não disponível. Instale com: pip install optuna")
    
    # Criar estudo
    study = optuna.create_study(
        direction='maximize',  # Maximizar approximation_ratio
        sampler=TPESampler(seed=42),
        pruner=MedianPruner(
            n_startup_trials=10,
            n_warmup_steps=5
        ),
        study_name=f"qaoa_optimization_{problem.n_nodes}nodes"
    )
    
    # Otimizar
    study.optimize(
        lambda trial: objetivo_optuna_qaoa(trial, problem, evaluate_fn, config),
        n_trials=n_trials,
        timeout=timeout,
        n_jobs=n_jobs,
        show_progress_bar=True
    )
    
    # Resultados
    best_params = study.best_params
    best_value = study.best_value
    best_trial = study.best_trial
    
    print(f"\n{'='*80}")
    print("OTIMIZAÇÃO BAYESIANA CONCLUÍDA")
    print(f"{'='*80}")
    print(f"Melhor approximation_ratio: {best_value:.4f}")
    print(f"\nMelhores hiperparâmetros:")
    for key, value in best_params.items():
        print(f"  {key}: {value}")
    print(f"{'='*80}\n")
    
    return {
        'best_params': best_params,
        'best_value': best_value,
        'best_trial': best_trial,
        'study': study,
        'n_trials': len(study.trials),
        'duration': (study.trials[-1].datetime_complete - study.trials[0].datetime_start).total_seconds()
    }


def criar_espaco_busca_customizado(
    min_p_layers: int = 1,
    max_p_layers: int = 5,
    noise_models: list = None,
    noise_range: tuple = (1e-4, 1e-2)
) -> Dict[str, Any]:
    """
    Cria espaço de busca customizado para Optuna.
    
    Args:
        min_p_layers: Profundidade mínima
        max_p_layers: Profundidade máxima
        noise_models: Lista de modelos de ruído a testar
        noise_range: Range de níveis de ruído (min, max)
    
    Returns:
        Dicionário com configuração do espaço de busca
    """
    if noise_models is None:
        noise_models = ['sem_ruido', 'depolarizing', 'phase_damping']
    
    return {
        'p_layers': {'type': 'int', 'range': (min_p_layers, max_p_layers)},
        'noise_model': {'type': 'categorical', 'choices': noise_models},
        'noise_level': {'type': 'float', 'range': noise_range, 'log': True},
        'noise_schedule': {'type': 'categorical', 'choices': ['constant', 'linear']},
        'optimizer': {'type': 'categorical', 'choices': ['COBYLA', 'SLSQP']},
        'init_strategy': {'type': 'categorical', 'choices': ['random', 'heuristic']}
    }


def analisar_importancia_parametros(study: 'optuna.Study') -> Dict[str, float]:
    """
    Analisa importância dos hiperparâmetros.
    
    Args:
        study: Estudo Optuna completo
    
    Returns:
        Dicionário com importâncias normalizadas
    """
    if not OPTUNA_AVAILABLE:
        return {}
    
    try:
        importances = optuna.importance.get_param_importances(study)
        return importances
    except Exception as e:
        print(f"Erro ao calcular importâncias: {e}")
        return {}


def gerar_relatorio_otimizacao(
    study: 'optuna.Study',
    output_file: Optional[str] = None
) -> str:
    """
    Gera relatório detalhado da otimização.
    
    Args:
        study: Estudo Optuna
        output_file: Arquivo de saída (opcional)
    
    Returns:
        String com relatório formatado
    """
    relatorio = []
    relatorio.append("="*80)
    relatorio.append("RELATÓRIO DE OTIMIZAÇÃO BAYESIANA - QAOA")
    relatorio.append("="*80)
    relatorio.append("")
    
    # Estatísticas gerais
    relatorio.append("ESTATÍSTICAS GERAIS")
    relatorio.append("-"*80)
    relatorio.append(f"Trials completos: {len([t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE])}")
    relatorio.append(f"Trials com pruning: {len([t for t in study.trials if t.state == optuna.trial.TrialState.PRUNED])}")
    relatorio.append(f"Trials falhados: {len([t for t in study.trials if t.state == optuna.trial.TrialState.FAIL])}")
    relatorio.append("")
    
    # Melhor resultado
    relatorio.append("MELHOR RESULTADO")
    relatorio.append("-"*80)
    relatorio.append(f"Approximation Ratio: {study.best_value:.4f}")
    relatorio.append(f"Trial number: {study.best_trial.number}")
    relatorio.append("")
    relatorio.append("Hiperparâmetros:")
    for key, value in study.best_params.items():
        relatorio.append(f"  {key}: {value}")
    relatorio.append("")
    
    # Importância dos parâmetros
    importances = analisar_importancia_parametros(study)
    if importances:
        relatorio.append("IMPORTÂNCIA DOS HIPERPARÂMETROS")
        relatorio.append("-"*80)
        for param, importance in sorted(importances.items(), key=lambda x: x[1], reverse=True):
            relatorio.append(f"  {param}: {importance:.4f}")
        relatorio.append("")
    
    # Top 5 trials
    relatorio.append("TOP 5 TRIALS")
    relatorio.append("-"*80)
    trials_completos = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
    trials_ordenados = sorted(trials_completos, key=lambda t: t.value, reverse=True)[:5]
    
    for i, trial in enumerate(trials_ordenados, 1):
        relatorio.append(f"{i}. Trial {trial.number}: {trial.value:.4f}")
        for key, value in trial.params.items():
            relatorio.append(f"     {key}: {value}")
        relatorio.append("")
    
    relatorio.append("="*80)
    
    texto_relatorio = "\n".join(relatorio)
    
    # Salvar em arquivo se especificado
    if output_file:
        with open(output_file, 'w') as f:
            f.write(texto_relatorio)
        print(f"Relatório salvo em: {output_file}")
    
    return texto_relatorio


if __name__ == "__main__":
    print("Testando módulo de otimização Optuna...")
    
    if OPTUNA_AVAILABLE:
        print("✓ Optuna disponível")
        
        # Teste de criação de espaço de busca
        print("\n1. Criando espaço de busca customizado...")
        espaco = criar_espaco_busca_customizado()
        print(f"   Parâmetros: {list(espaco.keys())}")
        
        # Teste de função objetivo (mock)
        print("\n2. Testando função objetivo (mock)...")
        
        def mock_evaluate(problem, config):
            """Mock de avaliação para teste."""
            return np.random.uniform(0.5, 0.95)
        
        # Criar problema mock
        class MockProblem:
            n_nodes = 20
        
        # Executar pequena otimização de teste
        try:
            resultado = otimizar_hiperparametros_qaoa(
                MockProblem(),
                mock_evaluate,
                {},
                n_trials=5
            )
            
            print(f"\n   Melhor valor: {resultado['best_value']:.4f}")
            print(f"   Tempo: {resultado['duration']:.2f}s")
            print(f"   ✓ Otimização de teste bem-sucedida")
        
        except Exception as e:
            print(f"   ⚠️  Erro no teste: {e}")
    
    else:
        print("❌ Optuna não disponível para testes")
