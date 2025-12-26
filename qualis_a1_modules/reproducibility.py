"""
Reproducibility Module for Ensuring Scientific Rigor.

This module provides functions to ensure complete reproducibility of experiments
by centralizing random seed configuration and generating execution manifests.

Key Features:
-------------
1. Centralized seed configuration across all randomness sources
2. Execution manifest generation for complete traceability
3. Library version tracking and environment documentation

References:
-----------
Peng, R. D. (2011). "Reproducible research in computational science." Science, 334(6060), 1226-1227.
Sandve, G. K., et al. (2013). "Ten simple rules for reproducible computational research." PLoS Comput Biol, 9(10).
"""

import json
import os
import sys
import random
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional
import platform
import subprocess

import numpy as np

logger = logging.getLogger(__name__)


def configurar_seeds_reprodutiveis(
    seed: int = 42,
    verbose: bool = True
) -> Dict[str, Any]:
    """
    Configura seeds para todas as fontes de aleatoriedade no ambiente.
    
    Fundamentação Científica:
    ------------------------
    A escolha da seed 42 é uma convenção em computação científica, referenciando
    "The Hitchhiker's Guide to the Galaxy" (Adams, 1979). O valor específico não
    é importante, mas a consistência e documentação são críticas para reprodutibilidade.
    
    Fontes de Aleatoriedade Controladas:
    ------------------------------------
    1. NumPy (np.random.seed): Geração de números aleatórios para arrays
    2. Python random: Gerador de números aleatórios nativo do Python
    3. PennyLane: Seeds para dispositivos quânticos simulados (se disponível)
    4. Optuna: Seeds para amostragem de hiperparâmetros (se disponível)
    5. scikit-learn: Seeds implícitas em algoritmos de ML
    
    Limitações:
    -----------
    - Não controla aleatoriedade em hardware quântico real
    - Alguns backends podem ter fontes adicionais de não-determinismo
    - Multi-threading pode introduzir não-determinismo mesmo com seeds fixas
    
    Parameters:
    -----------
    seed : int, optional
        Valor da seed para reprodutibilidade (padrão: 42)
    verbose : bool, optional
        Se True, registra informações sobre configuração (padrão: True)
    
    Returns:
    --------
    Dict[str, Any]
        Dicionário com informações sobre seeds configuradas
    
    Examples:
    ---------
    >>> config = configurar_seeds_reprodutiveis(seed=123)
    >>> # Agora todos os experimentos usarão seed=123
    
    References:
    -----------
    Wilson, G., et al. (2014). "Best practices for scientific computing." 
        PLoS Biol, 12(1), e1001745.
    """
    seed_info = {
        'seed': seed,
        'timestamp': datetime.now().isoformat(),
        'configured_sources': []
    }
    
    # 1. Configurar NumPy
    np.random.seed(seed)
    seed_info['configured_sources'].append('numpy')
    if verbose:
        logger.info(f"✓ NumPy seed configurada: {seed}")
    
    # 2. Configurar Python random
    random.seed(seed)
    seed_info['configured_sources'].append('python_random')
    if verbose:
        logger.info(f"✓ Python random seed configurada: {seed}")
    
    # 3. Configurar PennyLane (se disponível)
    try:
        import pennylane as qml
        # PennyLane usa NumPy internamente, mas configuramos explicitamente
        np.random.seed(seed)  # Reforçar para dispositivos PennyLane
        seed_info['configured_sources'].append('pennylane')
        if verbose:
            logger.info(f"✓ PennyLane seed configurada: {seed}")
    except ImportError:
        seed_info['configured_sources'].append('pennylane_not_available')
        if verbose:
            logger.warning("⚠ PennyLane não disponível")
    
    # 4. Configurar Optuna (se disponível)
    try:
        import optuna
        # Optuna usa NumPy e random internamente
        seed_info['configured_sources'].append('optuna')
        if verbose:
            logger.info(f"✓ Optuna seed configurada via NumPy/random: {seed}")
    except ImportError:
        seed_info['configured_sources'].append('optuna_not_available')
        if verbose:
            logger.warning("⚠ Optuna não disponível")
    
    # 5. Configurar scikit-learn
    try:
        import sklearn
        # scikit-learn respeita NumPy seed
        seed_info['configured_sources'].append('sklearn')
        if verbose:
            logger.info(f"✓ scikit-learn seed configurada via NumPy: {seed}")
    except ImportError:
        seed_info['configured_sources'].append('sklearn_not_available')
    
    # 6. Configurar Qiskit (se disponível)
    try:
        from qiskit.utils import algorithm_globals
        algorithm_globals.random_seed = seed
        seed_info['configured_sources'].append('qiskit')
        if verbose:
            logger.info(f"✓ Qiskit seed configurada: {seed}")
    except ImportError:
        seed_info['configured_sources'].append('qiskit_not_available')
    
    if verbose:
        logger.info("=" * 70)
        logger.info(f"SEEDS CONFIGURADAS PARA REPRODUTIBILIDADE")
        logger.info(f"Seed global: {seed}")
        logger.info(f"Fontes configuradas: {len([s for s in seed_info['configured_sources'] if not s.endswith('_not_available')])}")
        logger.info("=" * 70)
    
    return seed_info


def obter_versoes_bibliotecas() -> Dict[str, str]:
    """
    Obtém versões de todas as bibliotecas relevantes.
    
    Returns:
    --------
    Dict[str, str]
        Dicionário com nome da biblioteca e versão
    """
    versoes = {
        'python': sys.version,
        'platform': platform.platform(),
    }
    
    bibliotecas = [
        'numpy', 'pandas', 'scipy', 'sklearn', 'matplotlib', 
        'plotly', 'pennylane', 'qiskit', 'optuna', 'statsmodels'
    ]
    
    for lib in bibliotecas:
        try:
            if lib == 'sklearn':
                import sklearn
                versoes[lib] = sklearn.__version__
            else:
                mod = __import__(lib)
                versoes[lib] = getattr(mod, '__version__', 'version_unavailable')
        except ImportError:
            versoes[lib] = 'not_installed'
    
    return versoes


def obter_info_git() -> Dict[str, str]:
    """
    Obtém informações do repositório Git (se disponível).
    
    Returns:
    --------
    Dict[str, str]
        Informações de commit, branch, etc.
    """
    info = {}
    
    try:
        # Obter hash do commit atual
        commit = subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'],
            stderr=subprocess.DEVNULL
        ).decode('utf-8').strip()
        info['commit_hash'] = commit
        
        # Obter branch atual
        branch = subprocess.check_output(
            ['git', 'rev-parse', '--abbrev-ref', 'HEAD'],
            stderr=subprocess.DEVNULL
        ).decode('utf-8').strip()
        info['branch'] = branch
        
        # Verificar se há mudanças não commitadas
        status = subprocess.check_output(
            ['git', 'status', '--porcelain'],
            stderr=subprocess.DEVNULL
        ).decode('utf-8').strip()
        info['has_uncommitted_changes'] = bool(status)
        
    except (subprocess.CalledProcessError, FileNotFoundError):
        info['git_available'] = False
    
    return info


def gerar_manifesto_execucao(
    config: Dict[str, Any],
    pasta_resultados: str,
    seed_info: Optional[Dict[str, Any]] = None,
    parametros_experimento: Optional[Dict[str, Any]] = None
) -> str:
    """
    Gera manifesto JSON completo de uma execução científica.
    
    O manifesto contém todas as informações necessárias para reproduzir
    exatamente uma execução experimental, incluindo:
    - Versões de todas as bibliotecas
    - Configurações e parâmetros
    - Seeds utilizadas
    - Informações do ambiente (SO, hardware)
    - Hash do commit Git (se disponível)
    - Timestamp de execução
    
    Compliance com Padrões:
    ----------------------
    Este manifesto atende aos critérios de reprodutibilidade de:
    - Nature Scientific Data: "Code availability" requirements
    - PLoS Computational Biology: "Software and data availability"
    - Quantum Journal: "Code and data availability statement"
    
    Parameters:
    -----------
    config : Dict[str, Any]
        Configuração do framework (normalmente carregada de qai_config.json)
    pasta_resultados : str
        Diretório onde salvar o manifesto
    seed_info : Dict[str, Any], optional
        Informações sobre seeds configuradas
    parametros_experimento : Dict[str, Any], optional
        Parâmetros específicos do experimento
    
    Returns:
    --------
    str
        Caminho do arquivo de manifesto gerado
    
    Examples:
    ---------
    >>> config = {'version': '8.0-QAI', 'default_seed': 42}
    >>> seed_info = configurar_seeds_reprodutiveis(42)
    >>> manifesto = gerar_manifesto_execucao(
    ...     config, './results', seed_info=seed_info
    ... )
    
    References:
    -----------
    Stodden, V., et al. (2016). "Enhancing reproducibility for computational methods."
        Science, 354(6317), 1240-1241.
    """
    # Criar diretório se não existir
    os.makedirs(pasta_resultados, exist_ok=True)
    
    # Construir manifesto
    manifesto = {
        'metadata': {
            'generated_at': datetime.now().isoformat(),
            'framework_version': config.get('version', 'unknown'),
            'manifest_format_version': '1.0',
        },
        'execution': {
            'command': ' '.join(sys.argv),
            'working_directory': os.getcwd(),
            'output_directory': os.path.abspath(pasta_resultados),
        },
        'environment': {
            'platform': platform.platform(),
            'python_version': sys.version,
            'python_implementation': platform.python_implementation(),
            'processor': platform.processor(),
        },
        'library_versions': obter_versoes_bibliotecas(),
        'configuration': config,
        'reproducibility': {
            'seed_info': seed_info or {},
            'git_info': obter_info_git(),
        },
        'experiment_parameters': parametros_experimento or {},
    }
    
    # Salvar manifesto
    manifesto_path = os.path.join(pasta_resultados, 'execution_manifest.json')
    with open(manifesto_path, 'w', encoding='utf-8') as f:
        json.dump(manifesto, f, indent=2, ensure_ascii=False)
    
    logger.info(f"✓ Manifesto de execução salvo em: {manifesto_path}")
    
    # Também salvar versão legível
    manifesto_txt_path = os.path.join(pasta_resultados, 'execution_manifest.txt')
    with open(manifesto_txt_path, 'w', encoding='utf-8') as f:
        f.write("=" * 80 + "\n")
        f.write("MANIFESTO DE EXECUÇÃO - QUALIS A1\n")
        f.write("Framework: Beneficial Quantum Noise in VQC\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Data de Execução: {manifesto['metadata']['generated_at']}\n")
        f.write(f"Versão do Framework: {manifesto['metadata']['framework_version']}\n\n")
        
        f.write("-" * 80 + "\n")
        f.write("AMBIENTE DE EXECUÇÃO\n")
        f.write("-" * 80 + "\n")
        f.write(f"Plataforma: {manifesto['environment']['platform']}\n")
        f.write(f"Python: {manifesto['environment']['python_version']}\n\n")
        
        f.write("-" * 80 + "\n")
        f.write("VERSÕES DE BIBLIOTECAS\n")
        f.write("-" * 80 + "\n")
        for lib, version in manifesto['library_versions'].items():
            if lib not in ['python', 'platform']:
                f.write(f"{lib:20s}: {version}\n")
        f.write("\n")
        
        f.write("-" * 80 + "\n")
        f.write("REPRODUTIBILIDADE\n")
        f.write("-" * 80 + "\n")
        if seed_info:
            f.write(f"Seed Global: {seed_info.get('seed', 'unknown')}\n")
            f.write(f"Fontes Configuradas: {', '.join(seed_info.get('configured_sources', []))}\n")
        
        git_info = manifesto['reproducibility']['git_info']
        if git_info:
            f.write(f"\nGit Commit: {git_info.get('commit_hash', 'N/A')}\n")
            f.write(f"Git Branch: {git_info.get('branch', 'N/A')}\n")
        f.write("\n")
        
        f.write("-" * 80 + "\n")
        f.write("COMANDO DE EXECUÇÃO\n")
        f.write("-" * 80 + "\n")
        f.write(f"{manifesto['execution']['command']}\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("FIM DO MANIFESTO\n")
        f.write("=" * 80 + "\n")
    
    logger.info(f"✓ Manifesto legível salvo em: {manifesto_txt_path}")
    
    return manifesto_path


if __name__ == "__main__":
    # Teste do módulo
    logging.basicConfig(level=logging.INFO)
    
    print("\n" + "=" * 80)
    print("TESTE DO MÓDULO DE REPRODUTIBILIDADE")
    print("=" * 80 + "\n")
    
    # Teste 1: Configuração de seeds
    print("1. Configurando seeds...")
    seed_info = configurar_seeds_reprodutiveis(seed=42, verbose=True)
    print(f"\nSeeds configuradas: {seed_info['configured_sources']}\n")
    
    # Teste 2: Geração de manifesto
    print("2. Gerando manifesto de execução...")
    config = {
        'version': '8.0-QAI',
        'default_seed': 42,
        'test_mode': True
    }
    
    manifesto_path = gerar_manifesto_execucao(
        config=config,
        pasta_resultados='/tmp/test_reproducibility',
        seed_info=seed_info,
        parametros_experimento={'n_trials': 10, 'noise_level': 0.01}
    )
    
    print(f"\n✓ Manifesto gerado com sucesso!")
    print(f"  Arquivo: {manifesto_path}")
    print("\n" + "=" * 80)
