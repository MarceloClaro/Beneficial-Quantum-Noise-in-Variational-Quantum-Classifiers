#!/usr/bin/env python3
# =============================================================================
# EXECUTOR FRAMEWORK QUANTUM ADVANCED V8
# =============================================================================
"""
Script para executar experimentos com Framework Quantum Advanced V8.

Uso:
    python run_framework_quantum_advanced_v8.py
    python run_framework_quantum_advanced_v8.py --dataset hiv --n_qubits 6 --n_layers 3
    python run_framework_quantum_advanced_v8.py --framework qiskit --noise_level 0.05
"""

import argparse
import json
import logging
from pathlib import Path

from framework_quantum_advanced_v8 import (
    ExperimentConfig,
    QuantumCircuitConfig,
    NoiseConfig,
    OptimizationConfig,
    ErrorMitigationConfig,
    QuantumExperimentRunner,
    FrameworkType,
    NoiseModel,
    OptimizationMethod,
    ErrorMitigationTechnique,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def parse_arguments():
    """Parse argumentos de linha de comando."""
    parser = argparse.ArgumentParser(
        description="Framework Quantum Advanced V8 - Executor de Experimentos",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Framework
    parser.add_argument(
        "--framework",
        type=str,
        default="pennylane",
        choices=["pennylane", "qiskit", "cirq"],
        help="Framework quântico a usar"
    )

    # Configuração de circuito
    parser.add_argument(
        "--n_qubits",
        type=int,
        default=4,
        help="Número de qubits"
    )
    parser.add_argument(
        "--n_layers",
        type=int,
        default=2,
        help="Número de camadas do circuito"
    )
    parser.add_argument(
        "--entanglement",
        type=str,
        default="full",
        choices=["full", "linear"],
        help="Tipo de emaranhamento"
    )

    # Dataset
    parser.add_argument(
        "--dataset",
        type=str,
        default="iris",
        choices=["iris", "wine", "breast_cancer", "hiv", "malaria", "tb"],
        help="Dataset a usar"
    )

    # Ruído
    parser.add_argument(
        "--noise_model",
        type=str,
        default="depolarizing",
        choices=["depolarizing", "amplitude_damping", "phase_damping", "pauli", "none"],
        help="Modelo de ruído"
    )
    parser.add_argument(
        "--noise_level",
        type=float,
        default=0.01,
        help="Nível de ruído (probabilidade)"
    )

    # Otimização
    parser.add_argument(
        "--learning_rate",
        type=float,
        default=0.1,
        help="Taxa de aprendizado"
    )
    parser.add_argument(
        "--max_iterations",
        type=int,
        default=100,
        help="Número máximo de iterações"
    )
    parser.add_argument(
        "--opt_method",
        type=str,
        default="adam",
        choices=["adam", "spsa", "cobyla", "l_bfgs_b", "differential_evolution", "bayesian"],
        help="Método de otimização"
    )

    # Mitigação de erro
    parser.add_argument(
        "--error_mitigation",
        type=str,
        default="zne",
        choices=["zne", "trex", "auec", "readout_mitigation", "none"],
        help="Técnica de mitigação de erro"
    )

    # Shots
    parser.add_argument(
        "--n_shots",
        type=int,
        default=1024,
        help="Número de shots para simulação"
    )

    # Resultados
    parser.add_argument(
        "--results_dir",
        type=str,
        default="./results_quantum_v8",
        help="Diretório para salvar resultados"
    )

    # Seed
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Seed aleatória"
    )

    # Verbose
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Modo verbose"
    )

    return parser.parse_args()


def create_config_from_args(args) -> ExperimentConfig:
    """Cria configuração a partir dos argumentos."""

    framework_map = {
        "pennylane": FrameworkType.PENNYLANE,
        "qiskit": FrameworkType.QISKIT,
        "cirq": FrameworkType.CIRQ,
    }

    noise_model_map = {
        "depolarizing": NoiseModel.DEPOLARIZING,
        "amplitude_damping": NoiseModel.AMPLITUDE_DAMPING,
        "phase_damping": NoiseModel.PHASE_DAMPING,
        "pauli": NoiseModel.PAULI,
        "none": NoiseModel.NONE,
    }

    opt_method_map = {
        "adam": OptimizationMethod.ADAM,
        "spsa": OptimizationMethod.SPSA,
        "cobyla": OptimizationMethod.COBYLA,
        "l_bfgs_b": OptimizationMethod.L_BFGS_B,
        "differential_evolution": OptimizationMethod.DIFFERENTIAL_EVOLUTION,
        "bayesian": OptimizationMethod.BAYESIAN,
    }

    error_mitigation_map = {
        "zne": ErrorMitigationTechnique.ZNE,
        "trex": ErrorMitigationTechnique.TREX,
        "auec": ErrorMitigationTechnique.AUEC,
        "readout_mitigation": ErrorMitigationTechnique.READOUT_MITIGATION,
        "none": ErrorMitigationTechnique.NONE,
    }

    n_parameters = args.n_qubits * args.n_layers * 3

    config = ExperimentConfig(
        framework=framework_map[args.framework],
        circuit_config=QuantumCircuitConfig(
            n_qubits=args.n_qubits,
            n_layers=args.n_layers,
            n_parameters=n_parameters,
            entanglement=args.entanglement
        ),
        noise_config=NoiseConfig(
            noise_model=noise_model_map[args.noise_model],
            noise_level=args.noise_level
        ),
        optimization_config=OptimizationConfig(
            method=opt_method_map[args.opt_method],
            learning_rate=args.learning_rate,
            max_iterations=args.max_iterations
        ),
        error_mitigation_config=ErrorMitigationConfig(
            technique=error_mitigation_map[args.error_mitigation]
        ),
        dataset_name=args.dataset,
        n_shots=args.n_shots,
        seed=args.seed
    )

    return config


def main():
    """Função principal."""
    args = parse_arguments()

    logger.info("=" * 80)
    logger.info("FRAMEWORK QUANTUM ADVANCED V8")
    logger.info("=" * 80)
    logger.info(f"Args: {args}")

    # Cria configuração
    config = create_config_from_args(args)

    # Cria runner
    runner = QuantumExperimentRunner(config, results_dir=args.results_dir)

    # Executa experimento
    results = runner.run_full_experiment()

    # Salva resultados
    runner.save_results("results_quantum_v8.json")
    runner.save_plots()

    # Imprime resumo
    logger.info("\n" + "=" * 80)
    logger.info("RESULTADOS FINAIS")
    logger.info("=" * 80)

    summary = {
        k: v for k, v in results.items()
        if k not in ['inference_metrics', 'training_metrics']
    }

    print(json.dumps(summary, indent=2, default=str))

    # Imprime arquivo salvo
    results_file = Path(args.results_dir) / "results_quantum_v8.json"
    logger.info(f"\nResultados salvos em: {results_file}")


if __name__ == "__main__":
    main()
