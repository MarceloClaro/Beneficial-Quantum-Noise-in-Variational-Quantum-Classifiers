#!/usr/bin/env python3
# =============================================================================
# EXEMPLOS PRÁTICOS - Framework Quantum Advanced V8
# =============================================================================
"""
Exemplos detalhados de como usar o Framework Quantum Advanced V8
para diferentes cenários de pesquisa e produção.
"""

import json
import numpy as np
from pathlib import Path

from framework_quantum_advanced_v8 import (
    ExperimentConfig,
    QuantumCircuitConfig,
    NoiseConfig,
    OptimizationConfig,
    ErrorMitigationConfig,
    QuantumExperimentRunner,
    QuantumComplexityAnalyzer,
    ZeroNoiseExtrapolation,
    NoiseValidationFramework,
    QuantumAlgorithmBenchmark,
    DeepChemDatasetLoader,
    FrameworkType,
    NoiseModel,
    OptimizationMethod,
    ErrorMitigationTechnique,
)


# =============================================================================
# EXEMPLO 1: Experimento Básico com Iris
# =============================================================================
def exemplo_1_basico():
    """Experimento simples com dataset Iris para validação."""
    print("\n" + "="*80)
    print("EXEMPLO 1: Experimento Básico com Iris")
    print("="*80)
    
    config = ExperimentConfig(
        framework=FrameworkType.PENNYLANE,
        circuit_config=QuantumCircuitConfig(
            n_qubits=4,
            n_layers=2,
            n_parameters=24,
            architecture="strongly_entangling",
            entanglement="full"
        ),
        noise_config=NoiseConfig(
            noise_model=NoiseModel.DEPOLARIZING,
            noise_level=0.01,
            gate_error=0.001,
            readout_error=0.01
        ),
        optimization_config=OptimizationConfig(
            method=OptimizationMethod.ADAM,
            learning_rate=0.1,
            max_iterations=50,  # Reduzido para demo
            early_stopping=True
        ),
        error_mitigation_config=ErrorMitigationConfig(
            technique=ErrorMitigationTechnique.ZNE,
            zne_scale_factors=[1.0, 1.5, 2.0]
        ),
        dataset_name="iris",
        n_shots=1024,
        seed=42
    )
    
    runner = QuantumExperimentRunner(config, results_dir="./results_ex1_basico")
    results = runner.run_full_experiment()
    runner.save_results("results_ex1.json")
    runner.save_plots()
    
    print(f"\n✓ Experimento completado em {results['execution_time_seconds']:.2f}s")
    print(f"✓ Acurácia: {results.get('inference_metrics', {}).get('accuracy', 'N/A'):.4f}")
    
    return results


# =============================================================================
# EXEMPLO 2: Experimento com Dataset HIV (DeepChem)
# =============================================================================
def exemplo_2_hiv_dataset():
    """Experimento com dataset HIV para descoberta molecular."""
    print("\n" + "="*80)
    print("EXEMPLO 2: Classificação HIV com DeepChem Dataset")
    print("="*80)
    
    config = ExperimentConfig(
        framework=FrameworkType.PENNYLANE,
        circuit_config=QuantumCircuitConfig(
            n_qubits=6,
            n_layers=3,
            n_parameters=54,
            architecture="strongly_entangling",
            entanglement="full"
        ),
        noise_config=NoiseConfig(
            noise_model=NoiseModel.AMPLITUDE_DAMPING,
            noise_level=0.02,
            gate_error=0.002,
            readout_error=0.02,
            t1_time=100.0,
            t2_time=50.0
        ),
        optimization_config=OptimizationConfig(
            method=OptimizationMethod.ADAM,
            learning_rate=0.05,
            max_iterations=100,
            early_stopping=True,
            early_stopping_patience=20
        ),
        error_mitigation_config=ErrorMitigationConfig(
            technique=ErrorMitigationTechnique.ZNE,
            zne_scale_factors=[1.0, 1.5, 2.0, 2.5],
            zne_extrapolation_type="exponential"
        ),
        dataset_name="hiv",
        n_shots=2048,
        seed=42
    )
    
    runner = QuantumExperimentRunner(config, results_dir="./results_ex2_hiv")
    
    # Análise de complexidade antes de rodar
    complexity = runner.run_complexity_analysis()
    print("\n📊 Análise de Complexidade:")
    print(f"  Circuit Depth: {complexity['circuit_depth']}")
    print(f"  Total Gates: {complexity['gate_count']['total']}")
    print(f"  Barren Plateau Prob: {complexity['barren_plateau_probability']:.4f}")
    print(f"  Est. Time: {complexity['total_estimated_time_s']:.2f}s")
    
    # Executar experimento
    results = runner.run_full_experiment()
    runner.save_results("results_ex2_hiv.json")
    runner.save_plots()
    
    return results


# =============================================================================
# EXEMPLO 3: Validação de Fórmula de Ruído
# =============================================================================
def exemplo_3_validacao_ruido():
    """Validação da fórmula de predição de ruído."""
    print("\n" + "="*80)
    print("EXEMPLO 3: Validação de Fórmula de Predição de Ruído")
    print("="*80)
    
    validator = NoiseValidationFramework()
    
    # Simular dados experimentais
    noise_levels = [0.001, 0.005, 0.01, 0.02, 0.05]
    circuit_depths = [5, 10, 15, 20]
    
    print("\n🔍 Predizendo impacto de ruído em diferentes configurações:\n")
    
    predictions = {}
    for noise_level in noise_levels:
        predictions[noise_level] = []
        for depth in circuit_depths:
            fidelity = validator.predict_noise_impact(
                noise_level=noise_level,
                circuit_depth=depth,
                n_qubits=4
            )
            predictions[noise_level].append(fidelity)
            print(f"  Noise: {noise_level:.3f}, Depth: {depth:2d} → Predicted Fidelity: {fidelity:.4f}")
    
    # Simular dados "medidos"
    measured_errors = {}
    for noise_level in noise_levels:
        for depth in circuit_depths:
            # Simulando dados com pequeno ruído experimental
            predicted = validator.predict_noise_impact(noise_level, depth, 4)
            actual = predicted * (1 + np.random.normal(0, 0.05))  # 5% de incerteza
            measured_errors[(noise_level, depth)] = 1 - actual
    
    # Analisar escalamento
    print("\n📈 Análise de Escalamento:\n")
    scaling = validator.analyze_noise_scaling(
        noise_levels=noise_levels,
        circuit_depths=circuit_depths,
        measured_errors=measured_errors
    )
    
    if 'model' in scaling:
        print(f"  Modelo: {scaling['model']}")
        print(f"  Coeficiente de Ruído: {scaling['noise_coefficient']:.4f}")
        print(f"  Coeficiente de Profundidade: {scaling['depth_coefficient']:.4f}")
        print(f"  R²: {scaling['r_squared']:.4f}")
    
    return predictions, scaling


# =============================================================================
# EXEMPLO 4: Zero-Noise Extrapolation (ZNE)
# =============================================================================
def exemplo_4_zne():
    """Demonstração de Zero-Noise Extrapolation."""
    print("\n" + "="*80)
    print("EXEMPLO 4: Zero-Noise Extrapolation (ZNE)")
    print("="*80)
    
    zne_config = ErrorMitigationConfig(
        technique=ErrorMitigationTechnique.ZNE,
        zne_scale_factors=[1.0, 1.5, 2.0, 2.5, 3.0],
        zne_extrapolation_type="exponential"
    )
    
    zne = ZeroNoiseExtrapolation(zne_config)
    
    # Simular função observable que degrada com escala de ruído
    def observable_with_noise(scale):
        """Observable esperado que degrada exponencialmente com ruído."""
        noise_effect = np.exp(-scale)
        clean_value = 0.95
        noisy_value = clean_value * noise_effect + 0.05 * (1 - noise_effect)
        return noisy_value
    
    # Aplicar ZNE
    print("\n🔄 Aplicando Zero-Noise Extrapolation:\n")
    extrapolated, details = zne.apply_zne(observable_with_noise)
    
    print(f"  Valores medidos: {[f'{v:.4f}' for v in details['measured_values']]}")
    print(f"  Fatores de escala: {details['scales']}")
    print(f"  Tipo de extrapolação: {details['extrapolation_type']}")
    print(f"\n  ✓ Valor extrapolado (sem ruído): {extrapolated:.4f}")
    print(f"  ✓ Melhoria: {(extrapolated - details['measured_values'][0]) / details['measured_values'][0] * 100:.2f}%")
    
    return extrapolated, details


# =============================================================================
# EXEMPLO 5: Benchmarking vs Algoritmos Clássicos
# =============================================================================
def exemplo_5_benchmarking():
    """Benchmarking contra algoritmos clássicos."""
    print("\n" + "="*80)
    print("EXEMPLO 5: Benchmarking Quântico vs Clássico")
    print("="*80)
    
    benchmark = QuantumAlgorithmBenchmark()
    
    # Simular predições
    np.random.seed(42)
    y_true = np.random.binomial(1, 0.5, 100)
    
    # Predições VQC (melhor performance)
    vqc_predictions = y_true.astype(float) + np.random.normal(0, 0.15, 100)
    vqc_predictions = np.clip(vqc_predictions, 0, 1)
    
    # Predições clássicas (baseline)
    classical_predictions = y_true.astype(float) + np.random.normal(0, 0.25, 100)
    classical_predictions = np.clip(classical_predictions, 0, 1)
    
    # Benchmark
    comparison = benchmark.benchmark_against_classical(
        vqc_predictions,
        classical_predictions,
        y_true
    )
    
    print("\n📊 Comparação VQC vs Clássico:\n")
    print("Métrica              | VQC    | Clássico | Melhoria")
    print("-" * 55)
    
    for metric in ['accuracy', 'precision', 'recall', 'f1']:
        vqc_val = comparison['vqc_metrics'].get(metric, 0)
        classical_val = comparison['classical_metrics'].get(metric, 0)
        improvement = comparison['improvement_ratio'].get(metric, 0) * 100
        
        print(f"{metric.capitalize():20s} | {vqc_val:6.4f} | {classical_val:8.4f} | {improvement:+6.1f}%")
    
    print(f"\n✓ VQC ganhou em {comparison['vqc_wins']} métricas")
    
    return comparison


# =============================================================================
# EXEMPLO 6: Análise de Escalamento
# =============================================================================
def exemplo_6_escalamento():
    """Análise de escalamento com tamanho do sistema."""
    print("\n" + "="*80)
    print("EXEMPLO 6: Análise de Escalamento")
    print("="*80)
    
    benchmark = QuantumAlgorithmBenchmark()
    
    # Simular tempos de execução para diferentes tamanhos
    system_sizes = [2, 3, 4, 5, 6]
    # Escalamento exponencial típico de simulação clássica
    execution_times = [0.01 * (2 ** n) for n in system_sizes]
    
    scaling = benchmark.benchmark_scaling(
        system_sizes=system_sizes,
        execution_times=execution_times
    )
    
    print("\n📈 Escalamento de Recursos:\n")
    print("System Size | Hilbert Space | Exec Time (ms) | Scaling")
    print("-" * 60)
    
    for size, hilbert, time in zip(
        scaling['system_sizes'],
        scaling['quantum_state_dimension'],
        scaling['execution_times']
    ):
        print(f"{size:11d} | {hilbert:13d} | {time*1000:13.2f} | -")
    
    print(f"\n  Expoente de escalamento: {scaling['scaling_exponent']:.2f}")
    print(f"  Tipo: {scaling['scaling_type']}")
    print(f"  Modelo: {scaling['scaling_fit']}")
    
    return scaling


# =============================================================================
# EXEMPLO 7: Análise de Complexidade Detalhada
# =============================================================================
def exemplo_7_complexidade():
    """Análise detalhada de complexidade quântica."""
    print("\n" + "="*80)
    print("EXEMPLO 7: Análise de Complexidade Quântica Detalhada")
    print("="*80)
    
    analyzer = QuantumComplexityAnalyzer()
    
    configs = [
        ("Pequeno", 3, 2, "linear"),
        ("Médio", 5, 3, "full"),
        ("Grande", 8, 4, "full"),
    ]
    
    print("\n🔧 Comparação de Configurações:\n")
    print("Config    | Qubits | Layers | Entanglement | Depth | Gates | BP Prob")
    print("-" * 75)
    
    for name, n_qubits, n_layers, entanglement in configs:
        circuit_config = QuantumCircuitConfig(
            n_qubits=n_qubits,
            n_layers=n_layers,
            n_parameters=n_qubits * n_layers * 3,
            entanglement=entanglement
        )
        
        depth = analyzer.calculate_circuit_depth(n_qubits, n_layers, entanglement)
        gates = analyzer.calculate_gate_count(n_qubits, n_layers, entanglement)
        bp_prob = analyzer.estimate_barren_plateau_probability(n_qubits, depth, 0.01)
        
        print(f"{name:9s} | {n_qubits:6d} | {n_layers:6d} | {entanglement:12s} | {depth:5d} | {gates['total']:5d} | {bp_prob:7.4f}")
    
    return configs


# =============================================================================
# EXEMPLO 8: Carregamento de Datasets DeepChem
# =============================================================================
def exemplo_8_deepchem():
    """Demonstração de carregamento de datasets DeepChem."""
    print("\n" + "="*80)
    print("EXEMPLO 8: Datasets DeepChem")
    print("="*80)
    
    loader = DeepChemDatasetLoader()
    
    datasets_to_try = ["hiv", "malaria", "tb"]
    
    print("\n💊 Datasets Moleculares Disponíveis:\n")
    
    for dataset_name in datasets_to_try:
        try:
            X, y = loader.load_dataset(dataset_name)
            print(f"✓ {dataset_name.upper():10s}: shape={X.shape}, labels={np.unique(y)}")
        except Exception as e:
            print(f"✗ {dataset_name.upper():10s}: {str(e)[:50]}")
    
    return


# =============================================================================
# EXEMPLO 9: Experimento Completo com Todas as Técnicas
# =============================================================================
def exemplo_9_completo():
    """Experimento completo utilizando todas as técnicas."""
    print("\n" + "="*80)
    print("EXEMPLO 9: Experimento Completo com Todas as Técnicas")
    print("="*80)
    
    config = ExperimentConfig(
        framework=FrameworkType.PENNYLANE,
        circuit_config=QuantumCircuitConfig(
            n_qubits=5,
            n_layers=3,
            n_parameters=45,
            architecture="strongly_entangling",
            entanglement="full"
        ),
        noise_config=NoiseConfig(
            noise_model=NoiseModel.DEPOLARIZING,
            noise_level=0.02,
            gate_error=0.002,
            readout_error=0.02,
            t1_time=100.0,
            t2_time=50.0
        ),
        optimization_config=OptimizationConfig(
            method=OptimizationMethod.ADAM,
            learning_rate=0.08,
            max_iterations=100,
            early_stopping=True,
            early_stopping_patience=25
        ),
        error_mitigation_config=ErrorMitigationConfig(
            technique=ErrorMitigationTechnique.ZNE,
            zne_scale_factors=[1.0, 1.5, 2.0, 2.5],
            zne_extrapolation_type="exponential"
        ),
        dataset_name="wine",
        n_shots=1024,
        seed=42
    )
    
    print("\n⚙️ Configuração do Experimento:")
    print(f"  Framework: {config.framework.value}")
    print(f"  Qubits: {config.circuit_config.n_qubits}")
    print(f"  Layers: {config.circuit_config.n_layers}")
    print(f"  Ruído: {config.noise_config.noise_model.value} ({config.noise_config.noise_level:.4f})")
    print(f"  Mitigação: {config.error_mitigation_config.technique.value}")
    
    runner = QuantumExperimentRunner(config, results_dir="./results_ex9_completo")
    results = runner.run_full_experiment()
    runner.save_results("results_ex9.json")
    runner.save_plots()
    
    return results


# =============================================================================
# EXECUTOR PRINCIPAL
# =============================================================================
def main():
    """Executa todos os exemplos."""
    
    print("\n" + "="*80)
    print("EXEMPLOS PRÁTICOS - FRAMEWORK QUANTUM ADVANCED V8")
    print("="*80)
    
    # Escolher qual exemplo executar
    examples = {
        "1": ("Experimento Básico com Iris", exemplo_1_basico),
        "2": ("HIV Dataset (DeepChem)", exemplo_2_hiv_dataset),
        "3": ("Validação de Ruído", exemplo_3_validacao_ruido),
        "4": ("Zero-Noise Extrapolation", exemplo_4_zne),
        "5": ("Benchmarking", exemplo_5_benchmarking),
        "6": ("Escalamento", exemplo_6_escalamento),
        "7": ("Complexidade", exemplo_7_complexidade),
        "8": ("DeepChem Datasets", exemplo_8_deepchem),
        "9": ("Experimento Completo", exemplo_9_completo),
        "all": ("Executar Todos", None),
    }
    
    print("\nExemplos Disponíveis:\n")
    for key, (name, _) in examples.items():
        print(f"  {key}. {name}")
    
    choice = input("\nEscolha um exemplo (1-9, 'all' para todos, ou 'q' para sair): ").strip()
    
    if choice.lower() == 'q':
        print("Encerrando...")
        return
    
    if choice == 'all':
        # Executar todos os exemplos (exceto DeepChem se não disponível)
        for key in ['1', '3', '4', '5', '6', '7', '9']:
            if key in examples and examples[key][1] is not None:
                try:
                    print("\n")
                    examples[key][1]()
                except Exception as e:
                    print(f"❌ Erro ao executar exemplo {key}: {e}")
    else:
        if choice in examples and examples[choice][1] is not None:
            try:
                examples[choice][1]()
            except Exception as e:
                print(f"❌ Erro: {e}")
        else:
            print("Opção inválida")
    
    print("\n✓ Exemplos concluídos!")


if __name__ == "__main__":
    main()
