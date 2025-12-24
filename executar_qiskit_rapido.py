"""
Script otimizado para gerar visualizações rapidamente.
"""

import os
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

from framework_qiskit import (
    ClassificadorVQCQiskit,
    carregar_datasets,
    visualizar_bloch_sphere,
    visualizar_state_city_3d,
    visualizar_qsphere,
)

def main():
    pasta_resultados = 'resultados_qiskit_framework'
    os.makedirs(pasta_resultados, exist_ok=True)
    os.makedirs(os.path.join(pasta_resultados, 'visualizacoes'), exist_ok=True)
    
    logger.info("="*80)
    logger.info("FRAMEWORK QISKIT - GERAÇÃO RÁPIDA DE VISUALIZAÇÕES")
    logger.info("="*80)
    
    # Carregar datasets
    datasets = carregar_datasets()
    dataset = datasets['moons']
    
    # Experimento 1: Com ruído phase damping
    logger.info("\n1. Treinando VQC com Phase Damping...")
    vqc1 = ClassificadorVQCQiskit(
        n_qubits=4,
        n_camadas=2,
        arquitetura='strongly_entangling',
        tipo_ruido='phase_damping',
        nivel_ruido=0.005,
        n_epocas=5,  # Reduzido para velocidade
        seed=42
    )
    
    vqc1.fit(dataset['X_train'], dataset['y_train'])
    acc1 = vqc1.score(dataset['X_test'], dataset['y_test'])
    logger.info(f"✓ Acurácia: {acc1:.4f}")
    
    # Gerar visualizações
    x_exemplo = dataset['X_test'][0]
    
    logger.info("\n2. Gerando Visualizações...")
    
    # Circuito
    circuit_path = os.path.join(pasta_resultados, 'circuito_qiskit.png')
    vqc1.get_circuit_diagram(circuit_path)
    logger.info(f"✓ Circuito salvo: {circuit_path}")
    
    # Bloch sphere
    bloch_path = os.path.join(pasta_resultados, 'bloch_sphere_qiskit.png')
    visualizar_bloch_sphere(vqc1, x_exemplo, bloch_path)
    logger.info(f"✓ Bloch sphere salva: {bloch_path}")
    
    # State City 3D
    city_path = os.path.join(pasta_resultados, 'state_city_3d_qiskit.png')
    visualizar_state_city_3d(vqc1, x_exemplo, city_path)
    logger.info(f"✓ State City 3D salva: {city_path}")
    
    # Q-Sphere
    qsphere_path = os.path.join(pasta_resultados, 'qsphere_qiskit.png')
    visualizar_qsphere(vqc1, x_exemplo, qsphere_path)
    logger.info(f"✓ Q-Sphere salva: {qsphere_path}")
    
    # Experimento 2: Sem ruído (comparação)
    logger.info("\n3. Treinando VQC sem Ruído (comparação)...")
    vqc2 = ClassificadorVQCQiskit(
        n_qubits=4,
        n_camadas=2,
        arquitetura='basico',
        tipo_ruido='sem_ruido',
        n_epocas=5,
        seed=42
    )
    
    vqc2.fit(dataset['X_train'], dataset['y_train'])
    acc2 = vqc2.score(dataset['X_test'], dataset['y_test'])
    logger.info(f"✓ Acurácia: {acc2:.4f}")
    
    # Circuito sem ruído
    circuit_path_2 = os.path.join(pasta_resultados, 'visualizacoes', 'circuito_sem_ruido.png')
    vqc2.get_circuit_diagram(circuit_path_2)
    logger.info(f"✓ Circuito (sem ruído) salvo: {circuit_path_2}")
    
    # Bloch adicional
    bloch_path_2 = os.path.join(pasta_resultados, 'visualizacoes', 'bloch_sphere_sem_ruido.png')
    visualizar_bloch_sphere(vqc2, x_exemplo, bloch_path_2)
    logger.info(f"✓ Bloch sphere (sem ruído) salva: {bloch_path_2}")
    
    logger.info("\n" + "="*80)
    logger.info("EXECUÇÃO COMPLETA")
    logger.info("="*80)
    logger.info(f"✓ Total de visualizações: 6")
    logger.info(f"✓ Acurácia com ruído (γ=0.005): {acc1:.4f}")
    logger.info(f"✓ Acurácia sem ruído: {acc2:.4f}")
    logger.info(f"✓ Delta: {acc1 - acc2:+.4f}")
    logger.info(f"\nArquivos em: {pasta_resultados}/")
    logger.info("="*80)

if __name__ == "__main__":
    main()
