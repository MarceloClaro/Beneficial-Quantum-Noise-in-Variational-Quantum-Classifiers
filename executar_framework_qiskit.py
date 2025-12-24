"""
Script para executar o framework Qiskit completo e gerar visualizações.
Executa experimentos e salva resultados em formato adequado para documentação.
"""

import os
import sys
import logging

# Configurar logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Importar framework
from framework_qiskit import (
    ClassificadorVQCQiskit,
    carregar_datasets,
    executar_experimento_qiskit,
    visualizar_bloch_sphere,
    visualizar_state_city_3d,
    visualizar_qsphere,
)

def main():
    """Executa framework Qiskit completo."""
    
    # Criar diretório de resultados
    pasta_resultados = 'resultados_qiskit_framework'
    os.makedirs(pasta_resultados, exist_ok=True)
    os.makedirs(os.path.join(pasta_resultados, 'visualizacoes'), exist_ok=True)
    
    logger.info("="*80)
    logger.info("EXECUTANDO FRAMEWORK QISKIT COMPLETO")
    logger.info("="*80)
    
    # 1. Experimento básico com visualizações
    logger.info("\n1. Experimento Básico - Dataset Moons")
    logger.info("-" * 80)
    
    resultado = executar_experimento_qiskit(
        dataset_nome='moons',
        arquitetura='strongly_entangling',
        tipo_ruido='phase_damping',
        nivel_ruido=0.005,
        n_epocas=15,
        pasta_resultados=pasta_resultados,
        verbose=True
    )
    
    logger.info(f"✓ Acurácia teste: {resultado['acuracia_teste']:.4f}")
    logger.info(f"✓ Tempo de treino: {resultado['tempo_treino']:.2f}s")
    
    # 2. Experimento com ruído depolarizante
    logger.info("\n2. Experimento - Ruído Depolarizante")
    logger.info("-" * 80)
    
    datasets = carregar_datasets()
    dataset = datasets['circles']
    
    vqc = ClassificadorVQCQiskit(
        n_qubits=4,
        n_camadas=2,
        arquitetura='hardware_efficient',
        tipo_ruido='depolarizante',
        nivel_ruido=0.01,
        n_epocas=12,
        seed=42
    )
    
    vqc.fit(dataset['X_train'], dataset['y_train'])
    acc_circles = vqc.score(dataset['X_test'], dataset['y_test'])
    
    logger.info(f"✓ Acurácia (Circles): {acc_circles:.4f}")
    
    # Gerar visualizações adicionais
    logger.info("\n3. Gerando Visualizações Adicionais")
    logger.info("-" * 80)
    
    x_exemplo = dataset['X_test'][0]
    
    # Bloch sphere adicional
    bloch_path_2 = os.path.join(pasta_resultados, 'visualizacoes', 'bloch_sphere_depolarizing.png')
    visualizar_bloch_sphere(vqc, x_exemplo, bloch_path_2)
    logger.info(f"✓ Bloch sphere salva: {bloch_path_2}")
    
    # State city 3D adicional
    city_path_2 = os.path.join(pasta_resultados, 'visualizacoes', 'state_city_3d_depolarizing.png')
    visualizar_state_city_3d(vqc, x_exemplo, city_path_2)
    logger.info(f"✓ State City 3D salva: {city_path_2}")
    
    # Q-sphere adicional
    qsphere_path_2 = os.path.join(pasta_resultados, 'visualizacoes', 'qsphere_depolarizing.png')
    visualizar_qsphere(vqc, x_exemplo, qsphere_path_2)
    logger.info(f"✓ Q-Sphere salva: {qsphere_path_2}")
    
    # 3. Experimento comparativo - sem ruído vs com ruído
    logger.info("\n4. Comparação: Sem Ruído vs. Com Ruído")
    logger.info("-" * 80)
    
    dataset_iris = datasets['iris']
    
    # Sem ruído
    vqc_sem = ClassificadorVQCQiskit(
        n_qubits=4,
        n_camadas=2,
        arquitetura='basico',
        tipo_ruido='sem_ruido',
        n_epocas=10,
        seed=42
    )
    vqc_sem.fit(dataset_iris['X_train'], dataset_iris['y_train'])
    acc_sem = vqc_sem.score(dataset_iris['X_test'], dataset_iris['y_test'])
    
    # Com ruído phase damping
    vqc_com = ClassificadorVQCQiskit(
        n_qubits=4,
        n_camadas=2,
        arquitetura='basico',
        tipo_ruido='phase_damping',
        nivel_ruido=0.005,
        n_epocas=10,
        seed=42
    )
    vqc_com.fit(dataset_iris['X_train'], dataset_iris['y_train'])
    acc_com = vqc_com.score(dataset_iris['X_test'], dataset_iris['y_test'])
    
    logger.info(f"✓ Sem ruído: {acc_sem:.4f}")
    logger.info(f"✓ Com ruído (γ=0.005): {acc_com:.4f}")
    logger.info(f"✓ Delta: {acc_com - acc_sem:+.4f}")
    
    # Visualizações comparativas
    x_iris = dataset_iris['X_test'][0]
    
    # Circuito sem ruído
    circuit_path_sem = os.path.join(pasta_resultados, 'visualizacoes', 'circuito_sem_ruido.png')
    vqc_sem.get_circuit_diagram(circuit_path_sem)
    logger.info(f"✓ Circuito (sem ruído) salvo: {circuit_path_sem}")
    
    # Circuito com ruído
    circuit_path_com = os.path.join(pasta_resultados, 'visualizacoes', 'circuito_com_ruido.png')
    vqc_com.get_circuit_diagram(circuit_path_com)
    logger.info(f"✓ Circuito (com ruído) salvo: {circuit_path_com}")
    
    # Resumo final
    logger.info("\n" + "="*80)
    logger.info("EXECUÇÃO COMPLETA - RESUMO")
    logger.info("="*80)
    logger.info(f"✓ Pasta de resultados: {pasta_resultados}")
    logger.info(f"✓ Total de visualizações geradas: 8")
    logger.info(f"✓ Experimentos executados: 4")
    logger.info("\nArquivos gerados:")
    logger.info(f"  - {pasta_resultados}/circuito_qiskit.png")
    logger.info(f"  - {pasta_resultados}/bloch_sphere_qiskit.png")
    logger.info(f"  - {pasta_resultados}/state_city_3d_qiskit.png")
    logger.info(f"  - {pasta_resultados}/visualizacoes/bloch_sphere_depolarizing.png")
    logger.info(f"  - {pasta_resultados}/visualizacoes/state_city_3d_depolarizing.png")
    logger.info(f"  - {pasta_resultados}/visualizacoes/qsphere_depolarizing.png")
    logger.info(f"  - {pasta_resultados}/visualizacoes/circuito_sem_ruido.png")
    logger.info(f"  - {pasta_resultados}/visualizacoes/circuito_com_ruido.png")
    
    logger.info("\n✓ Framework Qiskit executado com sucesso!")
    logger.info("="*80)
    
    return {
        'pasta_resultados': pasta_resultados,
        'resultado_moons': resultado,
        'acuracia_circles': acc_circles,
        'acuracia_iris_sem_ruido': acc_sem,
        'acuracia_iris_com_ruido': acc_com
    }

if __name__ == "__main__":
    try:
        resultados = main()
        sys.exit(0)
    except Exception as e:
        logger.error(f"Erro durante execução: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
