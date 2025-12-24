"""
Script Rápido de Demonstração Qiskit com Visualizações
Gera visualizações de exemplo para documentação em ~10-15 minutos
"""

import sys
import os

# Verificar e instalar dependências se necessário
try:
    import numpy as np
    import pandas as pd
    from sklearn import datasets as sk_datasets
    from sklearn.model_selection import train_test_split
    import qiskit
    from qiskit_aer import AerSimulator
    print("✓ Todas as dependências estão instaladas")
except ImportError as e:
    print(f"❌ Faltam dependências: {e}")
    print("Instalando dependências...")
    os.system("pip install -q numpy pandas scikit-learn qiskit qiskit-aer matplotlib")
    print("✓ Dependências instaladas")
    # Reimportar
    import numpy as np
    import pandas as pd
    from sklearn import datasets as sk_datasets
    from sklearn.model_selection import train_test_split
    import qiskit
    from qiskit_aer import AerSimulator

from framework_qiskit import (
    ClassificadorVQCQiskit,
    visualizar_bloch_sphere,
    visualizar_state_city_3d,
    visualizar_qsphere,
)

print("\n" + "=" * 80)
print("DEMONSTRAÇÃO RÁPIDA - FRAMEWORK QISKIT")
print("=" * 80)

# Criar pasta de resultados
pasta_demo = "demo_qiskit_visualizacoes"
os.makedirs(pasta_demo, exist_ok=True)
print(f"\n📁 Pasta de resultados: {pasta_demo}/")

# Carregar dataset simples (moons)
print("\n1️⃣ Carregando dataset Moons...")
X, y = sk_datasets.make_moons(n_samples=100, noise=0.1, random_state=42)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
print(f"   ✓ Train: {len(X_train)} amostras, Test: {len(X_test)} amostras")

# Experimento 1: Sem ruído
print("\n2️⃣ Experimento 1: Circuito básico SEM ruído")
print("   Treinando VQC...")
vqc1 = ClassificadorVQCQiskit(
    n_qubits=4,
    n_camadas=2,
    arquitetura='basico',
    estrategia_init='aleatorio',
    tipo_ruido='sem_ruido',
    nivel_ruido=0.0,
    n_epocas=3,  # Apenas 3 épocas para velocidade
    batch_size=32,
    seed=42,
    shots=256  # Shots reduzidos para velocidade
)
vqc1.fit(X_train, y_train)
acc1 = vqc1.score(X_test, y_test)
print(f"   ✓ Acurácia: {acc1:.4f}")

print("   Gerando visualizações...")
try:
    x_sample = X_test[0]
    
    # Bloch Sphere
    bloch1 = os.path.join(pasta_demo, "exp1_sem_ruido_bloch.png")
    visualizar_bloch_sphere(vqc1, x_sample, bloch1)
    print(f"   ✓ Bloch sphere: {bloch1}")
    
    # State City 3D
    city1 = os.path.join(pasta_demo, "exp1_sem_ruido_city3d.png")
    visualizar_state_city_3d(vqc1, x_sample, city1)
    print(f"   ✓ State City 3D: {city1}")
    
    # Q-Sphere
    qsphere1 = os.path.join(pasta_demo, "exp1_sem_ruido_qsphere.png")
    visualizar_qsphere(vqc1, x_sample, qsphere1)
    print(f"   ✓ Q-Sphere: {qsphere1}")
    
    # Circuito
    circuit1 = os.path.join(pasta_demo, "exp1_sem_ruido_circuit.png")
    vqc1.get_circuit_diagram(circuit1)
    print(f"   ✓ Circuito: {circuit1}")
    
except Exception as e:
    print(f"   ⚠️ Erro ao gerar visualizações: {e}")

# Experimento 2: Com ruído phase damping (benéfico)
print("\n3️⃣ Experimento 2: Com ruído PHASE DAMPING (γ=0.005)")
print("   Treinando VQC...")
vqc2 = ClassificadorVQCQiskit(
    n_qubits=4,
    n_camadas=2,
    arquitetura='strongly_entangling',
    estrategia_init='matematico',
    tipo_ruido='phase_damping',
    nivel_ruido=0.005,
    n_epocas=3,
    batch_size=32,
    seed=42,
    shots=256
)
vqc2.fit(X_train, y_train)
acc2 = vqc2.score(X_test, y_test)
print(f"   ✓ Acurácia: {acc2:.4f}")

print("   Gerando visualizações...")
try:
    # Bloch Sphere
    bloch2 = os.path.join(pasta_demo, "exp2_phase_damping_bloch.png")
    visualizar_bloch_sphere(vqc2, x_sample, bloch2)
    print(f"   ✓ Bloch sphere: {bloch2}")
    
    # State City 3D
    city2 = os.path.join(pasta_demo, "exp2_phase_damping_city3d.png")
    visualizar_state_city_3d(vqc2, x_sample, city2)
    print(f"   ✓ State City 3D: {city2}")
    
    # Q-Sphere
    qsphere2 = os.path.join(pasta_demo, "exp2_phase_damping_qsphere.png")
    visualizar_qsphere(vqc2, x_sample, qsphere2)
    print(f"   ✓ Q-Sphere: {qsphere2}")
    
    # Circuito
    circuit2 = os.path.join(pasta_demo, "exp2_phase_damping_circuit.png")
    vqc2.get_circuit_diagram(circuit2)
    print(f"   ✓ Circuito: {circuit2}")
    
except Exception as e:
    print(f"   ⚠️ Erro ao gerar visualizações: {e}")

# Experimento 3: Com ruído amplitude damping
print("\n4️⃣ Experimento 3: Com ruído AMPLITUDE DAMPING (γ=0.005)")
print("   Treinando VQC...")
vqc3 = ClassificadorVQCQiskit(
    n_qubits=4,
    n_camadas=2,
    arquitetura='hardware_efficient',
    estrategia_init='quantico',
    tipo_ruido='amplitude_damping',
    nivel_ruido=0.005,
    n_epocas=3,
    batch_size=32,
    seed=42,
    shots=256
)
vqc3.fit(X_train, y_train)
acc3 = vqc3.score(X_test, y_test)
print(f"   ✓ Acurácia: {acc3:.4f}")

print("   Gerando visualizações...")
try:
    # Apenas Bloch e City para este
    bloch3 = os.path.join(pasta_demo, "exp3_amplitude_damping_bloch.png")
    visualizar_bloch_sphere(vqc3, x_sample, bloch3)
    print(f"   ✓ Bloch sphere: {bloch3}")
    
    city3 = os.path.join(pasta_demo, "exp3_amplitude_damping_city3d.png")
    visualizar_state_city_3d(vqc3, x_sample, city3)
    print(f"   ✓ State City 3D: {city3}")
    
except Exception as e:
    print(f"   ⚠️ Erro ao gerar visualizações: {e}")

# Resumo
print("\n" + "=" * 80)
print("RESUMO DA DEMONSTRAÇÃO")
print("=" * 80)
print(f"\n📊 Resultados:")
print(f"   Exp 1 (Sem ruído):          Acurácia = {acc1:.4f}")
print(f"   Exp 2 (Phase damping):      Acurácia = {acc2:.4f}")
print(f"   Exp 3 (Amplitude damping):  Acurácia = {acc3:.4f}")

# Contar imagens geradas
import glob
imagens = glob.glob(os.path.join(pasta_demo, "*.png"))
print(f"\n🎨 Total de visualizações geradas: {len(imagens)}")
for img in sorted(imagens):
    print(f"   • {os.path.basename(img)}")

print(f"\n✅ Demonstração completa! Visualizações em: {pasta_demo}/")
print("=" * 80)

# Salvar resumo em CSV
resumo = pd.DataFrame({
    'experimento': ['Exp1_Sem_Ruido', 'Exp2_Phase_Damping', 'Exp3_Amplitude_Damping'],
    'arquitetura': ['basico', 'strongly_entangling', 'hardware_efficient'],
    'tipo_ruido': ['sem_ruido', 'phase_damping', 'amplitude_damping'],
    'nivel_ruido': [0.0, 0.005, 0.005],
    'acuracia': [acc1, acc2, acc3],
})
csv_path = os.path.join(pasta_demo, 'resumo_experimentos.csv')
resumo.to_csv(csv_path, index=False)
print(f"\n💾 Resumo salvo: {csv_path}")
