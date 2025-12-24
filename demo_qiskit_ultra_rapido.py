"""
Demo Qiskit Ultra-Rápido com Timeout de 600s
Gera visualizações essenciais para documentação
"""

import sys
import os
import time
from pathlib import Path

print("Instalando dependências mínimas...")
os.system("pip install -q numpy scikit-learn qiskit qiskit-aer matplotlib 2>&1 | tail -5")

import numpy as np
from sklearn import datasets as sk_datasets
from sklearn.model_selection import train_test_split

# Imports Qiskit minimalistas
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit import Parameter
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error, phase_damping_error
from qiskit.quantum_info import Statevector
from qiskit.visualization import plot_bloch_multivector, plot_state_city, circuit_drawer
import matplotlib
matplotlib.use('Agg')  # Backend sem display
import matplotlib.pyplot as plt

print("\n" + "=" * 80)
print("DEMO QISKIT - VISUALIZAÇÕES COM TIMEOUT 600s")
print("=" * 80)

# Timeout global
TIMEOUT = 600  # 600 segundos = 10 minutos
inicio_global = time.time()

def check_timeout():
    """Verifica se excedeu o timeout."""
    elapsed = time.time() - inicio_global
    if elapsed > TIMEOUT:
        print(f"\n⏱️ TIMEOUT atingido ({elapsed:.1f}s > {TIMEOUT}s)")
        return True
    return False

# Criar pasta
pasta = "visualizacoes_qiskit_md"
os.makedirs(pasta, exist_ok=True)
print(f"\n📁 Pasta: {pasta}/\n")

# Dataset simples
X, y = sk_datasets.make_moons(n_samples=50, noise=0.1, random_state=42)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

def criar_circuito_simples(x):
    """Cria circuito VQC simples."""
    qc = QuantumCircuit(2)
    
    # Encoding
    qc.ry(x[0] * np.pi, 0)
    qc.ry(x[1] * np.pi, 1)
    
    # Variational
    qc.ry(0.5, 0)
    qc.ry(0.3, 1)
    qc.cx(0, 1)
    qc.ry(0.7, 0)
    qc.ry(0.4, 1)
    
    return qc

# Visualização 1: Circuito Básico
print("1️⃣ Gerando circuito básico...")
try:
    qc_basico = criar_circuito_simples(X_test[0])
    fig = circuit_drawer(qc_basico, output='mpl', style='iqp')
    path1 = os.path.join(pasta, "01_circuito_basico.png")
    fig.savefig(path1, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"   ✓ {path1}")
except Exception as e:
    print(f"   ✗ Erro: {e}")

if check_timeout():
    sys.exit(0)

# Visualização 2: Bloch Sphere (sem ruído)
print("\n2️⃣ Gerando Bloch sphere (sem ruído)...")
try:
    qc = criar_circuito_simples(X_test[0])
    state = Statevector.from_instruction(qc)
    fig = plot_bloch_multivector(state)
    path2 = os.path.join(pasta, "02_bloch_sem_ruido.png")
    fig.savefig(path2, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"   ✓ {path2}")
except Exception as e:
    print(f"   ✗ Erro: {e}")

if check_timeout():
    sys.exit(0)

# Visualização 3: State City 3D (sem ruído)
print("\n3️⃣ Gerando State City 3D (sem ruído)...")
try:
    state = Statevector.from_instruction(qc)
    fig = plot_state_city(state)
    path3 = os.path.join(pasta, "03_city3d_sem_ruido.png")
    fig.savefig(path3, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"   ✓ {path3}")
except Exception as e:
    print(f"   ✗ Erro: {e}")

if check_timeout():
    sys.exit(0)

# Visualização 4: Bloch com ruído Phase Damping
print("\n4️⃣ Gerando Bloch sphere (COM ruído phase damping)...")
rho = None  # Definir globalmente
try:
    qc_noise = criar_circuito_simples(X_test[1])
    
    # Criar noise model com erros separados para 1-qubit e 2-qubit
    noise_model = NoiseModel()
    error_1q = phase_damping_error(0.01)
    error_2q = error_1q.tensor(error_1q)  # 2-qubit error
    noise_model.add_all_qubit_quantum_error(error_1q, ['ry'])
    noise_model.add_all_qubit_quantum_error(error_2q, ['cx'])
    
    # Simular com ruído
    simulator = AerSimulator(noise_model=noise_model, method='density_matrix')
    qc_noise.save_density_matrix()
    result = simulator.run(qc_noise, shots=1).result()
    rho = result.data()['density_matrix']
    
    # Plotar apenas o primeiro qubit
    from qiskit.quantum_info import partial_trace
    rho_qubit0 = partial_trace(rho, [1])
    
    # Converter para vetor de Bloch manualmente
    from qiskit.visualization import plot_bloch_vector
    
    # Extrair coordenadas de Bloch
    x_coord = 2 * np.real(rho_qubit0[0, 1])
    y_coord = 2 * np.imag(rho_qubit0[1, 0])
    z_coord = np.real(rho_qubit0[0, 0] - rho_qubit0[1, 1])
    
    fig = plot_bloch_vector([x_coord, y_coord, z_coord], title="State with Phase Damping")
    path4 = os.path.join(pasta, "04_bloch_com_ruido.png")
    fig.savefig(path4, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"   ✓ {path4}")
except Exception as e:
    print(f"   ✗ Erro: {e}")

if check_timeout():
    sys.exit(0)

# Visualização 5: City 3D com ruído
print("\n5️⃣ Gerando State City 3D (COM ruído)...")
try:
    fig = plot_state_city(rho)
    path5 = os.path.join(pasta, "05_city3d_com_ruido.png")
    fig.savefig(path5, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"   ✓ {path5}")
except Exception as e:
    print(f"   ✗ Erro: {e}")

if check_timeout():
    sys.exit(0)

# Visualização 6: Circuito com mais camadas
print("\n6️⃣ Gerando circuito com entanglement...")
try:
    qc_ent = QuantumCircuit(4)
    
    # Encoding
    for i, val in enumerate([0.5, 0.3, 0.7, 0.4]):
        qc_ent.ry(val * np.pi, i)
    
    # Entanglement layers
    for i in range(3):
        qc_ent.cx(i, i+1)
    qc_ent.cx(3, 0)
    
    # Rotations
    for i in range(4):
        qc_ent.ry(0.2 * (i+1), i)
        qc_ent.rz(0.15 * (i+1), i)
    
    fig = circuit_drawer(qc_ent, output='mpl', style='iqp')
    path6 = os.path.join(pasta, "06_circuito_entangled.png")
    fig.savefig(path6, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"   ✓ {path6}")
except Exception as e:
    print(f"   ✗ Erro: {e}")

if check_timeout():
    sys.exit(0)

# Resumo
tempo_total = time.time() - inicio_global
print("\n" + "=" * 80)
print("RESUMO")
print("=" * 80)

import glob
imagens = sorted(glob.glob(os.path.join(pasta, "*.png")))
print(f"\n✅ {len(imagens)} visualizações geradas em {tempo_total:.1f}s")
print(f"📁 Pasta: {pasta}/\n")

for img in imagens:
    print(f"   • {os.path.basename(img)}")

print(f"\n⏱️ Tempo total: {tempo_total:.1f}s / {TIMEOUT}s")
print("=" * 80)
