# =============================================================================
# FRAMEWORK DE RUÍDO QUÂNTICO BENÉFICO - VERSÃO QISKIT
# =============================================================================
"""
Framework Investigativo Completo v7.2 - Versão Qiskit (IBM)

Este arquivo implementa o framework completo de análise de ruído quântico benéfico
usando Qiskit (IBM) em vez de PennyLane, mantendo a mesma funcionalidade e interface.

Referências:
- Qiskit Documentation: https://qiskit.org/documentation/
- IBM Quantum: https://quantum-computing.ibm.com/
- Schuld et al. (2020). "Circuit-centric quantum classifiers." Phys. Rev. A.
"""

import os
import json
import time
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional, Any, List, Tuple

import numpy as np
import pandas as pd

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.model_selection import train_test_split
from sklearn import datasets as sk_datasets
from sklearn.metrics import confusion_matrix
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier

# Estatística
from scipy.stats import f_oneway, ttest_ind
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm

# Visualização
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Qiskit imports
try:
    from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
    from qiskit.circuit import Parameter, ParameterVector
    from qiskit.circuit.library import RealAmplitudes, EfficientSU2, TwoLocal
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, depolarizing_error, amplitude_damping_error, phase_damping_error
    from qiskit.quantum_info import Statevector, DensityMatrix, state_fidelity
    from qiskit.visualization import plot_bloch_multivector, plot_state_city, plot_state_qsphere
    from qiskit.visualization import circuit_drawer
    QISKIT_AVAILABLE = True
except ImportError as e:
    QISKIT_AVAILABLE = False
    print(f"⚠️ Qiskit não disponível. Instale com: pip install qiskit qiskit-aer qiskit-ibm-runtime")
    print(f"   Erro: {e}")

# Optuna para otimização Bayesiana
try:
    import optuna
    from optuna.samplers import TPESampler
    from optuna.pruners import MedianPruner
    OPTUNA_AVAILABLE = True
except ImportError:
    OPTUNA_AVAILABLE = False

# Inicializar logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-8s | %(name)-20s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# ============================================================================
# MÓDULO 1: CONSTANTES FUNDAMENTAIS (Compatível com PennyLane)
# ============================================================================

class ConstantesFundamentais:
    """
    Constantes matemáticas e físicas para inicialização de parâmetros.
    
    Referências:
    - Constantes matemáticas: Weisstein, "MathWorld"
    - Constantes quânticas: CODATA 2018 (Mohr et al., 2019)
    - Normalização: Grant et al. (2019). Quantum.
    """
    
    # Constantes Matemáticas
    PI = np.pi
    E = np.e
    PHI = (1 + np.sqrt(5)) / 2  # Razão Áurea
    SQRT2 = np.sqrt(2)
    LN2 = np.log(2)
    GAMMA = 0.5772156649  # Euler-Mascheroni
    
    # Constantes Quânticas (CODATA 2018)
    HBAR = 1.054571817e-34  # ℏ [J·s]
    ALPHA = 7.2973525693e-3  # α [adimensional]
    RYDBERG = 10973731.568160  # R∞ [m⁻¹]
    
    @classmethod
    def normalizar(cls, valores):
        """Normaliza valores para [-π, π] usando escala logarítmica."""
        log_vals = np.log10(np.abs(valores) + 1e-10)
        norm = (log_vals - log_vals.min()) / (log_vals.max() - log_vals.min() + 1e-10)
        return -np.pi + norm * 2 * np.pi
    
    @classmethod
    def inicializar(cls, n_params, estrategia='aleatorio', seed=42):
        """
        Inicializa parâmetros com diferentes estratégias.
        
        Args:
            n_params: Número de parâmetros
            estrategia: 'matematico', 'quantico', 'fibonacci_spiral', etc.
            seed: Semente aleatória
            
        Returns:
            Array NumPy com parâmetros inicializados
        """
        np.random.seed(seed)
        
        if estrategia == 'matematico':
            const = np.array([cls.PI, cls.E, cls.PHI, cls.SQRT2, cls.LN2, cls.GAMMA])
            n_rep = int(np.ceil(n_params / len(const)))
            params = np.tile(const, n_rep)[:n_params]
            params += np.random.normal(0, 0.1, n_params)
            return cls.normalizar(params)
        
        elif estrategia == 'quantico':
            const = np.array([cls.HBAR, cls.ALPHA, cls.RYDBERG])
            n_rep = int(np.ceil(n_params / len(const)))
            params = np.tile(const, n_rep)[:n_params]
            params += np.random.normal(0, 0.1, n_params)
            return cls.normalizar(params)
        
        elif estrategia == 'fibonacci_spiral':
            phi = (1 + np.sqrt(5)) / 2
            angles = []
            for i in range(n_params):
                angle = (2 * np.pi * (i / phi)) % (2 * np.pi) - np.pi
                angles.append(angle)
            return np.array(angles)
        
        elif estrategia == 'quantum_harmonic':
            hbar_norm = 1.0
            omega = 1.0
            angles = []
            for n in range(n_params):
                E_n = hbar_norm * omega * (n + 0.5)
                angle = (E_n % (2 * np.pi)) - np.pi
                angles.append(angle)
            return np.array(angles)
        
        elif estrategia == 'primes':
            primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
                      59, 61, 67, 71, 73, 79, 83, 89]
            angles = []
            for i in range(n_params):
                p = primes[i % len(primes)]
                angle = (p * np.pi / 10) % (2 * np.pi) - np.pi
                angles.append(angle)
            return np.array(angles)
        
        elif estrategia == 'identity_blocks':
            # Grant et al. (2019) - Identity blocks inicialização
            params = np.random.normal(0, 0.01, n_params)
            return params
        
        else:  # aleatorio (baseline)
            params = np.random.uniform(-np.pi, np.pi, n_params)
            return params


# ============================================================================
# MÓDULO 2: MODELOS DE RUÍDO QUÂNTICO (Qiskit)
# ============================================================================

def criar_modelo_ruido_depolarizante(nivel_ruido: float) -> NoiseModel:
    """Cria modelo de ruído depolarizante usando Qiskit."""
    noise_model = NoiseModel()
    error = depolarizing_error(nivel_ruido, 1)
    noise_model.add_all_qubit_quantum_error(error, ['rx', 'ry', 'rz', 'h'])
    
    # Erro de 2 qubits
    error_2q = depolarizing_error(nivel_ruido * 2, 2)
    noise_model.add_all_qubit_quantum_error(error_2q, ['cx', 'cz'])
    
    return noise_model


def criar_modelo_ruido_amplitude_damping(nivel_ruido: float) -> NoiseModel:
    """Cria modelo de amplitude damping (perda de energia)."""
    noise_model = NoiseModel()
    
    # Erro de 1 qubit
    error_1q = amplitude_damping_error(nivel_ruido)
    noise_model.add_all_qubit_quantum_error(error_1q, ['rx', 'ry', 'rz', 'h'])
    
    # Erro de 2 qubits
    error_2q = amplitude_damping_error(nivel_ruido).tensor(amplitude_damping_error(nivel_ruido))
    noise_model.add_all_qubit_quantum_error(error_2q, ['cx', 'cz'])
    
    return noise_model


def criar_modelo_ruido_phase_damping(nivel_ruido: float) -> NoiseModel:
    """Cria modelo de phase damping (perda de coerência)."""
    noise_model = NoiseModel()
    
    # Erro de 1 qubit
    error_1q = phase_damping_error(nivel_ruido)
    noise_model.add_all_qubit_quantum_error(error_1q, ['rx', 'ry', 'rz', 'h'])
    
    # Erro de 2 qubits (aplicado após a porta de 2 qubits)
    error_2q = phase_damping_error(nivel_ruido).tensor(phase_damping_error(nivel_ruido))
    noise_model.add_all_qubit_quantum_error(error_2q, ['cx', 'cz'])
    
    return noise_model


def criar_modelo_ruido_combinado(nivel_ruido: float) -> NoiseModel:
    """Cria modelo combinando diferentes tipos de ruído."""
    noise_model = NoiseModel()
    
    # Combinar depolarizante + amplitude damping
    error_dep = depolarizing_error(nivel_ruido * 0.5, 1)
    error_amp = amplitude_damping_error(nivel_ruido * 0.5)
    
    noise_model.add_all_qubit_quantum_error(error_dep, ['rx', 'ry', 'rz'])
    noise_model.add_all_qubit_quantum_error(error_amp, ['h', 'cx', 'cz'])
    
    return noise_model


def criar_modelo_ruido_crosstalk(nivel_ruido: float) -> NoiseModel:
    """Cria modelo de crosstalk (interferência entre qubits vizinhos)."""
    noise_model = NoiseModel()
    
    # Erro de 1 qubit (pequeno)
    error_1q = depolarizing_error(nivel_ruido * 0.3, 1)
    noise_model.add_all_qubit_quantum_error(error_1q, ['rx', 'ry', 'rz', 'h'])
    
    # Erro de 2 qubits maior (simula crosstalk)
    error_2q = depolarizing_error(nivel_ruido * 1.5, 2)
    noise_model.add_all_qubit_quantum_error(error_2q, ['cx', 'cz'])
    
    return noise_model


def criar_modelo_ruido_correlacionado(nivel_ruido: float) -> NoiseModel:
    """Cria modelo de ruído correlacionado (afeta qubits de forma correlacionada)."""
    noise_model = NoiseModel()
    
    # Combinar múltiplos tipos de erro de forma correlacionada
    # Depolarizante como base
    error_dep_1q = depolarizing_error(nivel_ruido * 0.6, 1)
    noise_model.add_all_qubit_quantum_error(error_dep_1q, ['rx', 'ry', 'rz'])
    
    # Phase damping correlacionado
    error_phase_1q = phase_damping_error(nivel_ruido * 0.4)
    noise_model.add_all_qubit_quantum_error(error_phase_1q, ['h'])
    
    # 2-qubit errors com correlação
    error_2q_corr = depolarizing_error(nivel_ruido * 0.8, 2)
    noise_model.add_all_qubit_quantum_error(error_2q_corr, ['cx', 'cz'])
    
    return noise_model


# Dicionário de modelos de ruído
MODELOS_RUIDO_QISKIT = {
    'sem_ruido': lambda x: None,
    'depolarizante': criar_modelo_ruido_depolarizante,
    'amplitude_damping': criar_modelo_ruido_amplitude_damping,
    'phase_damping': criar_modelo_ruido_phase_damping,
    'combinado': criar_modelo_ruido_combinado,
    'crosstalk': criar_modelo_ruido_crosstalk,
    'correlacionado': criar_modelo_ruido_correlacionado,
    'correlated_noise': criar_modelo_ruido_correlacionado,  # Alias
}


# ============================================================================
# MÓDULO 3: ARQUITETURAS DE CIRCUITOS QUÂNTICOS (Qiskit)
# ============================================================================

def criar_circuito_basico(n_qubits: int, n_camadas: int) -> Tuple[QuantumCircuit, ParameterVector]:
    """
    Circuito básico com rotações RY e entanglement CNOT em cadeia.
    
    Args:
        n_qubits: Número de qubits
        n_camadas: Número de camadas
        
    Returns:
        Tupla (circuito, parameters)
    """
    n_params = n_qubits * n_camadas * 3  # RX, RY, RZ por qubit por camada
    params = ParameterVector('θ', n_params)
    
    qc = QuantumCircuit(n_qubits)
    
    param_idx = 0
    for camada in range(n_camadas):
        # Rotações
        for q in range(n_qubits):
            qc.rx(params[param_idx], q)
            param_idx += 1
            qc.ry(params[param_idx], q)
            param_idx += 1
            qc.rz(params[param_idx], q)
            param_idx += 1
        
        # Entanglement em cadeia
        if camada < n_camadas - 1:
            for q in range(n_qubits - 1):
                qc.cx(q, q + 1)
    
    return qc, params


def criar_circuito_strongly_entangling(n_qubits: int, n_camadas: int) -> Tuple[QuantumCircuit, ParameterVector]:
    """Circuito fortemente entrelaçado (all-to-all entanglement)."""
    n_params = n_qubits * n_camadas * 3
    params = ParameterVector('θ', n_params)
    
    qc = QuantumCircuit(n_qubits)
    
    param_idx = 0
    for camada in range(n_camadas):
        # Rotações
        for q in range(n_qubits):
            qc.rx(params[param_idx], q)
            param_idx += 1
            qc.ry(params[param_idx], q)
            param_idx += 1
            qc.rz(params[param_idx], q)
            param_idx += 1
        
        # Entanglement all-to-all
        if camada < n_camadas - 1:
            for q1 in range(n_qubits):
                for q2 in range(q1 + 1, n_qubits):
                    qc.cx(q1, q2)
    
    return qc, params


def criar_circuito_hardware_efficient(n_qubits: int, n_camadas: int) -> Tuple[QuantumCircuit, ParameterVector]:
    """Circuito eficiente para hardware (RY + CNOT linear)."""
    n_params = n_qubits * n_camadas
    params = ParameterVector('θ', n_params)
    
    qc = QuantumCircuit(n_qubits)
    
    param_idx = 0
    for camada in range(n_camadas):
        # Apenas RY (mais eficiente)
        for q in range(n_qubits):
            qc.ry(params[param_idx], q)
            param_idx += 1
        
        # CNOT linear
        if camada < n_camadas - 1:
            for q in range(n_qubits - 1):
                qc.cx(q, q + 1)
    
    return qc, params


def criar_circuito_alternating_layers(n_qubits: int, n_camadas: int) -> Tuple[QuantumCircuit, ParameterVector]:
    """Camadas alternadas de rotações e entanglement."""
    n_params = n_qubits * n_camadas * 2
    params = ParameterVector('θ', n_params)
    
    qc = QuantumCircuit(n_qubits)
    
    param_idx = 0
    for camada in range(n_camadas):
        # Camada par: RY + CNOT ímpar
        if camada % 2 == 0:
            for q in range(n_qubits):
                qc.ry(params[param_idx], q)
                param_idx += 1
                qc.rz(params[param_idx], q)
                param_idx += 1
            
            for q in range(0, n_qubits - 1, 2):
                qc.cx(q, q + 1)
        
        # Camada ímpar: RX + CNOT par
        else:
            for q in range(n_qubits):
                qc.rx(params[param_idx], q)
                param_idx += 1
                qc.ry(params[param_idx], q)
                param_idx += 1
            
            for q in range(1, n_qubits - 1, 2):
                qc.cx(q, q + 1)
    
    return qc, params


def criar_circuito_brickwork(n_qubits: int, n_camadas: int) -> Tuple[QuantumCircuit, ParameterVector]:
    """Padrão brickwork (tijolos) de entanglement."""
    n_params = n_qubits * n_camadas * 2
    params = ParameterVector('θ', n_params)
    
    qc = QuantumCircuit(n_qubits)
    
    param_idx = 0
    for camada in range(n_camadas):
        # Rotações
        for q in range(n_qubits):
            qc.ry(params[param_idx], q)
            param_idx += 1
            qc.rz(params[param_idx], q)
            param_idx += 1
        
        # Padrão brickwork alternado
        if camada % 2 == 0:
            for q in range(0, n_qubits - 1, 2):
                qc.cx(q, q + 1)
        else:
            for q in range(1, n_qubits - 1, 2):
                qc.cx(q, q + 1)
    
    return qc, params


def criar_circuito_random_entangling(n_qubits: int, n_camadas: int, seed: int = 42) -> Tuple[QuantumCircuit, ParameterVector]:
    """Circuito com entanglement aleatório."""
    np.random.seed(seed)
    
    n_params = n_qubits * n_camadas * 2
    params = ParameterVector('θ', n_params)
    
    qc = QuantumCircuit(n_qubits)
    
    param_idx = 0
    for camada in range(n_camadas):
        # Rotações
        for q in range(n_qubits):
            qc.ry(params[param_idx], q)
            param_idx += 1
            qc.rz(params[param_idx], q)
            param_idx += 1
        
        # Entanglement aleatório
        if camada < n_camadas - 1:
            n_entangle = np.random.randint(n_qubits // 2, n_qubits)
            pairs = []
            for _ in range(n_entangle):
                q1, q2 = np.random.choice(n_qubits, 2, replace=False)
                if (q1, q2) not in pairs and (q2, q1) not in pairs:
                    qc.cx(q1, q2)
                    pairs.append((q1, q2))
    
    return qc, params


def criar_circuito_tree(n_qubits: int, n_camadas: int) -> Tuple[QuantumCircuit, ParameterVector]:
    """Circuito com topologia de árvore."""
    n_params = n_qubits * n_camadas * 2
    params = ParameterVector('θ', n_params)
    
    qc = QuantumCircuit(n_qubits)
    
    param_idx = 0
    for camada in range(n_camadas):
        # Rotações
        for q in range(n_qubits):
            qc.ry(params[param_idx], q)
            param_idx += 1
            qc.rz(params[param_idx], q)
            param_idx += 1
        
        # Entanglement em árvore (binary tree)
        if camada < n_camadas - 1:
            for level in range(int(np.log2(n_qubits)) + 1):
                step = 2 ** (level + 1)
                for q in range(0, n_qubits, step):
                    if q + 2 ** level < n_qubits:
                        qc.cx(q, q + 2 ** level)
    
    return qc, params


def criar_circuito_star_entanglement(n_qubits: int, n_camadas: int) -> Tuple[QuantumCircuit, ParameterVector]:
    """Circuito com entanglement em estrela (hub central)."""
    n_params = n_qubits * n_camadas * 2
    params = ParameterVector('θ', n_params)
    
    qc = QuantumCircuit(n_qubits)
    
    param_idx = 0
    for camada in range(n_camadas):
        # Rotações
        for q in range(n_qubits):
            qc.ry(params[param_idx], q)
            param_idx += 1
            qc.rz(params[param_idx], q)
            param_idx += 1
        
        # Entanglement em estrela (qubit 0 como hub)
        if camada < n_camadas - 1:
            for q in range(1, n_qubits):
                qc.cx(0, q)
    
    return qc, params


def criar_circuito_qaoa(n_qubits: int, n_camadas: int) -> Tuple[QuantumCircuit, ParameterVector]:
    """Circuito inspirado em QAOA."""
    n_params = n_qubits * n_camadas * 2
    params = ParameterVector('θ', n_params)
    
    qc = QuantumCircuit(n_qubits)
    
    # Hadamard inicial
    for q in range(n_qubits):
        qc.h(q)
    
    param_idx = 0
    for camada in range(n_camadas):
        # Mixing layer (RX)
        for q in range(n_qubits):
            qc.rx(params[param_idx], q)
            param_idx += 1
        
        # Problem layer (RZZ entre pares adjacentes)
        for q in range(n_qubits - 1):
            qc.cx(q, q + 1)
            qc.rz(params[param_idx], q + 1)
            qc.cx(q, q + 1)
            param_idx += 1
        
        # Último parâmetro para fechar a camada
        if param_idx < len(params):
            qc.rz(params[param_idx], 0)
            param_idx += 1
    
    return qc, params


# Dicionário de arquiteturas
ARQUITETURAS_QISKIT = {
    'basico': criar_circuito_basico,
    'basic_entangler': criar_circuito_basico,
    'strongly_entangling': criar_circuito_strongly_entangling,
    'hardware_efficient': criar_circuito_hardware_efficient,
    'alternating_layers': criar_circuito_alternating_layers,
    'brickwork': criar_circuito_brickwork,
    'random_entangling': criar_circuito_random_entangling,
    'tree': criar_circuito_tree,
    'star_entanglement': criar_circuito_star_entanglement,
    'qaoa': criar_circuito_qaoa,
}


# ============================================================================
# MÓDULO 4: CLASSIFICADOR VQC COM QISKIT
# ============================================================================

class ClassificadorVQCQiskit(BaseEstimator, ClassifierMixin):
    """
    Classificador Quântico Variacional usando Qiskit.
    
    Implementa interface scikit-learn (BaseEstimator, ClassifierMixin) para
    compatibilidade com pipelines de ML clássico.
    
    Referências:
    - Schuld et al. (2020). "Circuit-centric quantum classifiers." Phys. Rev. A.
    - Qiskit Machine Learning: https://qiskit.org/documentation/machine-learning/
    """
    
    def __init__(self, n_qubits=4, n_camadas=2, arquitetura='basico',
                 estrategia_init='aleatorio', tipo_ruido='sem_ruido', nivel_ruido=0.01,
                 taxa_aprendizado=0.01, n_epocas=20, batch_size=32, seed=42,
                 shots=1024, otimizador='adam'):
        """
        Args:
            n_qubits: Número de qubits (2-20)
            n_camadas: Profundidade do circuito (1-10)
            arquitetura: Nome da arquitetura do circuito
            estrategia_init: Estratégia de inicialização de parâmetros
            tipo_ruido: Tipo de ruído quântico
            nivel_ruido: Taxa de erro (0.0-0.05)
            taxa_aprendizado: Learning rate
            n_epocas: Número de épocas de treinamento
            batch_size: Tamanho do mini-batch
            seed: Semente aleatória
            shots: Número de medições por execução
            otimizador: 'adam', 'sgd', etc.
        """
        if not QISKIT_AVAILABLE:
            raise ImportError("Qiskit não está disponível. Instale com: pip install qiskit qiskit-aer")
        
        self.n_qubits = n_qubits
        self.n_camadas = n_camadas
        self.arquitetura = arquitetura
        self.estrategia_init = estrategia_init
        self.tipo_ruido = tipo_ruido
        self.nivel_ruido = nivel_ruido
        self.taxa_aprendizado = taxa_aprendizado
        self.n_epocas = n_epocas
        self.batch_size = batch_size
        self.seed = seed
        self.shots = shots
        self.otimizador = otimizador
        
        # Histórico de treinamento
        self.historico_ = {
            'custo': [],
            'acuracia_treino': [],
            'epoca': [],
            'nivel_ruido': []
        }
        
        np.random.seed(seed)
    
    def _criar_circuito(self):
        """Cria o circuito quântico parametrizado."""
        # Criar circuito base
        criar_fn = ARQUITETURAS_QISKIT.get(self.arquitetura, criar_circuito_basico)
        self.circuito_base_, self.params_ = criar_fn(self.n_qubits, self.n_camadas)
        
        # Inicializar pesos
        n_params = len(self.params_)
        self.weights_ = ConstantesFundamentais.inicializar(
            n_params, self.estrategia_init, self.seed
        )
        
        # Bias
        self.bias_ = 0.0
        
        # Criar modelo de ruído
        if self.tipo_ruido != 'sem_ruido':
            criar_ruido_fn = MODELOS_RUIDO_QISKIT.get(self.tipo_ruido)
            if criar_ruido_fn:
                self.noise_model_ = criar_ruido_fn(self.nivel_ruido)
            else:
                self.noise_model_ = None
        else:
            self.noise_model_ = None
        
        # Criar simulador
        if self.noise_model_:
            self.simulador_ = AerSimulator(noise_model=self.noise_model_)
        else:
            self.simulador_ = AerSimulator()
    
    def _codificar_dados(self, qc: QuantumCircuit, x: np.ndarray):
        """
        Codifica dados clássicos no circuito quântico (amplitude encoding).
        
        Args:
            qc: Circuito quântico
            x: Vetor de features
        """
        # Normalizar entrada
        norm = np.linalg.norm(x)
        if norm > 0:
            x_norm = x / norm
        else:
            x_norm = x
        
        # Codificação por rotações RY (ângulos proporcionais aos dados)
        for i in range(min(len(x_norm), self.n_qubits)):
            qc.ry(x_norm[i] * np.pi, i)
    
    def _executar_circuito(self, x: np.ndarray, weights: np.ndarray) -> float:
        """
        Executa o circuito para uma entrada e retorna a expectação.
        
        Args:
            x: Vetor de entrada
            weights: Pesos do circuito
            
        Returns:
            Valor de expectação (predição)
        """
        # Criar circuito completo
        qc = QuantumCircuit(self.n_qubits, 1)
        
        # Codificar dados
        self._codificar_dados(qc, x)
        
        # Adicionar circuito parametrizado
        param_dict = {self.params_[i]: weights[i] for i in range(len(weights))}
        qc_parametrizado = self.circuito_base_.assign_parameters(param_dict)
        qc.compose(qc_parametrizado, inplace=True)
        
        # Medir primeiro qubit
        qc.measure(0, 0)
        
        # Transpilar e executar
        qc_transpiled = transpile(qc, self.simulador_)
        result = self.simulador_.run(qc_transpiled, shots=self.shots).result()
        counts = result.get_counts()
        
        # Calcular expectação: P(0) - P(1)
        p0 = counts.get('0', 0) / self.shots
        p1 = counts.get('1', 0) / self.shots
        expectation = p0 - p1
        
        return expectation
    
    def _gradiente_finito(self, x: np.ndarray, y: float, weights: np.ndarray, epsilon: float = 1e-3) -> np.ndarray:
        """
        Calcula gradiente por diferenças finitas (parameter shift rule).
        
        Args:
            x: Entrada
            y: Label verdadeiro
            weights: Pesos atuais
            epsilon: Tamanho do passo
            
        Returns:
            Gradiente do custo
        """
        grad = np.zeros_like(weights)
        
        for i in range(len(weights)):
            # Shift positivo
            weights_plus = weights.copy()
            weights_plus[i] += epsilon
            pred_plus = self._executar_circuito(x, weights_plus)
            
            # Shift negativo
            weights_minus = weights.copy()
            weights_minus[i] -= epsilon
            pred_minus = self._executar_circuito(x, weights_minus)
            
            # Gradiente do MSE: d/dw (y - pred)^2
            pred_atual = (pred_plus + pred_minus) / 2
            grad[i] = -2 * (y - pred_atual) * (pred_plus - pred_minus) / (2 * epsilon)
        
        return grad
    
    def fit(self, X, y):
        """
        Treina o classificador.
        
        Args:
            X: Dados de treinamento (n_samples, n_features)
            y: Labels (n_samples,)
            
        Returns:
            self
        """
        # Codificar labels como ±1
        self.label_encoder_ = LabelEncoder()
        y_le = self.label_encoder_.fit_transform(y)
        y_encoded = 2 * y_le - 1
        self.classes_ = self.label_encoder_.classes_
        
        # Criar circuito
        self._criar_circuito()
        
        X_arr = np.asarray(X)
        y_arr = np.asarray(y_encoded)
        
        # Treinamento
        logger.info(f"Iniciando treinamento com {self.n_epocas} épocas...")
        
        for epoca in range(self.n_epocas):
            # Embaralhar dados
            indices = np.random.permutation(len(X_arr))
            
            custos_epoca = []
            
            # Mini-batch gradient descent
            for i in range(0, len(X_arr), self.batch_size):
                batch_idx = indices[i:i + self.batch_size]
                X_batch = X_arr[batch_idx]
                y_batch = y_arr[batch_idx]
                
                # Calcular gradiente médio do batch
                grad_total = np.zeros_like(self.weights_)
                custo_batch = 0
                
                for x, y_true in zip(X_batch, y_batch):
                    pred = self._executar_circuito(x, self.weights_)
                    custo_batch += (y_true - pred) ** 2
                    
                    # Gradiente
                    grad = self._gradiente_finito(x, y_true, self.weights_)
                    grad_total += grad
                
                # Média do gradiente
                grad_total /= len(X_batch)
                custo_batch /= len(X_batch)
                custos_epoca.append(custo_batch)
                
                # Atualizar pesos (gradiente descendente)
                if self.otimizador == 'adam':
                    # Simplificação: usar SGD por agora
                    # Adam completo seria mais complexo
                    self.weights_ -= self.taxa_aprendizado * grad_total
                else:
                    self.weights_ -= self.taxa_aprendizado * grad_total
                
                # Atualizar bias
                erro_medio = np.mean([y_true - self._executar_circuito(x, self.weights_) 
                                     for x, y_true in zip(X_batch, y_batch)])
                self.bias_ += self.taxa_aprendizado * erro_medio
            
            # Registrar histórico
            custo_medio = np.mean(custos_epoca)
            acuracia = self.score(X_arr, y)
            
            self.historico_['custo'].append(custo_medio)
            self.historico_['acuracia_treino'].append(acuracia)
            self.historico_['epoca'].append(epoca)
            self.historico_['nivel_ruido'].append(self.nivel_ruido)
            
            if epoca % 5 == 0:
                logger.info(f"Época {epoca}/{self.n_epocas} - Custo: {custo_medio:.4f} - Acurácia: {acuracia:.4f}")
        
        return self
    
    def predict(self, X):
        """Predições binárias."""
        X_arr = np.asarray(X)
        predicoes = []
        
        for x in X_arr:
            pred = self._executar_circuito(x, self.weights_) + self.bias_
            # Converter ±1 de volta para classes
            classe = 1 if pred > 0 else 0
            predicoes.append(classe)
        
        return self.label_encoder_.inverse_transform(predicoes)
    
    def score(self, X, y):
        """Acurácia do modelo."""
        predicoes = self.predict(X)
        return np.mean(predicoes == y)
    
    def get_circuit_diagram(self, output_path: Optional[str] = None) -> str:
        """
        Gera diagrama do circuito.
        
        Args:
            output_path: Caminho para salvar a imagem
            
        Returns:
            Caminho do arquivo salvo ou string do diagrama
        """
        # Criar circuito exemplo
        param_dict = {self.params_[i]: self.weights_[i] for i in range(len(self.weights_))}
        qc_exemplo = self.circuito_base_.assign_parameters(param_dict)
        
        if output_path:
            qc_exemplo.draw('mpl', filename=output_path)
            return output_path
        else:
            return str(qc_exemplo.draw('text'))


# ============================================================================
# MÓDULO 5: VISUALIZAÇÕES AVANÇADAS (Qiskit + Plotly)
# ============================================================================

def visualizar_bloch_sphere(vqc: ClassificadorVQCQiskit, x: np.ndarray, 
                           output_path: Optional[str] = None) -> str:
    """
    Visualiza o estado quântico na esfera de Bloch após codificação.
    
    Args:
        vqc: Classificador VQC treinado
        x: Vetor de entrada para codificar
        output_path: Caminho para salvar a figura
        
    Returns:
        Caminho do arquivo salvo
    """
    # Criar circuito com dados codificados
    qc = QuantumCircuit(vqc.n_qubits)
    vqc._codificar_dados(qc, x)
    
    # Adicionar circuito parametrizado
    param_dict = {vqc.params_[i]: vqc.weights_[i] for i in range(len(vqc.weights_))}
    qc_parametrizado = vqc.circuito_base_.assign_parameters(param_dict)
    qc.compose(qc_parametrizado, inplace=True)
    
    # Obter statevector
    from qiskit.quantum_info import Statevector
    state = Statevector.from_instruction(qc)
    
    # Visualizar
    fig = plot_bloch_multivector(state)
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        return output_path
    else:
        return "bloch_sphere_temp.png"


def visualizar_state_city_3d(vqc: ClassificadorVQCQiskit, x: np.ndarray,
                             output_path: Optional[str] = None) -> str:
    """
    Visualização 3D do estado quântico (city plot).
    
    Args:
        vqc: Classificador VQC treinado
        x: Vetor de entrada
        output_path: Caminho para salvar
        
    Returns:
        Caminho do arquivo salvo
    """
    # Criar circuito
    qc = QuantumCircuit(vqc.n_qubits)
    vqc._codificar_dados(qc, x)
    
    param_dict = {vqc.params_[i]: vqc.weights_[i] for i in range(len(vqc.weights_))}
    qc_parametrizado = vqc.circuito_base_.assign_parameters(param_dict)
    qc.compose(qc_parametrizado, inplace=True)
    
    # Obter density matrix
    from qiskit.quantum_info import DensityMatrix
    state = DensityMatrix.from_instruction(qc)
    
    # Visualizar em 3D
    fig = plot_state_city(state)
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        return output_path
    else:
        return "state_city_3d_temp.png"


def visualizar_qsphere(vqc: ClassificadorVQCQiskit, x: np.ndarray,
                       output_path: Optional[str] = None) -> str:
    """
    Visualização Q-sphere do estado quântico.
    
    Args:
        vqc: Classificador VQC treinado
        x: Vetor de entrada
        output_path: Caminho para salvar
        
    Returns:
        Caminho do arquivo salvo
    """
    # Criar circuito
    qc = QuantumCircuit(vqc.n_qubits)
    vqc._codificar_dados(qc, x)
    
    param_dict = {vqc.params_[i]: vqc.weights_[i] for i in range(len(vqc.weights_))}
    qc_parametrizado = vqc.circuito_base_.assign_parameters(param_dict)
    qc.compose(qc_parametrizado, inplace=True)
    
    # Obter statevector
    from qiskit.quantum_info import Statevector
    state = Statevector.from_instruction(qc)
    
    # Visualizar
    fig = plot_state_qsphere(state)
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        return output_path
    else:
        return "qsphere_temp.png"


# ============================================================================
# MÓDULO 6: UTILIDADES E DATASETS
# ============================================================================

def carregar_datasets():
    """Carrega datasets de benchmark (compatível com framework PennyLane)."""
    datasets = {}
    
    # Iris
    iris = sk_datasets.load_iris()
    X_iris = iris.data[:, :2]  # Primeiras 2 features
    y_iris = (iris.target == 0).astype(int)  # Binário: setosa vs. resto
    X_train, X_test, y_train, y_test = train_test_split(
        X_iris, y_iris, test_size=0.3, random_state=42
    )
    datasets['iris'] = {
        'X_train': X_train, 'X_test': X_test,
        'y_train': y_train, 'y_test': y_test
    }
    
    # Moons
    X_moons, y_moons = sk_datasets.make_moons(n_samples=200, noise=0.1, random_state=42)
    X_train, X_test, y_train, y_test = train_test_split(
        X_moons, y_moons, test_size=0.3, random_state=42
    )
    datasets['moons'] = {
        'X_train': X_train, 'X_test': X_test,
        'y_train': y_train, 'y_test': y_test
    }
    
    # Circles
    X_circles, y_circles = sk_datasets.make_circles(n_samples=200, noise=0.1, factor=0.5, random_state=42)
    X_train, X_test, y_train, y_test = train_test_split(
        X_circles, y_circles, test_size=0.3, random_state=42
    )
    datasets['circles'] = {
        'X_train': X_train, 'X_test': X_test,
        'y_train': y_train, 'y_test': y_test
    }
    
    # Breast Cancer
    breast_cancer = sk_datasets.load_breast_cancer()
    # Usar PCA para reduzir para 4 features (compatível com 4 qubits)
    from sklearn.decomposition import PCA
    pca = PCA(n_components=4)
    X_bc = pca.fit_transform(breast_cancer.data)
    # Normalizar
    X_bc = (X_bc - X_bc.mean(axis=0)) / (X_bc.std(axis=0) + 1e-8)
    y_bc = breast_cancer.target
    X_train, X_test, y_train, y_test = train_test_split(
        X_bc, y_bc, test_size=0.3, random_state=42, stratify=y_bc
    )
    datasets['breast_cancer'] = {
        'X_train': X_train, 'X_test': X_test,
        'y_train': y_train, 'y_test': y_test
    }
    
    # Wine
    wine = sk_datasets.load_wine()
    # Binário: classe 0 vs. resto
    X_wine = wine.data[:, :4]  # Primeiras 4 features
    # Normalizar
    X_wine = (X_wine - X_wine.mean(axis=0)) / (X_wine.std(axis=0) + 1e-8)
    y_wine = (wine.target == 0).astype(int)
    X_train, X_test, y_train, y_test = train_test_split(
        X_wine, y_wine, test_size=0.3, random_state=42, stratify=y_wine
    )
    datasets['wine'] = {
        'X_train': X_train, 'X_test': X_test,
        'y_train': y_train, 'y_test': y_test
    }
    
    return datasets


# ============================================================================
# MÓDULO 7: EXECUÇÃO PRINCIPAL
# ============================================================================

def executar_experimento_qiskit(
    dataset_nome: str = 'moons',
    arquitetura: str = 'basico',
    tipo_ruido: str = 'sem_ruido',
    nivel_ruido: float = 0.01,
    n_epocas: int = 20,
    pasta_resultados: Optional[str] = None,
    verbose: bool = True
) -> Dict[str, Any]:
    """
    Executa um experimento completo com Qiskit VQC.
    
    Args:
        dataset_nome: Nome do dataset
        arquitetura: Arquitetura do circuito
        tipo_ruido: Tipo de ruído quântico
        nivel_ruido: Nível de ruído
        n_epocas: Número de épocas
        pasta_resultados: Pasta para salvar resultados
        verbose: Verbosidade
        
    Returns:
        Dicionário com resultados do experimento
    """
    if verbose:
        logger.info("="*80)
        logger.info("EXECUTANDO EXPERIMENTO QISKIT")
        logger.info("="*80)
        logger.info(f"Dataset: {dataset_nome}")
        logger.info(f"Arquitetura: {arquitetura}")
        logger.info(f"Tipo de ruído: {tipo_ruido}")
        logger.info(f"Nível de ruído: {nivel_ruido}")
    
    # Carregar dataset
    datasets = carregar_datasets()
    dataset = datasets[dataset_nome]
    
    # Criar classificador
    vqc = ClassificadorVQCQiskit(
        n_qubits=4,
        n_camadas=2,
        arquitetura=arquitetura,
        estrategia_init='aleatorio',
        tipo_ruido=tipo_ruido,
        nivel_ruido=nivel_ruido,
        taxa_aprendizado=0.01,
        n_epocas=n_epocas,
        batch_size=32,
        seed=42,
        shots=1024
    )
    
    # Treinar
    inicio = time.time()
    vqc.fit(dataset['X_train'], dataset['y_train'])
    tempo_treino = time.time() - inicio
    
    # Avaliar
    acuracia_treino = vqc.score(dataset['X_train'], dataset['y_train'])
    acuracia_teste = vqc.score(dataset['X_test'], dataset['y_test'])
    
    if verbose:
        logger.info(f"\n✓ Treinamento concluído em {tempo_treino:.2f}s")
        logger.info(f"  Acurácia treino: {acuracia_treino:.4f}")
        logger.info(f"  Acurácia teste: {acuracia_teste:.4f}")
    
    # Gerar visualizações se pasta fornecida
    if pasta_resultados:
        os.makedirs(pasta_resultados, exist_ok=True)
        
        # Diagrama de circuito
        circuit_path = os.path.join(pasta_resultados, 'circuito_qiskit.png')
        vqc.get_circuit_diagram(circuit_path)
        
        # Esfera de Bloch
        x_exemplo = dataset['X_test'][0]
        bloch_path = os.path.join(pasta_resultados, 'bloch_sphere_qiskit.png')
        visualizar_bloch_sphere(vqc, x_exemplo, bloch_path)
        
        # State city 3D
        city_path = os.path.join(pasta_resultados, 'state_city_3d_qiskit.png')
        visualizar_state_city_3d(vqc, x_exemplo, city_path)
        
        if verbose:
            logger.info(f"\n✓ Visualizações salvas em: {pasta_resultados}")
    
    # Resultados
    resultado = {
        'dataset': dataset_nome,
        'arquitetura': arquitetura,
        'tipo_ruido': tipo_ruido,
        'nivel_ruido': nivel_ruido,
        'acuracia_treino': acuracia_treino,
        'acuracia_teste': acuracia_teste,
        'tempo_treino': tempo_treino,
        'historico': vqc.historico_
    }
    
    return resultado


if __name__ == "__main__":
    # Exemplo de uso
    logger.info("Framework Qiskit v7.2 - Ruído Quântico Benéfico")
    logger.info("="*80)
    
    if not QISKIT_AVAILABLE:
        logger.error("Qiskit não está instalado!")
        logger.info("Instale com: pip install qiskit qiskit-aer qiskit-ibm-runtime")
    else:
        # Executar experimento de demonstração
        resultado = executar_experimento_qiskit(
            dataset_nome='moons',
            arquitetura='basico',
            tipo_ruido='depolarizante',
            nivel_ruido=0.01,
            n_epocas=10,
            pasta_resultados='resultados_qiskit_demo',
            verbose=True
        )
        
        logger.info("\n" + "="*80)
        logger.info("EXPERIMENTO CONCLUÍDO")
        logger.info("="*80)
        logger.info(f"Acurácia teste: {resultado['acuracia_teste']:.4f}")
