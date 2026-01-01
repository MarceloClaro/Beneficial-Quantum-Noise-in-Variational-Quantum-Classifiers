# =============================================================================
# FRAMEWORK DE RUÍDO QUÂNTICO BENÉFICO - VERSÃO CIRQ
# =============================================================================
"""
Framework Investigativo Completo v8.0 - Versão Cirq (Google)

Este arquivo implementa o framework completo de análise de ruído quântico benéfico
usando Cirq (Google) em vez de PennyLane, mantendo a mesma funcionalidade e interface.

Referências:
- Cirq Documentation: https://quantumai.google/cirq
- Google Quantum AI: https://quantumai.google/
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

# Cirq imports
try:
    import cirq
    import sympy
    CIRQ_AVAILABLE = True
except ImportError as e:
    CIRQ_AVAILABLE = False
    print(f"⚠️ Cirq não disponível. Instale com: pip install cirq")
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
# MÓDULO 1: CONSTANTES FUNDAMENTAIS (Compatível com PennyLane/Qiskit)
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
    
    # Constantes para Inicialização
    GLOROT_UNIFORM = lambda n_in, n_out: np.sqrt(6.0 / (n_in + n_out))
    HE_UNIFORM = lambda n_in: np.sqrt(6.0 / n_in)
    LECUN_UNIFORM = lambda n_in: np.sqrt(3.0 / n_in)

# ============================================================================
# MÓDULO 2: MODELOS DE RUÍDO QUÂNTICO PARA CIRQ
# ============================================================================

class ModeloRuidoCirq:
    """Classe base para modelos de ruído quântico no Cirq."""
    def __init__(self, nivel=0.01):
        self.nivel = nivel
    
    def aplicar(self, circuit: cirq.Circuit, qubits: List[cirq.Qid], nivel_override=None):
        """
        Apply quantum noise to circuit qubits.
        
        Args:
            circuit: Cirq circuit to apply noise to
            qubits: List of qubits to apply noise to
            nivel_override: Optional noise level override
        """
        raise NotImplementedError

class RuidoDepolarizanteCirq(ModeloRuidoCirq):
    """Canal de depolarização em Cirq."""
    def aplicar(self, circuit: cirq.Circuit, qubits: List[cirq.Qid], nivel_override=None):
        p = self.nivel if nivel_override is None else nivel_override
        for qubit in qubits:
            circuit.append(cirq.depolarize(p).on(qubit))

class RuidoAmplitudeDampingCirq(ModeloRuidoCirq):
    """Canal de amplitude damping em Cirq."""
    def aplicar(self, circuit: cirq.Circuit, qubits: List[cirq.Qid], nivel_override=None):
        gamma = self.nivel if nivel_override is None else nivel_override
        for qubit in qubits:
            circuit.append(cirq.amplitude_damp(gamma).on(qubit))

class RuidoPhaseDampingCirq(ModeloRuidoCirq):
    """Canal de phase damping em Cirq."""
    def aplicar(self, circuit: cirq.Circuit, qubits: List[cirq.Qid], nivel_override=None):
        gamma = self.nivel if nivel_override is None else nivel_override
        for qubit in qubits:
            circuit.append(cirq.phase_damp(gamma).on(qubit))

class RuidoBitFlipCirq(ModeloRuidoCirq):
    """Canal de bit flip em Cirq."""
    def aplicar(self, circuit: cirq.Circuit, qubits: List[cirq.Qid], nivel_override=None):
        p = self.nivel if nivel_override is None else nivel_override
        for qubit in qubits:
            circuit.append(cirq.bit_flip(p).on(qubit))

class RuidoPhaseFlipCirq(ModeloRuidoCirq):
    """Canal de phase flip em Cirq."""
    def aplicar(self, circuit: cirq.Circuit, qubits: List[cirq.Qid], nivel_override=None):
        p = self.nivel if nivel_override is None else nivel_override
        for qubit in qubits:
            circuit.append(cirq.phase_flip(p).on(qubit))

# Dicionário de modelos de ruído
MODELOS_RUIDO_CIRQ = {
    'sem_ruido': None,
    'depolarizante': RuidoDepolarizanteCirq,
    'amplitude_damping': RuidoAmplitudeDampingCirq,
    'phase_damping': RuidoPhaseDampingCirq,
    'bit_flip': RuidoBitFlipCirq,
    'phase_flip': RuidoPhaseFlipCirq,
}

# ============================================================================
# MÓDULO 3: ANSÄTZE QUÂNTICOS (9 TIPOS)
# ============================================================================

def criar_circuito_basico_cirq(
    n_qubits: int,
    n_camadas: int
) -> Tuple[cirq.Circuit, List[sympy.Symbol], List[cirq.Qid]]:
    """
    Circuito Básico: RY encoding + RY rotations + CNOT ring.
    
    Args:
        n_qubits: Number of qubits
        n_camadas: Number of layers
    
    Returns:
        Tuple of (circuit, parameters, qubits)
    """
    qubits = [cirq.LineQubit(i) for i in range(n_qubits)]
    circuit = cirq.Circuit()
    
    # Parâmetros: features (n_qubits) + variational (n_qubits * n_camadas)
    params = []
    
    # Feature parameters
    x_params = [sympy.Symbol(f'x_{i}') for i in range(n_qubits)]
    params.extend(x_params)
    
    # Encoding layer
    for i, qubit in enumerate(qubits):
        circuit.append(cirq.ry(np.pi * x_params[i])(qubit))
    
    # Variational layers
    for layer in range(n_camadas):
        # Rotation layer
        layer_params = [sympy.Symbol(f'theta_{layer}_{i}') for i in range(n_qubits)]
        params.extend(layer_params)
        
        for i, qubit in enumerate(qubits):
            circuit.append(cirq.ry(layer_params[i])(qubit))
        
        # Entanglement layer (CNOT ring)
        for i in range(n_qubits):
            circuit.append(cirq.CNOT(qubits[i], qubits[(i + 1) % n_qubits]))
    
    return circuit, params, qubits

def criar_circuito_hardware_efficient_cirq(
    n_qubits: int,
    n_camadas: int
) -> Tuple[cirq.Circuit, List[sympy.Symbol], List[cirq.Qid]]:
    """
    Hardware-Efficient Ansatz: RY encoding + (RX, RZ, CNOT) layers.
    
    Referência: Kandala et al. (2017). "Hardware-efficient variational quantum eigensolver"
    """
    qubits = [cirq.LineQubit(i) for i in range(n_qubits)]
    circuit = cirq.Circuit()
    params = []
    
    # Feature encoding
    x_params = [sympy.Symbol(f'x_{i}') for i in range(n_qubits)]
    params.extend(x_params)
    for i, qubit in enumerate(qubits):
        circuit.append(cirq.ry(np.pi * x_params[i])(qubit))
    
    # Variational layers
    for layer in range(n_camadas):
        # RY rotations
        ry_params = [sympy.Symbol(f'ry_{layer}_{i}') for i in range(n_qubits)]
        params.extend(ry_params)
        for i, qubit in enumerate(qubits):
            circuit.append(cirq.ry(ry_params[i])(qubit))
        
        # RZ rotations
        rz_params = [sympy.Symbol(f'rz_{layer}_{i}') for i in range(n_qubits)]
        params.extend(rz_params)
        for i, qubit in enumerate(qubits):
            circuit.append(cirq.rz(rz_params[i])(qubit))
        
        # CNOT entanglement
        for i in range(n_qubits - 1):
            circuit.append(cirq.CNOT(qubits[i], qubits[i + 1]))
    
    return circuit, params, qubits

def criar_circuito_strongly_entangling_cirq(
    n_qubits: int,
    n_camadas: int
) -> Tuple[cirq.Circuit, List[sympy.Symbol], List[cirq.Qid]]:
    """
    Strongly Entangling Ansatz: Dense entanglement pattern.
    
    Referência: Schuld et al. (2020). "Circuit-centric quantum classifiers"
    """
    qubits = [cirq.LineQubit(i) for i in range(n_qubits)]
    circuit = cirq.Circuit()
    params = []
    
    # Feature encoding
    x_params = [sympy.Symbol(f'x_{i}') for i in range(n_qubits)]
    params.extend(x_params)
    for i, qubit in enumerate(qubits):
        circuit.append(cirq.ry(np.pi * x_params[i])(qubit))
    
    # Variational layers
    for layer in range(n_camadas):
        # Three rotation parameters per qubit
        for rot_type, gate in [('rx', cirq.rx), ('ry', cirq.ry), ('rz', cirq.rz)]:
            rot_params = [sympy.Symbol(f'{rot_type}_{layer}_{i}') for i in range(n_qubits)]
            params.extend(rot_params)
            for i, qubit in enumerate(qubits):
                circuit.append(gate(rot_params[i])(qubit))
        
        # Full entanglement (all-to-all CNOT)
        for i in range(n_qubits):
            for j in range(i + 1, n_qubits):
                circuit.append(cirq.CNOT(qubits[i], qubits[j]))
    
    return circuit, params, qubits

def criar_circuito_alternating_layers_cirq(
    n_qubits: int,
    n_camadas: int
) -> Tuple[cirq.Circuit, List[sympy.Symbol], List[cirq.Qid]]:
    """
    Alternating Layers: Alternates between rotation and entanglement layers.
    
    Referência: McClean et al. (2018). "Barren plateaus in quantum neural networks"
    """
    qubits = [cirq.LineQubit(i) for i in range(n_qubits)]
    circuit = cirq.Circuit()
    params = []
    
    # Feature encoding
    x_params = [sympy.Symbol(f'x_{i}') for i in range(n_qubits)]
    params.extend(x_params)
    for i, qubit in enumerate(qubits):
        circuit.append(cirq.ry(np.pi * x_params[i])(qubit))
    
    # Variational layers
    for layer in range(n_camadas):
        # Rotation layer
        ry_params = [sympy.Symbol(f'ry_{layer}_{i}') for i in range(n_qubits)]
        params.extend(ry_params)
        for i, qubit in enumerate(qubits):
            circuit.append(cirq.ry(ry_params[i])(qubit))
        
        # Entanglement layer (even-odd pattern)
        # Even pairs
        for i in range(0, n_qubits - 1, 2):
            circuit.append(cirq.CNOT(qubits[i], qubits[i + 1]))
        # Odd pairs
        for i in range(1, n_qubits - 1, 2):
            circuit.append(cirq.CNOT(qubits[i], qubits[i + 1]))
    
    return circuit, params, qubits

def criar_circuito_brickwork_cirq(
    n_qubits: int,
    n_camadas: int
) -> Tuple[cirq.Circuit, List[sympy.Symbol], List[cirq.Qid]]:
    """
    Brickwork Ansatz: Alternating brick-layer pattern.
    
    Referência: Farhi & Neven (2018). "Classification with quantum neural networks"
    """
    qubits = [cirq.LineQubit(i) for i in range(n_qubits)]
    circuit = cirq.Circuit()
    params = []
    
    # Feature encoding
    x_params = [sympy.Symbol(f'x_{i}') for i in range(n_qubits)]
    params.extend(x_params)
    for i, qubit in enumerate(qubits):
        circuit.append(cirq.ry(np.pi * x_params[i])(qubit))
    
    # Variational layers
    for layer in range(n_camadas):
        # Rotation layer
        ry_params = [sympy.Symbol(f'ry_{layer}_{i}') for i in range(n_qubits)]
        params.extend(ry_params)
        for i, qubit in enumerate(qubits):
            circuit.append(cirq.ry(ry_params[i])(qubit))
        
        # Brickwork entanglement
        if layer % 2 == 0:
            # Even layer: 0-1, 2-3, 4-5, ...
            for i in range(0, n_qubits - 1, 2):
                circuit.append(cirq.CNOT(qubits[i], qubits[i + 1]))
        else:
            # Odd layer: 1-2, 3-4, 5-6, ...
            for i in range(1, n_qubits - 1, 2):
                circuit.append(cirq.CNOT(qubits[i], qubits[i + 1]))
    
    return circuit, params, qubits

def criar_circuito_random_entangling_cirq(
    n_qubits: int,
    n_camadas: int,
    seed: int = 42
) -> Tuple[cirq.Circuit, List[sympy.Symbol], List[cirq.Qid]]:
    """
    Random Entangling: Random CNOT pattern for each layer.
    
    Referência: Haug et al. (2021). "Capacity and quantum geometry of parametrized quantum circuits"
    """
    np.random.seed(seed)
    qubits = [cirq.LineQubit(i) for i in range(n_qubits)]
    circuit = cirq.Circuit()
    params = []
    
    # Feature encoding
    x_params = [sympy.Symbol(f'x_{i}') for i in range(n_qubits)]
    params.extend(x_params)
    for i, qubit in enumerate(qubits):
        circuit.append(cirq.ry(np.pi * x_params[i])(qubit))
    
    # Variational layers
    for layer in range(n_camadas):
        # Rotation layer
        ry_params = [sympy.Symbol(f'ry_{layer}_{i}') for i in range(n_qubits)]
        params.extend(ry_params)
        for i, qubit in enumerate(qubits):
            circuit.append(cirq.ry(ry_params[i])(qubit))
        
        # Random entanglement
        available_qubits = list(range(n_qubits))
        np.random.shuffle(available_qubits)
        for i in range(0, len(available_qubits) - 1, 2):
            circuit.append(cirq.CNOT(qubits[available_qubits[i]], qubits[available_qubits[i + 1]]))
    
    return circuit, params, qubits

def criar_circuito_tree_cirq(
    n_qubits: int,
    n_camadas: int
) -> Tuple[cirq.Circuit, List[sympy.Symbol], List[cirq.Qid]]:
    """
    Tree Tensor Network: Hierarchical binary tree entanglement.
    
    Referência: Grant et al. (2019). "Hierarchical quantum classifiers"
    """
    qubits = [cirq.LineQubit(i) for i in range(n_qubits)]
    circuit = cirq.Circuit()
    params = []
    
    # Feature encoding
    x_params = [sympy.Symbol(f'x_{i}') for i in range(n_qubits)]
    params.extend(x_params)
    for i, qubit in enumerate(qubits):
        circuit.append(cirq.ry(np.pi * x_params[i])(qubit))
    
    # Variational layers
    for layer in range(n_camadas):
        # Rotation layer
        ry_params = [sympy.Symbol(f'ry_{layer}_{i}') for i in range(n_qubits)]
        params.extend(ry_params)
        for i, qubit in enumerate(qubits):
            circuit.append(cirq.ry(ry_params[i])(qubit))
        
        # Tree entanglement (binary tree)
        step = 1
        while step < n_qubits:
            for i in range(0, n_qubits, step * 2):
                if i + step < n_qubits:
                    circuit.append(cirq.CNOT(qubits[i], qubits[i + step]))
            step *= 2
    
    return circuit, params, qubits

def criar_circuito_star_entanglement_cirq(
    n_qubits: int,
    n_camadas: int
) -> Tuple[cirq.Circuit, List[sympy.Symbol], List[cirq.Qid]]:
    """
    Star Entanglement: Central qubit entangles with all others.
    
    Referência: Verdon et al. (2019). "Learning to learn with quantum neural networks"
    """
    qubits = [cirq.LineQubit(i) for i in range(n_qubits)]
    circuit = cirq.Circuit()
    params = []
    
    # Feature encoding
    x_params = [sympy.Symbol(f'x_{i}') for i in range(n_qubits)]
    params.extend(x_params)
    for i, qubit in enumerate(qubits):
        circuit.append(cirq.ry(np.pi * x_params[i])(qubit))
    
    # Variational layers
    for layer in range(n_camadas):
        # Rotation layer
        ry_params = [sympy.Symbol(f'ry_{layer}_{i}') for i in range(n_qubits)]
        params.extend(ry_params)
        for i, qubit in enumerate(qubits):
            circuit.append(cirq.ry(ry_params[i])(qubit))
        
        # Star entanglement (central qubit = 0)
        for i in range(1, n_qubits):
            circuit.append(cirq.CNOT(qubits[0], qubits[i]))
    
    return circuit, params, qubits

def criar_circuito_qaoa_cirq(
    n_qubits: int,
    n_camadas: int
) -> Tuple[cirq.Circuit, List[sympy.Symbol], List[cirq.Qid]]:
    """
    QAOA-inspired Ansatz: Problem + mixer layers.
    
    Referência: Farhi et al. (2014). "A Quantum Approximate Optimization Algorithm"
    """
    qubits = [cirq.LineQubit(i) for i in range(n_qubits)]
    circuit = cirq.Circuit()
    params = []
    
    # Feature encoding
    x_params = [sympy.Symbol(f'x_{i}') for i in range(n_qubits)]
    params.extend(x_params)
    for i, qubit in enumerate(qubits):
        circuit.append(cirq.ry(np.pi * x_params[i])(qubit))
    
    # QAOA layers (problem + mixer)
    for layer in range(n_camadas):
        # Problem layer (ZZ interactions)
        gamma = sympy.Symbol(f'gamma_{layer}')
        params.append(gamma)
        for i in range(n_qubits - 1):
            circuit.append(cirq.ZZ(qubits[i], qubits[i + 1]) ** gamma)
        
        # Mixer layer (RX rotations)
        beta = sympy.Symbol(f'beta_{layer}')
        params.append(beta)
        for i, qubit in enumerate(qubits):
            circuit.append(cirq.rx(2 * beta)(qubit))
    
    return circuit, params, qubits

# Dicionário de ansätze disponíveis
ANSATZE_CIRQ = {
    'basico': criar_circuito_basico_cirq,
    'hardware_efficient': criar_circuito_hardware_efficient_cirq,
    'strongly_entangling': criar_circuito_strongly_entangling_cirq,
    'alternating_layers': criar_circuito_alternating_layers_cirq,
    'brickwork': criar_circuito_brickwork_cirq,
    'random_entangling': criar_circuito_random_entangling_cirq,
    'tree': criar_circuito_tree_cirq,
    'star_entanglement': criar_circuito_star_entanglement_cirq,
    'qaoa': criar_circuito_qaoa_cirq,
}

# ============================================================================
# MÓDULO 4: CLASSIFICADOR VQC CIRQ
# ============================================================================

class ClassificadorVQCCirq(BaseEstimator, ClassifierMixin):
    """
    Variational Quantum Classifier usando Cirq.
    
    Compatível com scikit-learn API.
    
    Parameters:
    -----------
    n_qubits : int
        Número de qubits (deve ser ≥ número de features)
    n_camadas : int
        Número de camadas variacionais
    ansatz : str
        Tipo de ansatz ('basico', 'hardware_efficient', etc.)
    modelo_ruido : str or None
        Tipo de ruído quântico a aplicar
    nivel_ruido : float
        Intensidade do ruído (0.0 a 1.0)
    n_shots : int
        Número de medições por execução
    otimizador : str
        Otimizador a usar ('adam', 'sgd', etc.)
    taxa_aprendizado : float
        Taxa de aprendizado do otimizador
    n_epocas : int
        Número de épocas de treinamento
    verbose : bool
        Se True, exibe progresso
    seed : int
        Seed para reprodutibilidade
    """
    
    def __init__(
        self,
        n_qubits: int = 4,
        n_camadas: int = 2,
        ansatz: str = 'hardware_efficient',
        modelo_ruido: Optional[str] = None,
        nivel_ruido: float = 0.01,
        n_shots: int = 1024,
        otimizador: str = 'adam',
        taxa_aprendizado: float = 0.01,
        n_epocas: int = 100,
        verbose: bool = True,
        seed: int = 42
    ):
        self.n_qubits = n_qubits
        self.n_camadas = n_camadas
        self.ansatz = ansatz
        self.modelo_ruido = modelo_ruido
        self.nivel_ruido = nivel_ruido
        self.n_shots = n_shots
        self.otimizador = otimizador
        self.taxa_aprendizado = taxa_aprendizado
        self.n_epocas = n_epocas
        self.verbose = verbose
        self.seed = seed
        
        self.weights_ = None
        self.classes_ = None
        self.n_features_in_ = None
        self.historico_ = {'custo': [], 'acuracia': []}
        
        # Configurar seed
        np.random.seed(seed)
    
    def _criar_circuito(self):
        """Cria o circuito quântico base."""
        if self.ansatz not in ANSATZE_CIRQ:
            raise ValueError(f"Ansatz '{self.ansatz}' não suportado. Opções: {list(ANSATZE_CIRQ.keys())}")
        
        circuit, params, qubits = ANSATZE_CIRQ[self.ansatz](self.n_qubits, self.n_camadas)
        
        # Adicionar medição
        circuit.append(cirq.measure(qubits[0], key='result'))
        
        return circuit, params, qubits
    
    def _executar_circuito(self, x: np.ndarray, weights: np.ndarray) -> float:
        """
        Executa o circuito com parâmetros dados.
        
        Returns:
            Expectation value ∈ [-1, 1]
        """
        # Criar circuito
        circuit, params, qubits = self._criar_circuito()
        
        # Criar resolver com parâmetros
        param_dict = {}
        
        # Features (primeiros n_features parâmetros)
        n_features = min(len(x), self.n_qubits)
        for i in range(n_features):
            param_dict[params[i]] = float(x[i])
        
        # Preencher features faltantes com zero
        for i in range(n_features, self.n_qubits):
            param_dict[params[i]] = 0.0
        
        # Weights variacionais (restante dos parâmetros)
        for i in range(self.n_qubits, len(params)):
            if i - self.n_qubits < len(weights):
                param_dict[params[i]] = float(weights[i - self.n_qubits])
            else:
                param_dict[params[i]] = 0.0
        
        # Resolver circuito
        resolved_circuit = cirq.resolve_parameters(circuit, param_dict)
        
        # Simular
        simulator = cirq.Simulator()
        result = simulator.run(resolved_circuit, repetitions=self.n_shots)
        
        # Calcular expectation value
        measurements = result.measurements['result']
        prob_1 = np.mean(measurements)
        prob_0 = 1 - prob_1
        
        # Mapear para [-1, 1]: <Z> = prob_0 - prob_1
        expectation = prob_0 - prob_1
        
        return expectation
    
    def _custo(self, weights: np.ndarray, X: np.ndarray, y: np.ndarray) -> float:
        """Função de custo (cross-entropy loss)."""
        custo_total = 0.0
        for xi, yi in zip(X, y):
            predicao = self._executar_circuito(xi, weights)
            # Cross-entropy loss
            prob = (predicao + 1) / 2  # Map [-1,1] to [0,1]
            prob = np.clip(prob, 1e-10, 1 - 1e-10)
            target = (yi + 1) / 2  # Map [-1,1] to [0,1]
            custo_total += -(target * np.log(prob) + (1 - target) * np.log(1 - prob))
        return custo_total / len(X)
    
    def _gradiente_numerico(self, weights: np.ndarray, X: np.ndarray, y: np.ndarray, epsilon: float = 1e-4) -> np.ndarray:
        """Calcula gradiente numérico."""
        grad = np.zeros_like(weights)
        for i in range(len(weights)):
            weights_plus = weights.copy()
            weights_plus[i] += epsilon
            weights_minus = weights.copy()
            weights_minus[i] -= epsilon
            
            custo_plus = self._custo(weights_plus, X, y)
            custo_minus = self._custo(weights_minus, X, y)
            
            grad[i] = (custo_plus - custo_minus) / (2 * epsilon)
        
        return grad
    
    def fit(self, X, y):
        """
        Treina o classificador VQC.
        
        Parameters:
        -----------
        X : array-like of shape (n_samples, n_features)
            Training data
        y : array-like of shape (n_samples,)
            Target values (-1 or +1)
        
        Returns:
        --------
        self : ClassificadorVQCCirq
            Fitted classifier
        """
        if not CIRQ_AVAILABLE:
            raise ImportError("Cirq não está instalado. Instale com: pip install cirq")
        
        X = np.asarray(X)
        y = np.asarray(y)
        
        # Armazenar classes e número de features
        self.classes_ = np.unique(y)
        self.n_features_in_ = X.shape[1]
        
        # Converter labels para -1/+1 se necessário
        if set(self.classes_) != {-1, 1}:
            y = 2 * (y == self.classes_[1]) - 1
        
        # Inicializar pesos
        n_params = self.n_qubits * self.n_camadas * 3  # Estimativa
        self.weights_ = np.random.uniform(-np.pi, np.pi, n_params)
        
        # Treinar com gradiente descendente
        for epoca in range(self.n_epocas):
            # Calcular gradiente
            grad = self._gradiente_numerico(self.weights_, X, y)
            
            # Atualizar pesos
            self.weights_ -= self.taxa_aprendizado * grad
            
            # Registrar histórico
            if epoca % 10 == 0:
                custo = self._custo(self.weights_, X, y)
                y_pred = self.predict(X)
                acuracia = np.mean(y_pred == y)
                
                self.historico_['custo'].append(custo)
                self.historico_['acuracia'].append(acuracia)
                
                if self.verbose:
                    logger.info(f"Época {epoca}/{self.n_epocas} | Custo: {custo:.4f} | Acurácia: {acuracia:.4f}")
        
        return self
    
    def predict(self, X):
        """
        Prediz classes para X.
        
        Parameters:
        -----------
        X : array-like of shape (n_samples, n_features)
            Test data
        
        Returns:
        --------
        y_pred : ndarray of shape (n_samples,)
            Predicted classes
        """
        X = np.asarray(X)
        predicoes = []
        
        for x in X:
            expectation = self._executar_circuito(x, self.weights_)
            # Map expectation to class: >0 → +1, ≤0 → -1
            pred_class = 1 if expectation > 0 else -1
            predicoes.append(pred_class)
        
        return np.array(predicoes)
    
    def get_circuit_diagram(self, output_path: Optional[str] = None) -> str:
        """
        Gera diagrama do circuito.
        
        Parameters:
        -----------
        output_path : str, optional
            Caminho para salvar o diagrama (SVG)
        
        Returns:
        --------
        str : Diagrama em texto ou caminho do arquivo salvo
        """
        circuit, _, _ = self._criar_circuito()
        
        if output_path:
            # Salvar como SVG
            svg = cirq.contrib.svg.circuit_to_svg(circuit)
            with open(output_path, 'w') as f:
                f.write(svg)
            return output_path
        else:
            # Retornar como texto
            return str(circuit)

# ============================================================================
# MÓDULO 5: UTILIDADES
# ============================================================================

def configurar_seeds_cirq(seed: int = 42):
    """Configura seeds para reprodutibilidade."""
    np.random.seed(seed)
    # Cirq usa NumPy internamente

def gerar_dataset_sintetico(n_samples: int = 100, n_features: int = 2, seed: int = 42):
    """Gera dataset sintético para testes."""
    np.random.seed(seed)
    from sklearn.datasets import make_moons
    X, y = make_moons(n_samples=n_samples, noise=0.1, random_state=seed)
    y = 2 * y - 1  # Convert to -1/+1
    return X, y

# ============================================================================
# EXEMPLO DE USO
# ============================================================================

if __name__ == "__main__":
    if not CIRQ_AVAILABLE:
        print("⚠️ Cirq não disponível. Instale com: pip install cirq")
    else:
        print("=" * 80)
        print("FRAMEWORK CIRQ - TESTE BÁSICO")
        print("=" * 80)
        
        # Gerar dados
        X, y = gerar_dataset_sintetico(n_samples=50, n_features=2)
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.3, random_state=42
        )
        
        # Criar classificador
        print("\nCriando classificador VQC com Cirq...")
        clf = ClassificadorVQCCirq(
            n_qubits=4,
            n_camadas=2,
            ansatz='hardware_efficient',
            n_epocas=50,
            verbose=True
        )
        
        # Treinar
        print("\nTreinando...")
        clf.fit(X_train, y_train)
        
        # Avaliar
        y_pred = clf.predict(X_test)
        acuracia = np.mean(y_pred == y_test)
        print(f"\n✓ Acurácia no teste: {acuracia:.2%}")
        
        # Gerar diagrama
        print("\nGerando diagrama do circuito...")
        diagram = clf.get_circuit_diagram()
        print(diagram)
        
        print("\n" + "=" * 80)
        print("✓ TESTE COMPLETO")
        print("=" * 80)
