# =============================================================================
# FRAMEWORK QAOA PARA 100 QUBITS COM ANÁLISE DE RUÍDO BENÉFICO
# =============================================================================
"""
Framework QAOA (Quantum Approximate Optimization Algorithm) para 100 Qubits
com Análise de Ruído Quântico Benéfico

Este módulo implementa QAOA escalável usando Qiskit, mantendo a metodologia
de análise de ruído benéfico do projeto original.

Referências:
- Farhi et al. (2014). "Quantum Approximate Optimization Algorithm." arXiv:1411.4028
- Qiskit Documentation: https://qiskit.org/documentation/
- Zhou et al. (2020). "Quantum approximate optimization algorithm: Performance, mechanism, and implementation on near-term devices." PRX Quantum, 1(2), 020319.

Autor: Framework adaptado para QAOA 100 qubits
Data: 2025-12-26
"""

import os
import json
import time
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional, Any, List, Tuple, Union
from dataclasses import dataclass

import numpy as np
import pandas as pd

from sklearn.base import BaseEstimator
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn import datasets as sk_datasets

# Estatística
from scipy.stats import f_oneway, ttest_ind
from scipy.optimize import minimize

# Visualização
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt

# Qiskit imports
try:
    from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
    from qiskit.circuit import Parameter, ParameterVector
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import (
        NoiseModel, depolarizing_error, amplitude_damping_error, 
        phase_damping_error, thermal_relaxation_error
    )
    from qiskit.quantum_info import Statevector, DensityMatrix, state_fidelity
    QISKIT_AVAILABLE = True
except ImportError as e:
    QISKIT_AVAILABLE = False
    print(f"⚠️ Qiskit não disponível. Instale com: pip install qiskit qiskit-aer")
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
# MÓDULO 1: DATACLASSES PARA CONFIGURAÇÃO
# ============================================================================

@dataclass
class ConfigQAOA:
    """Configuração para experimentos QAOA."""
    n_qubits: int = 100  # Número de qubits
    p_layers: int = 3    # Número de camadas QAOA (p)
    tipo_ruido: str = 'depolarizing'  # Tipo de ruído quântico
    nivel_ruido: float = 0.001  # Nível de ruído (0.0-0.05)
    shots: int = 1024    # Número de medições
    max_iter: int = 100  # Iterações máximas do otimizador
    seed: int = 42       # Semente aleatória
    problema: str = 'maxcut'  # Tipo de problema (maxcut, partition, etc)
    otimizador: str = 'COBYLA'  # Otimizador clássico


@dataclass
class ResultadoQAOA:
    """Resultado de uma execução QAOA."""
    energia_final: float
    parametros_otimos: np.ndarray
    historico_energia: List[float]
    tempo_execucao: float
    probabilidades: Dict[str, float]
    configuracao: ConfigQAOA
    iteracoes: int


# ============================================================================
# MÓDULO 2: CONSTRUTOR DE CIRCUITOS QAOA
# ============================================================================

class ConstrutorCircuitoQAOA:
    """
    Constrói circuitos QAOA para problemas de otimização com 100 qubits.
    
    O QAOA consiste em aplicar alternadamente:
    1. Hamiltoniano do problema (Problem Hamiltonian): U(C, γ) = e^(-iγC)
    2. Hamiltoniano de mixing (Mixer Hamiltonian): U(B, β) = e^(-iβB)
    
    Onde C é o operador de custo do problema e B é o mixing operator.
    """
    
    def __init__(self, n_qubits: int = 100, p_layers: int = 3, seed: int = 42):
        """
        Args:
            n_qubits: Número de qubits (padrão 100)
            p_layers: Número de camadas QAOA (profundidade p)
            seed: Semente aleatória
        """
        if not QISKIT_AVAILABLE:
            raise ImportError("Qiskit é necessário para QAOA")
        
        self.n_qubits = n_qubits
        self.p_layers = p_layers
        self.seed = seed
        np.random.seed(seed)
        
        logger.info(f"Construtor QAOA inicializado: {n_qubits} qubits, {p_layers} camadas")
    
    def criar_grafo_aleatorio(self, densidade: float = 0.5) -> np.ndarray:
        """
        Cria matriz de adjacência para problema MaxCut.
        
        Args:
            densidade: Densidade de conexões (0.0-1.0)
            
        Returns:
            Matriz de adjacência (n_qubits x n_qubits)
        """
        # Para 100 qubits, grafo completo seria muito denso
        # Usar grafo aleatório com densidade controlada
        matriz = np.zeros((self.n_qubits, self.n_qubits))
        
        for i in range(self.n_qubits):
            for j in range(i + 1, self.n_qubits):
                if np.random.random() < densidade:
                    peso = np.random.uniform(0.5, 1.5)  # Pesos aleatórios
                    matriz[i, j] = peso
                    matriz[j, i] = peso
        
        n_arestas = np.sum(matriz > 0) // 2
        logger.info(f"Grafo criado: {n_arestas} arestas (densidade: {densidade:.2%})")
        
        return matriz
    
    def criar_circuito_maxcut(
        self, 
        grafo: np.ndarray,
        gammas: np.ndarray, 
        betas: np.ndarray
    ) -> QuantumCircuit:
        """
        Cria circuito QAOA para problema MaxCut.
        
        MaxCut: Maximizar Σ_{(i,j) ∈ E} w_{ij}(1 - Z_i Z_j)/2
        
        Args:
            grafo: Matriz de adjacência do grafo
            gammas: Parâmetros γ para Hamiltoniano do problema [p_layers]
            betas: Parâmetros β para Hamiltoniano de mixing [p_layers]
            
        Returns:
            QuantumCircuit configurado para MaxCut
        """
        qc = QuantumCircuit(self.n_qubits, self.n_qubits)
        
        # Estado inicial: superposição uniforme |+⟩^⊗n
        qc.h(range(self.n_qubits))
        qc.barrier()
        
        # Aplicar p camadas QAOA
        for p in range(self.p_layers):
            # 1. Hamiltoniano do Problema: U(C, γ)
            # Para MaxCut: aplicar ZZ entre qubits conectados
            for i in range(self.n_qubits):
                for j in range(i + 1, self.n_qubits):
                    if grafo[i, j] > 0:
                        peso = grafo[i, j]
                        # ZZ gate: CNOT - RZ - CNOT
                        qc.cx(i, j)
                        qc.rz(2 * gammas[p] * peso, j)
                        qc.cx(i, j)
            
            qc.barrier()
            
            # 2. Hamiltoniano de Mixing: U(B, β)
            # Mixing padrão: Σ_i X_i
            for i in range(self.n_qubits):
                qc.rx(2 * betas[p], i)
            
            qc.barrier()
        
        # Medição em base computacional
        qc.measure(range(self.n_qubits), range(self.n_qubits))
        
        return qc
    
    def criar_circuito_parametrizado(self) -> Tuple[QuantumCircuit, ParameterVector]:
        """
        Cria circuito QAOA com parâmetros simbólicos (para otimização).
        
        Returns:
            (circuito, vetor_parametros)
        """
        # Criar parâmetros: 2p parâmetros (p gammas + p betas)
        n_params = 2 * self.p_layers
        params = ParameterVector('θ', n_params)
        
        qc = QuantumCircuit(self.n_qubits, self.n_qubits)
        
        # Estado inicial
        qc.h(range(self.n_qubits))
        qc.barrier()
        
        # Estrutura genérica (será customizada com grafo específico)
        for p in range(self.p_layers):
            # Placeholder para Hamiltoniano do problema
            # (será preenchido dinamicamente)
            gamma_idx = p
            
            # Mixing layer
            beta_idx = self.p_layers + p
            for i in range(self.n_qubits):
                qc.rx(2 * params[beta_idx], i)
            
            qc.barrier()
        
        qc.measure(range(self.n_qubits), range(self.n_qubits))
        
        return qc, params


# ============================================================================
# MÓDULO 3: MODELOS DE RUÍDO PARA 100 QUBITS
# ============================================================================

class ModeloRuidoQAOA:
    """
    Modelos de ruído quântico para QAOA em larga escala.
    
    Implementa diferentes canais de Lindblad:
    1. Depolarizing
    2. Amplitude Damping
    3. Phase Damping
    4. Thermal Relaxation
    """
    
    @staticmethod
    def criar_modelo_depolarizing(taxa_erro: float = 0.001) -> NoiseModel:
        """
        Ruído despolarizante: Canal que mistura estado com estado maximamente misto.
        
        ρ → (1-p)ρ + p·I/d
        
        Args:
            taxa_erro: Taxa de erro por porta (0.0-0.05)
            
        Returns:
            NoiseModel do Qiskit
        """
        noise_model = NoiseModel()
        
        # Erro em portas de 1 qubit
        error_1q = depolarizing_error(taxa_erro, 1)
        noise_model.add_all_qubit_quantum_error(error_1q, ['h', 'rx', 'ry', 'rz'])
        
        # Erro em portas de 2 qubits (maior)
        error_2q = depolarizing_error(taxa_erro * 10, 2)
        noise_model.add_all_qubit_quantum_error(error_2q, ['cx', 'cz'])
        
        logger.info(f"Modelo depolarizing criado: taxa_erro={taxa_erro}")
        return noise_model
    
    @staticmethod
    def criar_modelo_amplitude_damping(taxa_erro: float = 0.001) -> NoiseModel:
        """
        Amplitude damping: Simula perda de energia (T1 decay).
        
        Args:
            taxa_erro: Taxa de damping
            
        Returns:
            NoiseModel do Qiskit
        """
        noise_model = NoiseModel()
        
        # Erro em portas de 1 qubit
        error_1q = amplitude_damping_error(taxa_erro)
        noise_model.add_all_qubit_quantum_error(error_1q, ['h', 'rx', 'ry', 'rz'])
        
        # Erro em portas de 2 qubits
        error_2q = amplitude_damping_error(taxa_erro * 10).tensor(
            amplitude_damping_error(taxa_erro * 10)
        )
        noise_model.add_all_qubit_quantum_error(error_2q, ['cx', 'cz'])
        
        logger.info(f"Modelo amplitude damping criado: taxa_erro={taxa_erro}")
        return noise_model
    
    @staticmethod
    def criar_modelo_phase_damping(taxa_erro: float = 0.001) -> NoiseModel:
        """
        Phase damping: Simula perda de coerência (T2 decay).
        
        Args:
            taxa_erro: Taxa de dephasing
            
        Returns:
            NoiseModel do Qiskit
        """
        noise_model = NoiseModel()
        
        # Erro em portas de 1 qubit
        error_1q = phase_damping_error(taxa_erro)
        noise_model.add_all_qubit_quantum_error(error_1q, ['h', 'rx', 'ry', 'rz'])
        
        # Erro em portas de 2 qubits
        error_2q = phase_damping_error(taxa_erro * 10).tensor(
            phase_damping_error(taxa_erro * 10)
        )
        noise_model.add_all_qubit_quantum_error(error_2q, ['cx', 'cz'])
        
        logger.info(f"Modelo phase damping criado: taxa_erro={taxa_erro}")
        return noise_model
    
    @staticmethod
    def criar_modelo_thermal(
        T1: float = 50000.0,  # ns
        T2: float = 70000.0,  # ns
        tempo_porta: float = 100.0  # ns
    ) -> NoiseModel:
        """
        Thermal relaxation: Modelo realista com T1 e T2.
        
        Args:
            T1: Tempo de relaxação de amplitude (ns)
            T2: Tempo de decoerência (ns)
            tempo_porta: Duração da porta quântica (ns)
            
        Returns:
            NoiseModel do Qiskit
        """
        noise_model = NoiseModel()
        
        # Garantir T2 ≤ 2*T1
        T2 = min(T2, 2 * T1)
        
        # Erro em portas de 1 qubit
        error_1q = thermal_relaxation_error(T1, T2, tempo_porta)
        noise_model.add_all_qubit_quantum_error(error_1q, ['h', 'rx', 'ry', 'rz'])
        
        # Erro em portas de 2 qubits (tempo maior)
        error_2q = thermal_relaxation_error(T1, T2, tempo_porta * 2).tensor(
            thermal_relaxation_error(T1, T2, tempo_porta * 2)
        )
        noise_model.add_all_qubit_quantum_error(error_2q, ['cx', 'cz'])
        
        logger.info(f"Modelo thermal criado: T1={T1}ns, T2={T2}ns")
        return noise_model


# Dicionário de modelos disponíveis
MODELOS_RUIDO_QAOA = {
    'depolarizing': ModeloRuidoQAOA.criar_modelo_depolarizing,
    'amplitude_damping': ModeloRuidoQAOA.criar_modelo_amplitude_damping,
    'phase_damping': ModeloRuidoQAOA.criar_modelo_phase_damping,
    'thermal': ModeloRuidoQAOA.criar_modelo_thermal,
    'sem_ruido': lambda *args, **kwargs: None
}


# ============================================================================
# MÓDULO 4: OTIMIZADOR QAOA
# ============================================================================

class OtimizadorQAOA:
    """
    Otimizador para QAOA: gerencia loop quântico-clássico.
    
    Implementa otimização variacional dos parâmetros γ e β para
    minimizar a energia do Hamiltoniano do problema.
    """
    
    def __init__(self, config: ConfigQAOA):
        """
        Args:
            config: Configuração QAOA
        """
        self.config = config
        self.construtor = ConstrutorCircuitoQAOA(
            n_qubits=config.n_qubits,
            p_layers=config.p_layers,
            seed=config.seed
        )
        
        # Criar modelo de ruído
        if config.tipo_ruido in MODELOS_RUIDO_QAOA:
            criar_ruido_fn = MODELOS_RUIDO_QAOA[config.tipo_ruido]
            if config.tipo_ruido == 'thermal':
                self.noise_model = criar_ruido_fn()
            else:
                self.noise_model = criar_ruido_fn(config.nivel_ruido)
        else:
            self.noise_model = None
        
        # Criar simulador
        if self.noise_model:
            self.simulador = AerSimulator(noise_model=self.noise_model)
            logger.info(f"Simulador com ruído {config.tipo_ruido} criado")
        else:
            self.simulador = AerSimulator()
            logger.info("Simulador sem ruído criado")
        
        # Histórico
        self.historico_energia = []
        self.iteracao = 0
    
    def calcular_energia_maxcut(
        self, 
        contagens: Dict[str, int], 
        grafo: np.ndarray
    ) -> float:
        """
        Calcula energia esperada para MaxCut a partir das contagens.
        
        E = Σ_{(i,j)} w_{ij}(1 - ⟨Z_i Z_j⟩)/2
        
        Args:
            contagens: Dicionário {bitstring: count}
            grafo: Matriz de adjacência
            
        Returns:
            Energia esperada (valor a minimizar)
        """
        energia_total = 0.0
        total_shots = sum(contagens.values())
        
        for bitstring, count in contagens.items():
            # Calcular energia para esta configuração
            energia_config = 0.0
            
            # Converter bitstring para array
            config = np.array([int(b) for b in bitstring])
            
            # Somar contribuições de arestas
            for i in range(self.config.n_qubits):
                for j in range(i + 1, self.config.n_qubits):
                    if grafo[i, j] > 0:
                        # Z_i Z_j = (-1)^(s_i + s_j)
                        zi_zj = 1 if config[i] == config[j] else -1
                        energia_config += grafo[i, j] * (1 - zi_zj) / 2
            
            # Média ponderada
            prob = count / total_shots
            energia_total += prob * energia_config
        
        return energia_total
    
    def funcao_objetivo(
        self, 
        params: np.ndarray, 
        grafo: np.ndarray
    ) -> float:
        """
        Função objetivo para otimização: calcula energia esperada.
        
        Args:
            params: Parâmetros [γ_1, ..., γ_p, β_1, ..., β_p]
            grafo: Matriz de adjacência
            
        Returns:
            Energia esperada (negativa para maximização)
        """
        # Separar parâmetros
        gammas = params[:self.config.p_layers]
        betas = params[self.config.p_layers:]
        
        # Criar circuito
        qc = self.construtor.criar_circuito_maxcut(grafo, gammas, betas)
        
        # Executar
        transpiled = transpile(qc, self.simulador)
        job = self.simulador.run(transpiled, shots=self.config.shots)
        result = job.result()
        contagens = result.get_counts()
        
        # Calcular energia
        energia = self.calcular_energia_maxcut(contagens, grafo)
        
        # Armazenar histórico
        self.historico_energia.append(energia)
        self.iteracao += 1
        
        if self.iteracao % 10 == 0:
            logger.info(f"Iteração {self.iteracao}: Energia = {energia:.4f}")
        
        return energia
    
    def otimizar(
        self, 
        grafo: np.ndarray,
        params_iniciais: Optional[np.ndarray] = None
    ) -> ResultadoQAOA:
        """
        Executa otimização QAOA.
        
        Args:
            grafo: Matriz de adjacência do problema
            params_iniciais: Parâmetros iniciais (opcional)
            
        Returns:
            ResultadoQAOA com resultados otimizados
        """
        inicio = time.time()
        
        # Inicializar parâmetros
        if params_iniciais is None:
            # Estratégia padrão: valores pequenos aleatórios
            params_iniciais = np.random.uniform(
                0.0, 0.5, 
                size=2 * self.config.p_layers
            )
        
        logger.info(f"Iniciando otimização QAOA: {self.config.p_layers} camadas, "
                   f"{self.config.n_qubits} qubits")
        logger.info(f"Tipo ruído: {self.config.tipo_ruido}, "
                   f"Nível: {self.config.nivel_ruido}")
        
        # Resetar histórico
        self.historico_energia = []
        self.iteracao = 0
        
        # Otimizar usando scipy
        resultado_opt = minimize(
            fun=lambda p: self.funcao_objetivo(p, grafo),
            x0=params_iniciais,
            method=self.config.otimizador,
            options={'maxiter': self.config.max_iter}
        )
        
        tempo_total = time.time() - inicio
        
        # Executar circuito com parâmetros ótimos para obter probabilidades
        gammas_opt = resultado_opt.x[:self.config.p_layers]
        betas_opt = resultado_opt.x[self.config.p_layers:]
        
        qc_final = self.construtor.criar_circuito_maxcut(grafo, gammas_opt, betas_opt)
        transpiled = transpile(qc_final, self.simulador)
        job = self.simulador.run(transpiled, shots=self.config.shots)
        result = job.result()
        contagens = result.get_counts()
        
        # Normalizar contagens para probabilidades
        total = sum(contagens.values())
        probs = {k: v/total for k, v in contagens.items()}
        
        logger.info(f"Otimização concluída em {tempo_total:.2f}s")
        logger.info(f"Energia final: {resultado_opt.fun:.4f}")
        logger.info(f"Iterações: {resultado_opt.nit}")
        
        return ResultadoQAOA(
            energia_final=resultado_opt.fun,
            parametros_otimos=resultado_opt.x,
            historico_energia=self.historico_energia,
            tempo_execucao=tempo_total,
            probabilidades=probs,
            configuracao=self.config,
            iteracoes=resultado_opt.nit
        )


# ============================================================================
# MÓDULO 5: ANÁLISE DE HIPERPARÂMETROS
# ============================================================================

class AnalisadorHiperparametrosQAOA:
    """
    Busca e análise de hiperparâmetros para QAOA com ruído benéfico.
    
    Investiga:
    - Diferentes níveis de ruído
    - Diferentes profundidades (p-layers)
    - Diferentes tipos de ruído
    """
    
    def __init__(self, pasta_resultados: str = 'resultados_qaoa_100qubits'):
        """
        Args:
            pasta_resultados: Diretório para salvar resultados
        """
        self.pasta_resultados = Path(pasta_resultados)
        self.pasta_resultados.mkdir(parents=True, exist_ok=True)
        
        self.resultados_experimentos = []
        
        logger.info(f"Analisador inicializado: {self.pasta_resultados}")
    
    def grid_search_ruido(
        self,
        grafo: np.ndarray,
        niveis_ruido: List[float],
        tipos_ruido: List[str],
        p_layers: int = 3,
        n_repeticoes: int = 3
    ) -> pd.DataFrame:
        """
        Grid search sobre níveis e tipos de ruído.
        
        Args:
            grafo: Matriz de adjacência do problema
            niveis_ruido: Lista de níveis de ruído a testar
            tipos_ruido: Lista de tipos de ruído
            p_layers: Profundidade QAOA
            n_repeticoes: Repetições por configuração
            
        Returns:
            DataFrame com resultados
        """
        resultados = []
        
        total_exp = len(niveis_ruido) * len(tipos_ruido) * n_repeticoes
        logger.info(f"Iniciando grid search: {total_exp} experimentos")
        
        exp_count = 0
        
        for tipo_ruido in tipos_ruido:
            for nivel_ruido in niveis_ruido:
                for rep in range(n_repeticoes):
                    exp_count += 1
                    
                    # Configuração
                    config = ConfigQAOA(
                        n_qubits=grafo.shape[0],
                        p_layers=p_layers,
                        tipo_ruido=tipo_ruido,
                        nivel_ruido=nivel_ruido,
                        seed=42 + rep
                    )
                    
                    logger.info(f"\n[{exp_count}/{total_exp}] Tipo: {tipo_ruido}, "
                              f"Nível: {nivel_ruido:.4f}, Rep: {rep+1}")
                    
                    # Executar QAOA
                    otimizador = OtimizadorQAOA(config)
                    resultado = otimizador.otimizar(grafo)
                    
                    # Armazenar
                    resultados.append({
                        'tipo_ruido': tipo_ruido,
                        'nivel_ruido': nivel_ruido,
                        'p_layers': p_layers,
                        'repeticao': rep,
                        'energia_final': resultado.energia_final,
                        'tempo_execucao': resultado.tempo_execucao,
                        'iteracoes': resultado.iteracoes,
                        'convergiu': resultado.iteracoes < config.max_iter
                    })
        
        df = pd.DataFrame(resultados)
        
        # Salvar
        arquivo = self.pasta_resultados / f'grid_search_ruido_p{p_layers}.csv'
        df.to_csv(arquivo, index=False)
        logger.info(f"Resultados salvos: {arquivo}")
        
        return df
    
    def otimizacao_bayesiana(
        self,
        grafo: np.ndarray,
        n_trials: int = 50
    ) -> Dict[str, Any]:
        """
        Otimização Bayesiana de hiperparâmetros usando Optuna.
        
        Args:
            grafo: Matriz de adjacência
            n_trials: Número de trials
            
        Returns:
            Dicionário com melhores parâmetros
        """
        if not OPTUNA_AVAILABLE:
            logger.warning("Optuna não disponível. Pulando otimização bayesiana.")
            return {}
        
        def objetivo_optuna(trial):
            # Sugerir hiperparâmetros
            tipo_ruido = trial.suggest_categorical(
                'tipo_ruido',
                ['depolarizing', 'amplitude_damping', 'phase_damping', 'sem_ruido']
            )
            
            if tipo_ruido != 'sem_ruido':
                nivel_ruido = trial.suggest_float('nivel_ruido', 0.0001, 0.01, log=True)
            else:
                nivel_ruido = 0.0
            
            p_layers = trial.suggest_int('p_layers', 1, 5)
            
            # Configuração
            config = ConfigQAOA(
                n_qubits=grafo.shape[0],
                p_layers=p_layers,
                tipo_ruido=tipo_ruido,
                nivel_ruido=nivel_ruido,
                max_iter=50  # Reduzir para trials mais rápidos
            )
            
            # Executar
            otimizador = OtimizadorQAOA(config)
            resultado = otimizador.otimizar(grafo)
            
            # Retornar energia (minimizar)
            return resultado.energia_final
        
        # Criar estudo
        study = optuna.create_study(
            direction='minimize',
            sampler=TPESampler(seed=42),
            pruner=MedianPruner()
        )
        
        logger.info(f"Iniciando otimização Bayesiana: {n_trials} trials")
        study.optimize(objetivo_optuna, n_trials=n_trials)
        
        # Melhores parâmetros
        melhores = study.best_params
        logger.info(f"Melhores hiperparâmetros: {melhores}")
        logger.info(f"Melhor energia: {study.best_value:.4f}")
        
        # Salvar resultados
        df_trials = study.trials_dataframe()
        arquivo = self.pasta_resultados / 'otimizacao_bayesiana.csv'
        df_trials.to_csv(arquivo, index=False)
        
        return {
            'best_params': melhores,
            'best_value': study.best_value,
            'study': study
        }


# ============================================================================
# MÓDULO 6: VISUALIZAÇÕES
# ============================================================================

class VisualizadorQAOA:
    """Visualizações para resultados QAOA."""
    
    @staticmethod
    def plotar_convergencia(resultado: ResultadoQAOA, salvar: Optional[str] = None):
        """
        Plota convergência da energia durante otimização.
        
        Args:
            resultado: ResultadoQAOA
            salvar: Caminho para salvar figura (opcional)
        """
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(
            y=resultado.historico_energia,
            mode='lines+markers',
            name='Energia',
            line=dict(color='blue', width=2),
            marker=dict(size=4)
        ))
        
        fig.update_layout(
            title=f'Convergência QAOA - {resultado.configuracao.tipo_ruido} '
                  f'(nível: {resultado.configuracao.nivel_ruido:.4f})',
            xaxis_title='Iteração',
            yaxis_title='Energia',
            template='plotly_white',
            width=800,
            height=500
        )
        
        if salvar:
            fig.write_html(salvar)
            logger.info(f"Gráfico salvo: {salvar}")
        else:
            fig.show()
    
    @staticmethod
    def plotar_comparacao_ruido(
        df: pd.DataFrame,
        salvar: Optional[str] = None
    ):
        """
        Plota comparação de diferentes tipos e níveis de ruído.
        
        Args:
            df: DataFrame com resultados
            salvar: Caminho para salvar
        """
        fig = go.Figure()
        
        for tipo in df['tipo_ruido'].unique():
            df_tipo = df[df['tipo_ruido'] == tipo]
            
            # Agrupar por nível de ruído
            grouped = df_tipo.groupby('nivel_ruido').agg({
                'energia_final': ['mean', 'std']
            }).reset_index()
            
            fig.add_trace(go.Scatter(
                x=grouped['nivel_ruido'],
                y=grouped['energia_final']['mean'],
                error_y=dict(
                    type='data',
                    array=grouped['energia_final']['std'],
                    visible=True
                ),
                mode='lines+markers',
                name=tipo,
                line=dict(width=2),
                marker=dict(size=8)
            ))
        
        fig.update_layout(
            title='Impacto do Ruído na Energia QAOA (100 Qubits)',
            xaxis_title='Nível de Ruído',
            yaxis_title='Energia Final (média ± std)',
            xaxis_type='log',
            template='plotly_white',
            width=900,
            height=600,
            legend=dict(x=0.02, y=0.98)
        )
        
        if salvar:
            fig.write_html(salvar)
            logger.info(f"Comparação salva: {salvar}")
        else:
            fig.show()


# ============================================================================
# MÓDULO 7: FUNÇÕES DE UTILIDADE
# ============================================================================

def demo_qaoa_100qubits(
    densidade_grafo: float = 0.1,
    p_layers: int = 3,
    tipo_ruido: str = 'depolarizing',
    nivel_ruido: float = 0.001
) -> ResultadoQAOA:
    """
    Demonstração rápida de QAOA com 100 qubits.
    
    Args:
        densidade_grafo: Densidade de conexões do grafo
        p_layers: Profundidade QAOA
        tipo_ruido: Tipo de ruído quântico
        nivel_ruido: Nível de ruído
        
    Returns:
        ResultadoQAOA
    """
    logger.info("="*80)
    logger.info("DEMONSTRAÇÃO QAOA 100 QUBITS")
    logger.info("="*80)
    
    # Criar grafo
    construtor = ConstrutorCircuitoQAOA(n_qubits=100, p_layers=p_layers)
    grafo = construtor.criar_grafo_aleatorio(densidade=densidade_grafo)
    
    # Configurar QAOA
    config = ConfigQAOA(
        n_qubits=100,
        p_layers=p_layers,
        tipo_ruido=tipo_ruido,
        nivel_ruido=nivel_ruido,
        max_iter=100
    )
    
    # Executar
    otimizador = OtimizadorQAOA(config)
    resultado = otimizador.otimizar(grafo)
    
    # Visualizar
    visualizador = VisualizadorQAOA()
    visualizador.plotar_convergencia(resultado)
    
    logger.info("="*80)
    logger.info("DEMONSTRAÇÃO CONCLUÍDA")
    logger.info(f"Energia final: {resultado.energia_final:.4f}")
    logger.info(f"Tempo: {resultado.tempo_execucao:.2f}s")
    logger.info("="*80)
    
    return resultado


def experimento_completo_ruido_benefico(
    n_qubits: int = 100,
    densidade_grafo: float = 0.1,
    p_layers: int = 3
) -> Dict[str, Any]:
    """
    Experimento completo de análise de ruído benéfico em QAOA.
    
    Args:
        n_qubits: Número de qubits
        densidade_grafo: Densidade do grafo MaxCut
        p_layers: Profundidade QAOA
        
    Returns:
        Dicionário com resultados completos
    """
    logger.info("="*80)
    logger.info(f"EXPERIMENTO COMPLETO: RUÍDO BENÉFICO EM QAOA ({n_qubits} QUBITS)")
    logger.info("="*80)
    
    # Criar grafo
    construtor = ConstrutorCircuitoQAOA(n_qubits=n_qubits, p_layers=p_layers)
    grafo = construtor.criar_grafo_aleatorio(densidade=densidade_grafo)
    
    # Analisador
    analisador = AnalisadorHiperparametrosQAOA()
    
    # Grid search
    niveis_ruido = [0.0, 0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01]
    tipos_ruido = ['sem_ruido', 'depolarizing', 'amplitude_damping', 'phase_damping']
    
    df_grid = analisador.grid_search_ruido(
        grafo=grafo,
        niveis_ruido=niveis_ruido,
        tipos_ruido=tipos_ruido,
        p_layers=p_layers,
        n_repeticoes=3
    )
    
    # Visualizar
    visualizador = VisualizadorQAOA()
    visualizador.plotar_comparacao_ruido(
        df_grid,
        salvar=str(analisador.pasta_resultados / 'comparacao_ruido.html')
    )
    
    # Otimização Bayesiana (opcional)
    if OPTUNA_AVAILABLE:
        resultado_bayes = analisador.otimizacao_bayesiana(grafo, n_trials=30)
    else:
        resultado_bayes = {}
    
    logger.info("="*80)
    logger.info("EXPERIMENTO CONCLUÍDO")
    logger.info("="*80)
    
    return {
        'grid_search': df_grid,
        'bayesian_opt': resultado_bayes,
        'grafo': grafo
    }


# ============================================================================
# MAIN: EXEMPLO DE USO
# ============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Framework QAOA 100 Qubits')
    parser.add_argument('--modo', type=str, default='demo',
                       choices=['demo', 'completo', 'grid'],
                       help='Modo de execução')
    parser.add_argument('--n_qubits', type=int, default=100,
                       help='Número de qubits')
    parser.add_argument('--p_layers', type=int, default=3,
                       help='Profundidade QAOA')
    parser.add_argument('--densidade', type=float, default=0.1,
                       help='Densidade do grafo')
    parser.add_argument('--tipo_ruido', type=str, default='depolarizing',
                       help='Tipo de ruído')
    parser.add_argument('--nivel_ruido', type=float, default=0.001,
                       help='Nível de ruído')
    
    args = parser.parse_args()
    
    if not QISKIT_AVAILABLE:
        logger.error("Qiskit não está disponível. Instale com:")
        logger.error("pip install qiskit qiskit-aer")
        exit(1)
    
    if args.modo == 'demo':
        demo_qaoa_100qubits(
            densidade_grafo=args.densidade,
            p_layers=args.p_layers,
            tipo_ruido=args.tipo_ruido,
            nivel_ruido=args.nivel_ruido
        )
    
    elif args.modo == 'completo':
        experimento_completo_ruido_benefico(
            n_qubits=args.n_qubits,
            densidade_grafo=args.densidade,
            p_layers=args.p_layers
        )
    
    elif args.modo == 'grid':
        # Grid search customizado
        construtor = ConstrutorCircuitoQAOA(
            n_qubits=args.n_qubits,
            p_layers=args.p_layers
        )
        grafo = construtor.criar_grafo_aleatorio(densidade=args.densidade)
        
        analisador = AnalisadorHiperparametrosQAOA()
        df = analisador.grid_search_ruido(
            grafo=grafo,
            niveis_ruido=[0.0, 0.001, 0.005, 0.01],
            tipos_ruido=['sem_ruido', 'depolarizing', 'phase_damping'],
            p_layers=args.p_layers,
            n_repeticoes=3
        )
        
        print("\nResultados:")
        print(df.groupby(['tipo_ruido', 'nivel_ruido'])['energia_final'].mean())
