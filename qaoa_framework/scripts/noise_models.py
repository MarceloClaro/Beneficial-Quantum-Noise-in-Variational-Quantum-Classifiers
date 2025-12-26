"""
Tarefa 1: Modelos de Ruído para Qiskit
Migração dos modelos de ruído do VQC para QAOA usando Qiskit.

Referências:
- Qiskit Aer Noise: https://qiskit.org/documentation/apidoc/aer_noise.html
- Nielsen & Chuang (2010). Quantum Computation and Quantum Information.
"""

from typing import Optional, Dict, Any
import numpy as np

try:
    from qiskit_aer.noise import (
        NoiseModel,
        depolarizing_error,
        amplitude_damping_error,
        phase_damping_error,
        thermal_relaxation_error,
        pauli_error
    )
    QISKIT_AVAILABLE = True
except ImportError:
    QISKIT_AVAILABLE = False
    print("⚠️  Qiskit Aer não disponível")


def criar_noise_model_qiskit(
    tipo_ruido: str,
    nivel_ruido: float = 0.01,
    params: Optional[Dict[str, Any]] = None
) -> Optional[NoiseModel]:
    """
    Cria modelo de ruído do Qiskit para QAOA.
    
    Args:
        tipo_ruido: Tipo de ruído ('depolarizing', 'amplitude_damping', 
                    'phase_damping', 'thermal', 'pauli', 'sem_ruido')
        nivel_ruido: Taxa de erro (0.0-0.05)
        params: Parâmetros adicionais (ex: T1, T2 para thermal)
    
    Returns:
        NoiseModel do Qiskit ou None se sem_ruido
        
    Exemplo:
        >>> noise_model = criar_noise_model_qiskit('depolarizing', 0.01)
        >>> # Usar em AerSimulator(noise_model=noise_model)
    """
    if not QISKIT_AVAILABLE:
        raise ImportError("Qiskit Aer não está disponível")
    
    if tipo_ruido == 'sem_ruido' or nivel_ruido == 0.0:
        return None
    
    params = params or {}
    noise_model = NoiseModel()
    
    if tipo_ruido == 'depolarizing':
        # Canal despolarizante: ρ → (1-p)ρ + p·I/d
        # Aplica erro em portas de 1 e 2 qubits
        
        # Portas de 1 qubit
        error_1q = depolarizing_error(nivel_ruido, 1)
        noise_model.add_all_qubit_quantum_error(
            error_1q, 
            ['h', 'rx', 'ry', 'rz', 's', 't', 'x', 'y', 'z']
        )
        
        # Portas de 2 qubits (erro maior)
        error_2q = depolarizing_error(nivel_ruido * 10, 2)
        noise_model.add_all_qubit_quantum_error(
            error_2q,
            ['cx', 'cz', 'cy', 'swap', 'rzz']
        )
    
    elif tipo_ruido == 'amplitude_damping':
        # Amplitude damping: Simula perda de energia (T1 decay)
        # |1⟩ → |0⟩ com probabilidade p
        
        error_1q = amplitude_damping_error(nivel_ruido)
        noise_model.add_all_qubit_quantum_error(
            error_1q,
            ['h', 'rx', 'ry', 'rz', 's', 't', 'x', 'y', 'z']
        )
        
        # 2-qubit: produto tensorial
        error_2q = amplitude_damping_error(nivel_ruido * 10).tensor(
            amplitude_damping_error(nivel_ruido * 10)
        )
        noise_model.add_all_qubit_quantum_error(
            error_2q,
            ['cx', 'cz', 'cy', 'swap', 'rzz']
        )
    
    elif tipo_ruido == 'phase_damping':
        # Phase damping: Simula perda de coerência (T2 decay)
        # Preserva população, mas destrói coerência
        
        error_1q = phase_damping_error(nivel_ruido)
        noise_model.add_all_qubit_quantum_error(
            error_1q,
            ['h', 'rx', 'ry', 'rz', 's', 't', 'x', 'y', 'z']
        )
        
        error_2q = phase_damping_error(nivel_ruido * 10).tensor(
            phase_damping_error(nivel_ruido * 10)
        )
        noise_model.add_all_qubit_quantum_error(
            error_2q,
            ['cx', 'cz', 'cy', 'swap', 'rzz']
        )
    
    elif tipo_ruido == 'thermal':
        # Thermal relaxation: Modelo realista com T1 e T2
        # T1: amplitude decay, T2: phase decay (T2 ≤ 2·T1)
        
        T1 = params.get('T1', 50000.0)  # ns
        T2 = params.get('T2', 70000.0)  # ns
        T2 = min(T2, 2 * T1)  # Garantir T2 ≤ 2·T1
        
        gate_time_1q = params.get('gate_time_1q', 100.0)  # ns
        gate_time_2q = params.get('gate_time_2q', 200.0)  # ns
        
        error_1q = thermal_relaxation_error(T1, T2, gate_time_1q)
        noise_model.add_all_qubit_quantum_error(
            error_1q,
            ['h', 'rx', 'ry', 'rz', 's', 't', 'x', 'y', 'z']
        )
        
        error_2q = thermal_relaxation_error(T1, T2, gate_time_2q).tensor(
            thermal_relaxation_error(T1, T2, gate_time_2q)
        )
        noise_model.add_all_qubit_quantum_error(
            error_2q,
            ['cx', 'cz', 'cy', 'swap', 'rzz']
        )
    
    elif tipo_ruido == 'pauli':
        # Pauli error channel: Combinação de X, Y, Z errors
        px = params.get('px', nivel_ruido / 3)
        py = params.get('py', nivel_ruido / 3)
        pz = params.get('pz', nivel_ruido / 3)
        pi = 1 - (px + py + pz)
        
        error_1q = pauli_error([
            ('X', px),
            ('Y', py),
            ('Z', pz),
            ('I', pi)
        ])
        noise_model.add_all_qubit_quantum_error(
            error_1q,
            ['h', 'rx', 'ry', 'rz', 's', 't', 'x', 'y', 'z']
        )
        
        # 2-qubit: produto tensorial
        error_2q = error_1q.tensor(error_1q)
        noise_model.add_all_qubit_quantum_error(
            error_2q,
            ['cx', 'cz', 'cy', 'swap', 'rzz']
        )
    
    else:
        raise ValueError(f"Tipo de ruído '{tipo_ruido}' não reconhecido. "
                        f"Opções: depolarizing, amplitude_damping, phase_damping, "
                        f"thermal, pauli, sem_ruido")
    
    return noise_model


def aplicar_schedule_ruido(
    nivel_base: float,
    schedule: str,
    iteracao: int,
    max_iter: int
) -> float:
    """
    Aplica schedule dinâmico ao nível de ruído.
    
    Args:
        nivel_base: Nível de ruído inicial
        schedule: Tipo de schedule ('constant', 'linear', 'exponential')
        iteracao: Iteração atual
        max_iter: Máximo de iterações
    
    Returns:
        Nível de ruído ajustado
        
    Schedules:
        - constant: p(t) = p₀
        - linear: p(t) = p₀ · (1 - t/T)
        - exponential: p(t) = p₀ · exp(-λt/T)
    """
    if schedule == 'constant':
        return nivel_base
    
    elif schedule == 'linear':
        # Decay linear
        fator = 1 - (iteracao / max_iter)
        return nivel_base * max(fator, 0.0)
    
    elif schedule == 'exponential':
        # Decay exponencial (λ=3)
        lambda_decay = 3.0
        fator = np.exp(-lambda_decay * iteracao / max_iter)
        return nivel_base * fator
    
    else:
        raise ValueError(f"Schedule '{schedule}' não reconhecido. "
                        f"Opções: constant, linear, exponential")


# Dicionário de modelos disponíveis
MODELOS_RUIDO = {
    'depolarizing': lambda p, params=None: criar_noise_model_qiskit('depolarizing', p, params),
    'amplitude_damping': lambda p, params=None: criar_noise_model_qiskit('amplitude_damping', p, params),
    'phase_damping': lambda p, params=None: criar_noise_model_qiskit('phase_damping', p, params),
    'thermal': lambda p, params=None: criar_noise_model_qiskit('thermal', p, params),
    'pauli': lambda p, params=None: criar_noise_model_qiskit('pauli', p, params),
    'sem_ruido': lambda p, params=None: None
}


if __name__ == "__main__":
    # Teste dos modelos de ruído
    print("Testando modelos de ruído...")
    
    if QISKIT_AVAILABLE:
        for tipo in ['depolarizing', 'amplitude_damping', 'phase_damping', 'thermal']:
            print(f"\n{tipo}:")
            noise_model = criar_noise_model_qiskit(tipo, 0.01)
            if noise_model:
                print(f"  ✓ Modelo criado: {len(noise_model.to_dict()['errors'])} erros definidos")
            else:
                print(f"  ✓ Sem ruído")
        
        # Teste de schedule
        print("\n\nTestando schedules de ruído:")
        for schedule in ['constant', 'linear', 'exponential']:
            niveis = [aplicar_schedule_ruido(0.01, schedule, i, 10) for i in range(11)]
            print(f"  {schedule}: {niveis[0]:.4f} → {niveis[-1]:.4f}")
    
    else:
        print("❌ Qiskit Aer não disponível para testes")
