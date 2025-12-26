"""
Mathematical Validation Module for Quantum Operators.

This module provides functions to validate the mathematical properties
of quantum operators, particularly Kraus operators used in noise models.

Validates compliance with quantum channel axioms:
1. Trace preservation: Σ_k K_k† K_k = I
2. Complete positivity (implicit via Kraus representation)

References:
-----------
Nielsen, M. A., & Chuang, I. L. (2010). Quantum Computation and Quantum Information.
    Cambridge University Press. Chapter 8: Quantum Noise and Quantum Operations.
"""

import numpy as np
from typing import List, Union
import logging

logger = logging.getLogger(__name__)


def validar_operadores_kraus(
    operadores: List[np.ndarray], 
    tol: float = 1e-10,
    verbose: bool = True
) -> bool:
    """
    Valida se os operadores de Kraus satisfazem a condição de completude.
    
    A completude garante que o canal quântico preserva traço:
    
    $$\\sum_{k} K_k^\\dagger K_k = \\mathbb{I}$$
    
    Mathematical Foundation:
    -----------------------
    Para um canal quântico válido ℰ(ρ) = Σ_k K_k ρ K_k†, os operadores de Kraus
    devem satisfazer a condição de completude acima. Isso garante que:
    - Tr(ℰ(ρ)) = Tr(ρ) (preservação de traço)
    - ℰ é um mapa completamente positivo e que preserva traço (CPTP)
    
    Parameters:
    -----------
    operadores : List[np.ndarray]
        Lista de operadores de Kraus (matrizes complexas)
    tol : float, optional
        Tolerância para erro numérico (padrão: 1e-10)
    verbose : bool, optional
        Se True, registra informações detalhadas (padrão: True)
    
    Returns:
    --------
    bool
        True se os operadores satisfazem a condição de completude
    
    Raises:
    -------
    ValueError
        Se os operadores não satisfazem a completude dentro da tolerância
    TypeError
        Se os operadores não são matrizes numpy válidas
    
    Examples:
    ---------
    >>> # Operadores de Kraus para canal de depolarização com p=0.1
    >>> import numpy as np
    >>> p = 0.1
    >>> I = np.eye(2)
    >>> X = np.array([[0, 1], [1, 0]])
    >>> Y = np.array([[0, -1j], [1j, 0]])
    >>> Z = np.array([[1, 0], [0, -1]])
    >>> K0 = np.sqrt(1 - p) * I
    >>> K1 = np.sqrt(p/3) * X
    >>> K2 = np.sqrt(p/3) * Y
    >>> K3 = np.sqrt(p/3) * Z
    >>> validar_operadores_kraus([K0, K1, K2, K3])
    True
    
    References:
    -----------
    Nielsen & Chuang (2010), Section 8.2.3: "Operator-Sum Representation"
    Preskill, J. (2018). "Quantum Computing in the NISQ era and beyond." Quantum, 2, 79.
    """
    # Validar entrada
    if not operadores:
        raise ValueError("Lista de operadores está vazia")
    
    if not all(isinstance(K, np.ndarray) for K in operadores):
        raise TypeError("Todos os operadores devem ser arrays numpy")
    
    # Verificar dimensões consistentes
    dim = operadores[0].shape[0]
    if not all(K.shape == (dim, dim) for K in operadores):
        raise ValueError(
            f"Todos os operadores devem ter a mesma dimensão. "
            f"Encontradas formas: {[K.shape for K in operadores]}"
        )
    
    # Calcular Σ_k K_k† K_k
    soma = sum(np.conj(K).T @ K for K in operadores)
    identidade = np.eye(dim)
    
    # Calcular erro usando norma de Frobenius
    erro = np.linalg.norm(soma - identidade, ord='fro')
    
    if verbose:
        logger.info(f"Validação de Kraus: {len(operadores)} operadores, dimensão {dim}x{dim}")
        logger.info(f"Erro de completude (norma Frobenius): {erro:.2e}")
        logger.info(f"Tolerância: {tol:.2e}")
    
    # Verificar se erro está dentro da tolerância
    if erro > tol:
        erro_msg = (
            f"Operadores de Kraus não satisfazem completude: "
            f"||Σ K_k† K_k - I||_F = {erro:.2e} > {tol:.2e}"
        )
        logger.error(erro_msg)
        raise ValueError(erro_msg)
    
    if verbose:
        logger.info("✓ Operadores de Kraus validados com sucesso")
    
    return True


def obter_operadores_kraus_depolarizante(p: float) -> List[np.ndarray]:
    """
    Retorna os operadores de Kraus para o canal de depolarização.
    
    Mathematical Definition:
    -----------------------
    O canal de depolarização é definido por:
    
    $$\\mathcal{E}(\\rho) = (1-p)\\rho + \\frac{p}{3}(X\\rho X + Y\\rho Y + Z\\rho Z)$$
    
    Operadores de Kraus:
    -------------------
    - K₀ = √(1-p) I
    - K₁ = √(p/3) X
    - K₂ = √(p/3) Y
    - K₃ = √(p/3) Z
    
    Parameters:
    -----------
    p : float
        Probabilidade de erro (0 ≤ p ≤ 1)
    
    Returns:
    --------
    List[np.ndarray]
        Lista de 4 operadores de Kraus (matrizes 2x2 complexas)
    
    Examples:
    ---------
    >>> ops = obter_operadores_kraus_depolarizante(0.1)
    >>> validar_operadores_kraus(ops)
    True
    
    References:
    -----------
    Nielsen & Chuang (2010), Section 8.3.4: "The depolarizing channel"
    """
    if not (0 <= p <= 1):
        raise ValueError(f"Probabilidade p deve estar em [0, 1], recebido: {p}")
    
    # Matrizes de Pauli
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    
    # Operadores de Kraus
    K0 = np.sqrt(1 - p) * I
    K1 = np.sqrt(p / 3) * X
    K2 = np.sqrt(p / 3) * Y
    K3 = np.sqrt(p / 3) * Z
    
    return [K0, K1, K2, K3]


def obter_operadores_kraus_amplitude_damping(gamma: float) -> List[np.ndarray]:
    """
    Retorna os operadores de Kraus para o canal de amplitude damping.
    
    Mathematical Definition:
    -----------------------
    O canal de amplitude damping modela perda de energia (relaxamento T1).
    
    Operadores de Kraus:
    -------------------
    - K₀ = [[1, 0], [0, √(1-γ)]]
    - K₁ = [[0, √γ], [0, 0]]
    
    Parameters:
    -----------
    gamma : float
        Taxa de decaimento (0 ≤ γ ≤ 1)
    
    Returns:
    --------
    List[np.ndarray]
        Lista de 2 operadores de Kraus (matrizes 2x2 complexas)
    
    Examples:
    ---------
    >>> ops = obter_operadores_kraus_amplitude_damping(0.1)
    >>> validar_operadores_kraus(ops)
    True
    
    References:
    -----------
    Nielsen & Chuang (2010), Section 8.3.5: "The amplitude damping channel"
    Clerk et al. (2010). "Introduction to quantum noise." Rev. Mod. Phys., 82, 1155.
    """
    if not (0 <= gamma <= 1):
        raise ValueError(f"Taxa gamma deve estar em [0, 1], recebido: {gamma}")
    
    K0 = np.array([
        [1, 0],
        [0, np.sqrt(1 - gamma)]
    ], dtype=complex)
    
    K1 = np.array([
        [0, np.sqrt(gamma)],
        [0, 0]
    ], dtype=complex)
    
    return [K0, K1]


def obter_operadores_kraus_phase_damping(lmb: float) -> List[np.ndarray]:
    """
    Retorna os operadores de Kraus para o canal de phase damping.
    
    Mathematical Definition:
    -----------------------
    O canal de phase damping modela perda de coerência de fase (decoerência T2).
    
    Operadores de Kraus:
    -------------------
    - K₀ = [[1, 0], [0, √(1-λ)]]
    - K₁ = [[0, 0], [0, √λ]]
    
    Parameters:
    -----------
    lmb : float
        Taxa de perda de fase (0 ≤ λ ≤ 1)
    
    Returns:
    --------
    List[np.ndarray]
        Lista de 2 operadores de Kraus (matrizes 2x2 complexas)
    
    Examples:
    ---------
    >>> ops = obter_operadores_kraus_phase_damping(0.1)
    >>> validar_operadores_kraus(ops)
    True
    
    References:
    -----------
    Nielsen & Chuang (2010), Section 8.3.5: "The phase damping channel"
    Schlosshauer, M. (2007). "Decoherence and the Quantum-to-Classical Transition."
    """
    if not (0 <= lmb <= 1):
        raise ValueError(f"Taxa lambda deve estar em [0, 1], recebido: {lmb}")
    
    K0 = np.array([
        [1, 0],
        [0, np.sqrt(1 - lmb)]
    ], dtype=complex)
    
    K1 = np.array([
        [0, 0],
        [0, np.sqrt(lmb)]
    ], dtype=complex)
    
    return [K0, K1]


if __name__ == "__main__":
    # Testes de validação
    logging.basicConfig(level=logging.INFO)
    
    print("=" * 70)
    print("VALIDAÇÃO DE OPERADORES DE KRAUS")
    print("=" * 70)
    
    # Teste 1: Canal de depolarização
    print("\n1. Canal de Depolarização (p=0.1)")
    print("-" * 70)
    ops_depol = obter_operadores_kraus_depolarizante(0.1)
    validar_operadores_kraus(ops_depol)
    
    # Teste 2: Canal de amplitude damping
    print("\n2. Canal de Amplitude Damping (γ=0.2)")
    print("-" * 70)
    ops_amp = obter_operadores_kraus_amplitude_damping(0.2)
    validar_operadores_kraus(ops_amp)
    
    # Teste 3: Canal de phase damping
    print("\n3. Canal de Phase Damping (λ=0.15)")
    print("-" * 70)
    ops_phase = obter_operadores_kraus_phase_damping(0.15)
    validar_operadores_kraus(ops_phase)
    
    print("\n" + "=" * 70)
    print("✓ TODOS OS TESTES PASSARAM")
    print("=" * 70)
