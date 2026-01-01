"""
Noise-Tailored Ansatz (v10.1-A1)
- Extração T1/T2 (qiskit-ibm-runtime quando disponível)
- Roteamento curto em qubits com maior T1/T2
- Dynamical Decoupling automático (XYXY) em qubits ociosos
- Auditoria: coleta de métricas e hash SHA-256 do layout
"""
from __future__ import annotations
import hashlib
import logging
from typing import Dict, List, Tuple
import numpy as np

logger = logging.getLogger("V10-NoiseAnsatz")

DEFAULT_NOISE = {
    0: {"T1": 50e-6, "T2": 70e-6, "cx_err": 1e-3},
    1: {"T1": 48e-6, "T2": 68e-6, "cx_err": 1.2e-3},
    2: {"T1": 46e-6, "T2": 65e-6, "cx_err": 1.5e-3},
    3: {"T1": 45e-6, "T2": 63e-6, "cx_err": 1.8e-3},
}


def extract_backend_noise(backend_name: str = "ibmq_jakarta") -> Dict[int, Dict[str, float]]:
    """Retorna T1, T2 e erro de CX por qubit.

    Usa qiskit-ibm-provider quando disponível; caso contrário, retorna um
    dicionário sintético estável para permitir testes off-line.
    """
    try:
        from qiskit_ibm_provider import IBMProvider  # type: ignore
        provider = IBMProvider()
        backend = provider.get_backend(backend_name)
        props = backend.properties()
        n_qubits = backend.configuration().n_qubits
        noise = {}
        for q in range(n_qubits):
            try:
                cx_err = props.gate_error("cx", (q, (q + 1) % n_qubits))
            except Exception:
                cx_err = 1e-3
            noise[q] = {
                "T1": props.t1(q),
                "T2": props.t2(q),
                "cx_err": cx_err,
            }
        logger.info("Noise metrics collected from backend %s", backend_name)
        return noise
    except Exception as exc:  # pragma: no cover - fallback para ambiente sem qiskit
        logger.debug("Using synthetic noise profile: %s", exc)
        return DEFAULT_NOISE.copy()


def _rank_qubits(noise: Dict[int, Dict[str, float]]) -> List[int]:
    rank = {q: noise[q]["T1"] / 50e-6 + (1 - noise[q]["cx_err"]) for q in noise}
    ordered = sorted(rank.items(), key=lambda kv: kv[1], reverse=True)
    return [k for k, _ in ordered]


def _insert_dd(qubits: List[int], active: List[int]):
    """Sequência XYXY em qubits ociosos."""
    import pennylane as qml  # lazy import

    idle = [q for q in qubits if q not in active]
    for q in idle:
        qml.RX(float(np.pi), wires=int(q))  # type: ignore[arg-type]
        qml.RY(float(np.pi), wires=int(q))  # type: ignore[arg-type]
        qml.RX(float(np.pi), wires=int(q))  # type: ignore[arg-type]
        qml.RY(float(np.pi), wires=int(q))  # type: ignore[arg-type]


def noise_aware_ansatz(n_qubits: int, n_layers: int, backend_name: str = "ibmq_jakarta"):
    """Constrói função de circuito PennyLane adaptada ao ruído.

    - Prioriza qubits com maior T1/T2
    - Evita pares com maior erro de CX
    - Insere Dynamical Decoupling em janelas ociosas
    """
    import pennylane as qml  # lazy import

    noise = extract_backend_noise(backend_name)
    ranked = _rank_qubits(noise)[:n_qubits]

    def circuit(params, x):
        # Encoding
        active = []
        for i, q in enumerate(ranked):
            if i < len(x):
                qml.RY(np.pi * x[i], wires=q)
                active.append(q)

        # Layers
        idx = 0
        for _ in range(n_layers):
            for q in ranked:
                if idx < len(params):
                    qml.RY(params[idx], wires=q)
                    idx += 1
            # Entangling com filtragem por erro
            for i in range(len(ranked) - 1):
                q_a, q_b = ranked[i], ranked[i + 1]
                if noise.get(q_a, {}).get("cx_err", 1.0) < 0.002:
                    qml.CNOT(wires=[q_a, q_b])
            _insert_dd(ranked, active)
        return qml.expval(qml.PauliZ(ranked[0]))

    return circuit


def ansatz_sha256(noise: Dict[int, Dict[str, float]]) -> str:
    """Hash SHA-256 do perfil de ruído para rastreabilidade A1."""
    payload = "".join([
        f"{q}:{vals['T1']:.6e}:{vals['T2']:.6e}:{vals['cx_err']:.3e}" for q, vals in sorted(noise.items())
    ])
    return hashlib.sha256(payload.encode()).hexdigest()


__all__ = ["extract_backend_noise", "noise_aware_ansatz", "ansatz_sha256"]

