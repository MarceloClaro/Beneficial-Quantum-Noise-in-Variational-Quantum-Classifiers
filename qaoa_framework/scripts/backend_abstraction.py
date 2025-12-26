"""
Backend Abstraction Layer para QAOA Framework

Este módulo fornece uma camada de abstração unificada para execução de circuitos
quânticos em diferentes backends Qiskit (simuladores e hardware real).

Autor: QAOA Framework Team
Tarefa: Task 7 - Backend Abstraction Layer
"""

from typing import Dict, Any, Optional, List, Union
from dataclasses import dataclass
import warnings

try:
    from qiskit import QuantumCircuit, transpile
    from qiskit_aer import AerSimulator
    from qiskit.providers import Backend, JobV1 as Job
    from qiskit.result import Result
    from qiskit_aer.noise import NoiseModel
except ImportError:
    raise ImportError(
        "Qiskit não está instalado. Execute: pip install qiskit qiskit-aer"
    )


@dataclass
class BackendInfo:
    """Informações sobre o backend."""
    name: str
    backend_type: str  # 'simulator' ou 'hardware'
    max_qubits: int
    max_shots: int
    coupling_map: Optional[List] = None
    basis_gates: Optional[List] = None
    supports_noise: bool = True
    supports_mps: bool = False


class QuantumBackend:
    """
    Camada de abstração para backends quânticos.
    
    Suporta:
    - Simuladores Aer (statevector, qasm, matrix_product_state)
    - Hardware IBMQ (requer conta configurada)
    
    Exemplo:
        >>> backend = QuantumBackend('aer_simulator', method='matrix_product_state')
        >>> result = backend.execute(circuit, shots=4096)
        >>> counts = result.get_counts()
    """
    
    def __init__(
        self,
        backend_name: str = 'aer_simulator',
        method: Optional[str] = None,
        shots: int = 4096,
        seed_simulator: Optional[int] = None,
        optimization_level: int = 1,
        noise_model: Optional[NoiseModel] = None,
        provider: Optional[str] = None,
        hub: Optional[str] = None,
        group: Optional[str] = None,
        project: Optional[str] = None,
        **kwargs
    ):
        """
        Inicializa o backend.
        
        Args:
            backend_name: Nome do backend ('aer_simulator', 'ibmq_manila', etc.)
            method: Método de simulação ('automatic', 'statevector', 'matrix_product_state')
            shots: Número de shots para execução
            seed_simulator: Seed para reprodutibilidade
            optimization_level: Nível de otimização da transpilação (0-3)
            noise_model: Modelo de ruído (para simuladores)
            provider: Provider IBMQ (ex: 'ibm-q')
            hub: Hub IBMQ
            group: Grupo IBMQ
            project: Projeto IBMQ
        """
        self.backend_name = backend_name
        self.method = method
        self.shots = shots
        self.seed_simulator = seed_simulator
        self.optimization_level = optimization_level
        self.noise_model = noise_model
        self.provider_name = provider
        self.hub = hub
        self.group = group
        self.project = project
        self.kwargs = kwargs
        
        # Cache do backend
        self._backend: Optional[Backend] = None
        self._backend_info: Optional[BackendInfo] = None
        
        # Inicializa o backend
        self._initialize_backend()
    
    def _initialize_backend(self):
        """Inicializa o backend apropriado."""
        if self.backend_name.startswith('ibmq_') or self.provider_name:
            self._initialize_ibmq_backend()
        else:
            self._initialize_aer_backend()
    
    def _initialize_aer_backend(self):
        """Inicializa um backend Aer."""
        try:
            # Cria simulador Aer
            self._backend = AerSimulator()
            
            # Configura método se especificado
            if self.method:
                backend_options = {'method': self.method}
                if self.seed_simulator is not None:
                    backend_options['seed_simulator'] = self.seed_simulator
                self._backend.set_options(**backend_options)
            
            # Determina capacidades
            supports_mps = self.method == 'matrix_product_state'
            max_qubits = 100 if supports_mps else 30
            
            self._backend_info = BackendInfo(
                name=self.backend_name,
                backend_type='simulator',
                max_qubits=max_qubits,
                max_shots=1_000_000,
                supports_noise=True,
                supports_mps=supports_mps
            )
            
        except Exception as e:
            raise RuntimeError(f"Erro ao inicializar backend Aer: {e}")
    
    def _initialize_ibmq_backend(self):
        """Inicializa um backend IBMQ."""
        try:
            from qiskit_ibm_runtime import QiskitRuntimeService
            
            # Carrega provider
            if self.hub and self.group and self.project:
                service = QiskitRuntimeService(
                    channel="ibm_quantum",
                    instance=f"{self.hub}/{self.group}/{self.project}"
                )
            else:
                service = QiskitRuntimeService()
            
            # Obtém backend
            self._backend = service.backend(self.backend_name)
            
            # Obtém informações do backend
            config = self._backend.configuration()
            self._backend_info = BackendInfo(
                name=config.backend_name,
                backend_type='hardware',
                max_qubits=config.n_qubits,
                max_shots=config.max_shots,
                coupling_map=config.coupling_map if hasattr(config, 'coupling_map') else None,
                basis_gates=config.basis_gates if hasattr(config, 'basis_gates') else None,
                supports_noise=False,  # Hardware já tem ruído real
                supports_mps=False
            )
            
        except ImportError:
            raise ImportError(
                "qiskit-ibm-runtime não está instalado. "
                "Execute: pip install qiskit-ibm-runtime"
            )
        except Exception as e:
            raise RuntimeError(f"Erro ao inicializar backend IBMQ: {e}")
    
    @classmethod
    def from_config(cls, config: Dict[str, Any]) -> 'QuantumBackend':
        """
        Cria backend a partir de configuração YAML.
        
        Args:
            config: Dicionário de configuração do YAML
            
        Returns:
            Instância de QuantumBackend
            
        Exemplo:
            >>> config = yaml.safe_load(open('experiment_qaoa.yaml'))
            >>> backend = QuantumBackend.from_config(config['frameworks']['qiskit'])
        """
        return cls(**config)
    
    def execute(
        self,
        circuit: QuantumCircuit,
        noise_model: Optional[NoiseModel] = None,
        shots: Optional[int] = None,
        **kwargs
    ) -> Result:
        """
        Executa um circuito quântico.
        
        Args:
            circuit: Circuito a ser executado
            noise_model: Modelo de ruído (sobrescreve o padrão)
            shots: Número de shots (sobrescreve o padrão)
            **kwargs: Argumentos adicionais para transpile/execute
            
        Returns:
            Resultado da execução
            
        Exemplo:
            >>> circuit = QuantumCircuit(5)
            >>> result = backend.execute(circuit, shots=8192)
        """
        if self._backend is None:
            raise RuntimeError("Backend não foi inicializado")
        
        # Usa shots padrão se não especificado
        if shots is None:
            shots = self.shots
        
        # Usa noise model padrão se não especificado e for simulador
        if noise_model is None and self._backend_info.backend_type == 'simulator':
            noise_model = self.noise_model
        
        # Transpila o circuito
        transpiled_circuit = transpile(
            circuit,
            backend=self._backend,
            optimization_level=self.optimization_level,
            **kwargs
        )
        
        # Prepara opções de execução
        run_options = {'shots': shots}
        
        # Adiciona seed se especificado
        if self.seed_simulator is not None and self._backend_info.backend_type == 'simulator':
            run_options['seed_simulator'] = self.seed_simulator
        
        # Adiciona noise model se especificado
        if noise_model is not None and self._backend_info.backend_type == 'simulator':
            run_options['noise_model'] = noise_model
        
        # Executa
        try:
            job = self._backend.run(transpiled_circuit, **run_options)
            result = job.result()
            return result
        except Exception as e:
            raise RuntimeError(f"Erro ao executar circuito: {e}")
    
    def get_backend_info(self) -> Dict[str, Any]:
        """
        Retorna informações sobre o backend.
        
        Returns:
            Dicionário com informações do backend
        """
        if self._backend_info is None:
            return {}
        
        return {
            'name': self._backend_info.name,
            'type': self._backend_info.backend_type,
            'max_qubits': self._backend_info.max_qubits,
            'max_shots': self._backend_info.max_shots,
            'coupling_map': self._backend_info.coupling_map,
            'basis_gates': self._backend_info.basis_gates,
            'supports_noise': self._backend_info.supports_noise,
            'supports_mps': self._backend_info.supports_mps
        }
    
    def get_noise_model(self) -> Optional[NoiseModel]:
        """
        Extrai o modelo de ruído do backend (apenas para hardware).
        
        Returns:
            Modelo de ruído do hardware ou None
        """
        if self._backend_info.backend_type != 'hardware':
            return None
        
        try:
            from qiskit_aer.noise import NoiseModel
            noise_model = NoiseModel.from_backend(self._backend)
            return noise_model
        except Exception as e:
            warnings.warn(f"Não foi possível extrair noise model do hardware: {e}")
            return None
    
    def __repr__(self) -> str:
        """Representação em string."""
        info = self.get_backend_info()
        return (
            f"QuantumBackend(name='{info['name']}', "
            f"type='{info['type']}', "
            f"max_qubits={info['max_qubits']}, "
            f"shots={self.shots})"
        )


def criar_backend_padrao(config: Optional[Dict[str, Any]] = None) -> QuantumBackend:
    """
    Cria um backend com configuração padrão.
    
    Args:
        config: Configuração opcional (usa padrões se None)
        
    Returns:
        Backend configurado
        
    Exemplo:
        >>> backend = criar_backend_padrao()
        >>> # Ou com config personalizada
        >>> backend = criar_backend_padrao({'backend': 'aer_simulator', 'shots': 8192})
    """
    default_config = {
        'backend_name': 'aer_simulator',
        'method': 'matrix_product_state',
        'shots': 4096,
        'optimization_level': 1
    }
    
    if config:
        default_config.update(config)
    
    return QuantumBackend(**default_config)


# Exemplo de uso
if __name__ == "__main__":
    print("=" * 80)
    print("BACKEND ABSTRACTION LAYER - QAOA Framework")
    print("=" * 80)
    print()
    
    # 1. Backend simulador com MPS
    print("1. Criando backend simulador com MPS (100 qubits)...")
    backend_mps = QuantumBackend(
        backend_name='aer_simulator',
        method='matrix_product_state',
        shots=4096
    )
    print(f"   {backend_mps}")
    print(f"   Info: {backend_mps.get_backend_info()}")
    print()
    
    # 2. Backend padrão
    print("2. Criando backend padrão...")
    backend_default = criar_backend_padrao()
    print(f"   {backend_default}")
    print()
    
    # 3. Teste de execução (circuito simples)
    print("3. Testando execução de circuito...")
    from qiskit import QuantumCircuit
    
    circuit = QuantumCircuit(3)
    circuit.h(0)
    circuit.cx(0, 1)
    circuit.cx(1, 2)
    circuit.measure_all()
    
    result = backend_mps.execute(circuit, shots=1000)
    counts = result.get_counts()
    print(f"   Counts: {counts}")
    print()
    
    print("✅ Backend abstraction funcionando corretamente!")
