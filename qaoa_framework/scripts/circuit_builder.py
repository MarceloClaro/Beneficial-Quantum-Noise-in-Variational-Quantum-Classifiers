"""
Tarefa 4: Construtor de Circuito QAOA
Implementa ansatz QAOA escalável para 100 qubits com suporte MPS.

Referências:
- Farhi et al. (2014). "Quantum Approximate Optimization Algorithm"
- Zhou et al. (2020). "Quantum approximate optimization algorithm: Performance, mechanism"

Ansatz QAOA:
    |ψ(γ,β)⟩ = U(B,β_p)U(C,γ_p)...U(B,β_1)U(C,γ_1)|+⟩^⊗n
    
    onde:
    - U(C,γ) = exp(-iγC): Hamiltoniano do problema
    - U(B,β) = exp(-iβB): Hamiltoniano de mixing
"""

from typing import List, Tuple, Dict, Optional
import numpy as np

try:
    from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
    from qiskit.circuit import Parameter, ParameterVector
    QISKIT_AVAILABLE = True
except ImportError:
    QISKIT_AVAILABLE = False
    print("⚠️  Qiskit não disponível")


class QAOACircuitBuilder:
    """
    Construtor de circuitos QAOA para problemas MaxCut.
    
    Suporta:
    - Escalabilidade até 100 qubits
    - Profundidade p configurável
    - Hamiltonianos customizáveis
    """
    
    def __init__(
        self,
        n_qubits: int,
        p_layers: int,
        hamiltonian: Dict[Tuple[int, int], float]
    ):
        """
        Args:
            n_qubits: Número de qubits
            p_layers: Profundidade do circuito (número de camadas)
            hamiltonian: Hamiltoniano do problema {(i,j): coef}
        """
        if not QISKIT_AVAILABLE:
            raise ImportError("Qiskit não disponível")
        
        self.n_qubits = n_qubits
        self.p_layers = p_layers
        self.hamiltonian = hamiltonian
        
        # Parâmetros variacionais
        self.n_params = 2 * p_layers
        self.params = ParameterVector('θ', self.n_params)
        
        # Separar em γ (problem) e β (mixing)
        self.gammas = [self.params[i] for i in range(p_layers)]
        self.betas = [self.params[p_layers + i] for i in range(p_layers)]
    
    def build(self, measurements: bool = True) -> QuantumCircuit:
        """
        Constrói circuito QAOA completo.
        
        Args:
            measurements: Se True, adiciona medições
        
        Returns:
            QuantumCircuit parametrizado
        """
        # Criar registros
        qr = QuantumRegister(self.n_qubits, 'q')
        if measurements:
            cr = ClassicalRegister(self.n_qubits, 'c')
            qc = QuantumCircuit(qr, cr)
        else:
            qc = QuantumCircuit(qr)
        
        # Estado inicial: |+⟩^⊗n
        qc.h(range(self.n_qubits))
        qc.barrier(label='Initial State')
        
        # p camadas QAOA
        for layer in range(self.p_layers):
            # 1. Hamiltoniano do Problema: U(C, γ)
            self._apply_problem_hamiltonian(qc, self.gammas[layer])
            qc.barrier(label=f'Problem {layer+1}')
            
            # 2. Hamiltoniano de Mixing: U(B, β)
            self._apply_mixing_hamiltonian(qc, self.betas[layer])
            qc.barrier(label=f'Mixing {layer+1}')
        
        # Medições
        if measurements:
            qc.measure(qr, cr)
        
        return qc
    
    def _apply_problem_hamiltonian(
        self,
        qc: QuantumCircuit,
        gamma: Parameter
    ):
        """
        Aplica U(C, γ) = exp(-iγC) onde C = Σ_{(i,j)} c_{ij} Z_i Z_j
        
        Para cada aresta (i,j) com peso c_{ij}:
            CNOT(i,j) - RZ(2·γ·c_{ij}, j) - CNOT(i,j)
        
        Isso implementa exp(-iγ·c_{ij}·Z_i·Z_j)
        """
        for (i, j), coef in self.hamiltonian.items():
            # Aplicar ZZ interaction
            qc.cx(i, j)
            qc.rz(2 * gamma * coef, j)
            qc.cx(i, j)
    
    def _apply_mixing_hamiltonian(
        self,
        qc: QuantumCircuit,
        beta: Parameter
    ):
        """
        Aplica U(B, β) = exp(-iβB) onde B = Σ_i X_i
        
        Para cada qubit i:
            RX(2·β, i)
        
        Isso implementa exp(-iβ·X_i) = RX(2β)
        """
        for i in range(self.n_qubits):
            qc.rx(2 * beta, i)
    
    def get_parameter_bounds(self) -> List[Tuple[float, float]]:
        """
        Retorna bounds para otimização dos parâmetros.
        
        Returns:
            Lista de (lower, upper) para cada parâmetro
            
        Heurística:
            γ ∈ [0, π]  (problem Hamiltonian)
            β ∈ [0, π/2] (mixing Hamiltonian)
        """
        bounds = []
        
        # Bounds para γ (gammas)
        for _ in range(self.p_layers):
            bounds.append((0.0, np.pi))
        
        # Bounds para β (betas)
        for _ in range(self.p_layers):
            bounds.append((0.0, np.pi / 2))
        
        return bounds
    
    def initialize_parameters(
        self,
        strategy: str = 'random',
        seed: Optional[int] = None
    ) -> np.ndarray:
        """
        Inicializa parâmetros do QAOA.
        
        Args:
            strategy: 'random', 'zeros', 'heuristic'
            seed: Semente aleatória
        
        Returns:
            Array de parâmetros iniciais [γ_1,...,γ_p, β_1,...,β_p]
        """
        if seed is not None:
            np.random.seed(seed)
        
        if strategy == 'random':
            # Valores aleatórios nos bounds
            gammas = np.random.uniform(0, np.pi, self.p_layers)
            betas = np.random.uniform(0, np.pi/2, self.p_layers)
            return np.concatenate([gammas, betas])
        
        elif strategy == 'zeros':
            # Começar de zero
            return np.zeros(self.n_params)
        
        elif strategy == 'heuristic':
            # Heurística baseada em literatura
            # γ ≈ 0.5, β ≈ π/4 para p=1
            # Para p>1, interpolar
            gammas = np.linspace(0.3, 0.7, self.p_layers)
            betas = np.linspace(np.pi/4, np.pi/6, self.p_layers)
            return np.concatenate([gammas, betas])
        
        else:
            raise ValueError(f"Estratégia '{strategy}' não reconhecida. "
                           f"Opções: random, zeros, heuristic")
    
    def get_circuit_depth(self) -> int:
        """Retorna profundidade do circuito (número de camadas)."""
        qc = self.build(measurements=False)
        return qc.depth()
    
    def get_circuit_stats(self) -> Dict[str, int]:
        """
        Retorna estatísticas do circuito.
        
        Returns:
            Dict com qubits, depth, gates, cx_gates
        """
        qc = self.build(measurements=False)
        
        # Contar gates
        gate_counts = qc.count_ops()
        cx_count = gate_counts.get('cx', 0)
        total_gates = sum(gate_counts.values())
        
        return {
            'qubits': self.n_qubits,
            'depth': qc.depth(),
            'total_gates': total_gates,
            'cx_gates': cx_count,
            'parameters': self.n_params
        }


def criar_circuito_qaoa_from_config(config: Dict) -> QAOACircuitBuilder:
    """
    Cria circuito QAOA a partir de configuração YAML.
    
    Args:
        config: Dicionário com configuração (de experiment_qaoa.yaml)
    
    Returns:
        QAOACircuitBuilder configurado
        
    Exemplo:
        >>> config = yaml.safe_load(open('experiment_qaoa.yaml'))
        >>> builder = criar_circuito_qaoa_from_config(config)
    """
    from .problem_generator import gerar_problema_benchmark
    
    # Gerar problema
    problem_config = config['problem']
    problem = gerar_problema_benchmark(
        graph_type=problem_config['graph_type'],
        n_nodes=problem_config['n_nodes'],
        edge_probability=problem_config.get('edge_probability', 0.5)
    )
    
    # Obter Hamiltoniano
    hamiltonian = problem.to_hamiltonian()
    
    # Criar builder
    model_config = config['model']
    builder = QAOACircuitBuilder(
        n_qubits=problem_config['n_nodes'],
        p_layers=model_config['p_layers'],
        hamiltonian=hamiltonian
    )
    
    return builder


if __name__ == "__main__":
    print("Testando construtor de circuito QAOA...")
    
    if QISKIT_AVAILABLE:
        # Teste 1: Circuito pequeno
        print("\n1. Circuito QAOA pequeno (5 qubits, p=2):")
        hamiltonian = {
            (0, 1): -0.5,
            (1, 2): -0.5,
            (2, 3): -0.5,
            (3, 4): -0.5,
            (0, 4): -0.5
        }
        
        builder = QAOACircuitBuilder(n_qubits=5, p_layers=2, hamiltonian=hamiltonian)
        qc = builder.build()
        
        stats = builder.get_circuit_stats()
        print(f"   Qubits: {stats['qubits']}")
        print(f"   Depth: {stats['depth']}")
        print(f"   Gates: {stats['total_gates']} (CX: {stats['cx_gates']})")
        print(f"   Parâmetros: {stats['parameters']}")
        
        # Teste 2: Inicialização de parâmetros
        print("\n2. Estratégias de inicialização:")
        for strategy in ['random', 'zeros', 'heuristic']:
            params = builder.initialize_parameters(strategy, seed=42)
            print(f"   {strategy}: {params}")
        
        # Teste 3: Circuito grande (100 qubits)
        print("\n3. Circuito QAOA grande (100 qubits, p=3):")
        # Hamiltoniano de grafo esparso (10% densidade)
        hamiltonian_large = {}
        np.random.seed(42)
        for i in range(100):
            for j in range(i+1, 100):
                if np.random.random() < 0.1:
                    hamiltonian_large[(i, j)] = -np.random.uniform(0.5, 1.5) / 2
        
        builder_large = QAOACircuitBuilder(n_qubits=100, p_layers=3, hamiltonian=hamiltonian_large)
        stats_large = builder_large.get_circuit_stats()
        
        print(f"   Qubits: {stats_large['qubits']}")
        print(f"   Depth: {stats_large['depth']}")
        print(f"   Gates: {stats_large['total_gates']} (CX: {stats_large['cx_gates']})")
        print(f"   Parâmetros: {stats_large['parameters']}")
        print(f"   Arestas do grafo: {len(hamiltonian_large)}")
    
    else:
        print("❌ Qiskit não disponível para testes")
