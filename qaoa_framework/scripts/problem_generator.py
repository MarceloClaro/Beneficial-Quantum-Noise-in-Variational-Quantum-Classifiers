"""
Tarefa 3: Gerador de Problemas de Benchmark
Gera instâncias de problemas de otimização combinatória para QAOA.

Referências:
- Farhi et al. (2014). "Quantum Approximate Optimization Algorithm." arXiv:1411.4028
- Goemans & Williamson (1995). "Improved Approximation Algorithms for Maximum Cut."
"""

from typing import Tuple, List, Dict, Optional
from dataclasses import dataclass
import numpy as np
import networkx as nx


@dataclass
class MaxCutProblem:
    """
    Problema MaxCut para QAOA.
    
    Objetivo: Particionar vértices em dois conjuntos S e T para maximizar
             o número de arestas entre S e T.
    
    Formulação QAOA:
        C(z) = Σ_{(i,j)∈E} w_{ij} · (1 - z_i·z_j) / 2
        
        onde z_i ∈ {-1, +1} é a atribuição do vértice i
    """
    n_nodes: int
    edges: List[Tuple[int, int]]
    weights: Dict[Tuple[int, int], float]
    graph: nx.Graph
    optimal_value: Optional[float] = None
    optimal_cut: Optional[List[int]] = None
    
    def __post_init__(self):
        """Calcula solução ótima (para grafos pequenos) ou aproximação."""
        if self.n_nodes <= 20:
            # Força bruta para grafos pequenos
            self._compute_optimal_exact()
        else:
            # Aproximação gulosa para grafos grandes
            self._compute_optimal_greedy()
    
    def _compute_optimal_exact(self):
        """Calcula solução ótima por força bruta."""
        best_value = -np.inf
        best_cut = None
        
        # Enumerar todas as 2^n partições
        for i in range(2 ** self.n_nodes):
            cut = [(i >> j) & 1 for j in range(self.n_nodes)]
            value = self.evaluate_cut(cut)
            
            if value > best_value:
                best_value = value
                best_cut = cut
        
        self.optimal_value = best_value
        self.optimal_cut = best_cut
    
    def _compute_optimal_greedy(self):
        """Aproximação gulosa (Goemans-Williamson)."""
        # Heurística: alocar vértices alternadamente
        cut = [i % 2 for i in range(self.n_nodes)]
        value = self.evaluate_cut(cut)
        
        # Melhorar localmente
        improved = True
        while improved:
            improved = False
            for i in range(self.n_nodes):
                # Testar flip do vértice i
                cut[i] = 1 - cut[i]
                new_value = self.evaluate_cut(cut)
                
                if new_value > value:
                    value = new_value
                    improved = True
                else:
                    # Desfazer flip
                    cut[i] = 1 - cut[i]
        
        self.optimal_value = value
        self.optimal_cut = cut
    
    def evaluate_cut(self, cut: List[int]) -> float:
        """
        Avalia qualidade de um cut.
        
        Args:
            cut: Lista binária [0,1,0,1,...] indicando partição
        
        Returns:
            Valor objetivo (número de arestas cortadas)
        """
        value = 0.0
        for (i, j), w in self.weights.items():
            # Se vértices estão em partições diferentes
            if cut[i] != cut[j]:
                value += w
        return value
    
    def to_hamiltonian(self) -> Dict[Tuple[int, int], float]:
        """
        Converte para Hamiltoniano QAOA.
        
        Returns:
            Dicionário {(i,j): coef} representando Σ coef·Z_i·Z_j
            
        Formulação:
            H = Σ_{(i,j)} w_{ij} · (1 - Z_i·Z_j) / 2
              = constante + Σ_{(i,j)} (-w_{ij}/2) · Z_i·Z_j
        """
        hamiltonian = {}
        for (i, j), w in self.weights.items():
            hamiltonian[(i, j)] = -w / 2.0
        return hamiltonian


def gerar_grafo_erdos_renyi(
    n_nodes: int,
    edge_probability: float,
    seed: Optional[int] = None,
    weighted: bool = True
) -> MaxCutProblem:
    """
    Gera grafo aleatório Erdős-Rényi.
    
    Args:
        n_nodes: Número de vértices
        edge_probability: Probabilidade de cada aresta existir
        seed: Semente aleatória
        weighted: Se True, pesos aleatórios [0.5, 1.5]
    
    Returns:
        MaxCutProblem
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Criar grafo
    graph = nx.erdos_renyi_graph(n_nodes, edge_probability, seed=seed)
    
    # Adicionar pesos
    edges = list(graph.edges())
    weights = {}
    
    for (i, j) in edges:
        if weighted:
            w = np.random.uniform(0.5, 1.5)
        else:
            w = 1.0
        weights[(i, j)] = w
        weights[(j, i)] = w  # Simétrico
        graph[i][j]['weight'] = w
    
    return MaxCutProblem(
        n_nodes=n_nodes,
        edges=edges,
        weights=weights,
        graph=graph
    )


def gerar_grafo_regular(
    n_nodes: int,
    degree: int,
    seed: Optional[int] = None,
    weighted: bool = False
) -> MaxCutProblem:
    """
    Gera grafo d-regular (cada vértice tem grau d).
    
    Args:
        n_nodes: Número de vértices
        degree: Grau de cada vértice
        seed: Semente aleatória
        weighted: Se True, pesos aleatórios
    
    Returns:
        MaxCutProblem
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Criar grafo d-regular
    graph = nx.random_regular_graph(degree, n_nodes, seed=seed)
    
    edges = list(graph.edges())
    weights = {}
    
    for (i, j) in edges:
        if weighted:
            w = np.random.uniform(0.5, 1.5)
        else:
            w = 1.0
        weights[(i, j)] = w
        weights[(j, i)] = w
        graph[i][j]['weight'] = w
    
    return MaxCutProblem(
        n_nodes=n_nodes,
        edges=edges,
        weights=weights,
        graph=graph
    )


def gerar_grafo_completo(
    n_nodes: int,
    weighted: bool = False
) -> MaxCutProblem:
    """
    Gera grafo completo (todos os pares conectados).
    
    Args:
        n_nodes: Número de vértices
        weighted: Se True, pesos aleatórios
    
    Returns:
        MaxCutProblem
    """
    graph = nx.complete_graph(n_nodes)
    
    edges = list(graph.edges())
    weights = {}
    
    for (i, j) in edges:
        if weighted:
            w = np.random.uniform(0.5, 1.5)
        else:
            w = 1.0
        weights[(i, j)] = w
        weights[(j, i)] = w
        graph[i][j]['weight'] = w
    
    return MaxCutProblem(
        n_nodes=n_nodes,
        edges=edges,
        weights=weights,
        graph=graph
    )


def gerar_problema_benchmark(
    graph_type: str,
    n_nodes: int,
    seed: Optional[int] = None,
    **kwargs
) -> MaxCutProblem:
    """
    Gera problema de benchmark conforme especificação.
    
    Args:
        graph_type: 'erdos_renyi', 'regular', 'complete'
        n_nodes: Número de nós
        seed: Semente aleatória
        **kwargs: Parâmetros adicionais (edge_probability, degree, etc.)
    
    Returns:
        MaxCutProblem
    """
    if graph_type == 'erdos_renyi':
        edge_prob = kwargs.get('edge_probability', 0.5)
        weighted = kwargs.get('weighted', True)
        return gerar_grafo_erdos_renyi(n_nodes, edge_prob, seed, weighted)
    
    elif graph_type == 'regular':
        degree = kwargs.get('degree', 3)
        weighted = kwargs.get('weighted', False)
        return gerar_grafo_regular(n_nodes, degree, seed, weighted)
    
    elif graph_type == 'complete':
        weighted = kwargs.get('weighted', False)
        return gerar_grafo_completo(n_nodes, weighted)
    
    else:
        raise ValueError(f"Tipo de grafo '{graph_type}' não reconhecido. "
                        f"Opções: erdos_renyi, regular, complete")


# Biblioteca de instâncias de benchmark conhecidas
BENCHMARK_INSTANCES = {
    'small_erdos_renyi': {
        'graph_type': 'erdos_renyi',
        'n_nodes': 10,
        'edge_probability': 0.5,
        'weighted': False
    },
    'medium_erdos_renyi': {
        'graph_type': 'erdos_renyi',
        'n_nodes': 50,
        'edge_probability': 0.3,
        'weighted': True
    },
    'large_erdos_renyi': {
        'graph_type': 'erdos_renyi',
        'n_nodes': 100,
        'edge_probability': 0.2,
        'weighted': True
    },
    'small_regular': {
        'graph_type': 'regular',
        'n_nodes': 10,
        'degree': 3,
        'weighted': False
    },
    'large_regular': {
        'graph_type': 'regular',
        'n_nodes': 100,
        'degree': 3,
        'weighted': False
    }
}


def carregar_benchmark(nome: str, seed: Optional[int] = None) -> MaxCutProblem:
    """Carrega instância de benchmark por nome."""
    if nome not in BENCHMARK_INSTANCES:
        raise ValueError(f"Benchmark '{nome}' não encontrado. "
                        f"Disponíveis: {list(BENCHMARK_INSTANCES.keys())}")
    
    config = BENCHMARK_INSTANCES[nome]
    return gerar_problema_benchmark(seed=seed, **config)


if __name__ == "__main__":
    # Teste do gerador de problemas
    print("Testando gerador de problemas MaxCut...")
    
    # Teste 1: Grafo pequeno (força bruta)
    print("\n1. Grafo pequeno (10 nós, Erdős-Rényi):")
    problem = gerar_grafo_erdos_renyi(10, 0.5, seed=42)
    print(f"   Nós: {problem.n_nodes}")
    print(f"   Arestas: {len(problem.edges)}")
    print(f"   Valor ótimo: {problem.optimal_value:.2f}")
    print(f"   Cut ótimo: {problem.optimal_cut}")
    
    # Teste 2: Grafo grande (heurística)
    print("\n2. Grafo grande (100 nós, Erdős-Rényi):")
    problem = gerar_grafo_erdos_renyi(100, 0.2, seed=42)
    print(f"   Nós: {problem.n_nodes}")
    print(f"   Arestas: {len(problem.edges)}")
    print(f"   Valor aproximado: {problem.optimal_value:.2f}")
    
    # Teste 3: Hamiltoniano
    print("\n3. Hamiltoniano QAOA (grafo pequeno):")
    problem = gerar_grafo_erdos_renyi(5, 0.5, seed=42, weighted=False)
    hamiltonian = problem.to_hamiltonian()
    print(f"   Termos: {len(hamiltonian)}")
    print(f"   Exemplo: {list(hamiltonian.items())[:3]}")
    
    # Teste 4: Benchmarks
    print("\n4. Benchmarks disponíveis:")
    for nome in BENCHMARK_INSTANCES:
        problem = carregar_benchmark(nome, seed=42)
        print(f"   {nome}: {problem.n_nodes} nós, {len(problem.edges)} arestas")
