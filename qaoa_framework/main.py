"""
QAOA Framework: Main Execution Script
Integra todos os módulos para execução completa do experimento QAOA.

Autor: Projeto Beneficial Quantum Noise
Data: 2025-12-26
"""

import os
import sys
import yaml
import json
import time
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional

import numpy as np
import pandas as pd

# Adicionar diretório de scripts ao path
sys.path.insert(0, str(Path(__file__).parent / 'scripts'))

# Importar módulos do framework
from scripts.noise_models import criar_noise_model_qiskit, aplicar_schedule_ruido
from scripts.problem_generator import gerar_problema_benchmark
from scripts.circuit_builder import QAOACircuitBuilder

try:
    from qiskit import transpile
    from qiskit_aer import AerSimulator
    QISKIT_AVAILABLE = True
except ImportError:
    QISKIT_AVAILABLE = False
    print("❌ Qiskit não disponível. Instale com: pip install qiskit qiskit-aer")
    sys.exit(1)

# Configurar logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-8s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class QAOAExperiment:
    """Gerenciador de experimentos QAOA com análise de ruído benéfico."""
    
    def __init__(self, config_path: str):
        """
        Args:
            config_path: Caminho para arquivo experiment_qaoa.yaml
        """
        self.config_path = Path(config_path)
        self.config = self._load_config()
        self.run_id = self._generate_run_id()
        self.results = []
        
        # Criar diretório de resultados
        self.output_dir = Path(self.config['output']['output_dir']) / self.run_id
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Experimento QAOA inicializado: {self.run_id}")
        logger.info(f"Resultados serão salvos em: {self.output_dir}")
    
    def _load_config(self) -> Dict:
        """Carrega configuração do YAML."""
        if not self.config_path.exists():
            raise FileNotFoundError(f"Arquivo de configuração não encontrado: {self.config_path}")
        
        with open(self.config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        logger.info(f"Configuração carregada: {self.config_path}")
        return config
    
    def _generate_run_id(self) -> str:
        """Gera ID único para a execução."""
        run_id = self.config['run'].get('run_id', 'AUTO')
        
        if run_id == 'AUTO':
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            run_id = f"qaoa_run_{timestamp}"
        
        return run_id
    
    def executar(self):
        """Executa experimento completo."""
        logger.info("="*80)
        logger.info("INICIANDO EXPERIMENTO QAOA")
        logger.info("="*80)
        
        inicio_total = time.time()
        
        # 1. Gerar problema
        logger.info("\n[1/5] Gerando problema de benchmark...")
        problem = self._gerar_problema()
        
        # 2. Construir circuito
        logger.info("\n[2/5] Construindo circuito QAOA...")
        builder = self._construir_circuito(problem)
        
        # 3. Executar experimentos com diferentes níveis de ruído
        logger.info("\n[3/5] Executando experimentos com ruído...")
        self._executar_grid_search(builder, problem)
        
        # 4. Salvar resultados
        logger.info("\n[4/5] Salvando resultados...")
        self._salvar_resultados()
        
        # 5. Gerar manifesto
        logger.info("\n[5/5] Gerando manifesto de execução...")
        self._gerar_manifesto(problem, builder, time.time() - inicio_total)
        
        logger.info("\n" + "="*80)
        logger.info(f"EXPERIMENTO CONCLUÍDO EM {time.time() - inicio_total:.1f}s")
        logger.info(f"Resultados salvos em: {self.output_dir}")
        logger.info("="*80 + "\n")
    
    def _gerar_problema(self):
        """Gera problema MaxCut conforme configuração."""
        config_prob = self.config['problem']
        
        problem = gerar_problema_benchmark(
            graph_type=config_prob['graph_type'],
            n_nodes=config_prob['n_nodes'],
            edge_probability=config_prob.get('edge_probability', 0.5),
            seed=self.config['reproducibility']['seeds'][0]
        )
        
        logger.info(f"  Tipo: {config_prob['type']}")
        logger.info(f"  Nós: {problem.n_nodes}")
        logger.info(f"  Arestas: {len(problem.edges)}")
        logger.info(f"  Valor ótimo estimado: {problem.optimal_value:.2f}")
        
        return problem
    
    def _construir_circuito(self, problem):
        """Constrói circuito QAOA."""
        hamiltonian = problem.to_hamiltonian()
        
        builder = QAOACircuitBuilder(
            n_qubits=problem.n_nodes,
            p_layers=self.config['model']['p_layers'],
            hamiltonian=hamiltonian
        )
        
        stats = builder.get_circuit_stats()
        logger.info(f"  Qubits: {stats['qubits']}")
        logger.info(f"  P-layers: {self.config['model']['p_layers']}")
        logger.info(f"  Depth: {stats['depth']}")
        logger.info(f"  Gates: {stats['total_gates']} (CX: {stats['cx_gates']})")
        
        return builder
    
    def _executar_grid_search(self, builder, problem):
        """Executa grid search sobre níveis de ruído."""
        seeds = self.config['reproducibility']['seeds']
        niveis_ruido = self.config['noise']['params']['p']
        tipo_ruido = self.config['noise']['model']
        
        total_exp = len(seeds) * len(niveis_ruido)
        logger.info(f"  Total de experimentos: {total_exp}")
        
        exp_count = 0
        
        for seed in seeds:
            for nivel_ruido in niveis_ruido:
                exp_count += 1
                logger.info(f"\n  [{exp_count}/{total_exp}] Seed={seed}, Ruído={nivel_ruido:.4f}")
                
                resultado = self._executar_single_run(
                    builder, problem, tipo_ruido, nivel_ruido, seed
                )
                
                self.results.append(resultado)
    
    def _executar_single_run(self, builder, problem, tipo_ruido, nivel_ruido, seed):
        """Executa uma única execução QAOA."""
        np.random.seed(seed)
        
        # Criar modelo de ruído
        if self.config['noise']['enabled'] and nivel_ruido > 0:
            noise_model = criar_noise_model_qiskit(tipo_ruido, nivel_ruido)
        else:
            noise_model = None
        
        # Criar simulador
        backend_config = self.config['frameworks']['qiskit']
        if noise_model:
            simulator = AerSimulator(
                noise_model=noise_model,
                method=backend_config['method']
            )
        else:
            simulator = AerSimulator(method=backend_config['method'])
        
        # Construir circuito
        qc = builder.build()
        
        # Inicializar parâmetros
        params_init = builder.initialize_parameters(
            strategy=self.config['optimization']['initial_params'],
            seed=seed
        )
        
        # Otimizar (simplificado - versão completa em cost_function.py)
        from scipy.optimize import minimize
        
        def objective(params):
            # Bind parameters
            bound_circuit = qc.assign_parameters({builder.params[i]: params[i] for i in range(len(params))})
            
            # Transpile e executar
            transpiled = transpile(bound_circuit, simulator)
            job = simulator.run(transpiled, shots=backend_config['shots'])
            result = job.result()
            counts = result.get_counts()
            
            # Calcular energia
            energia = self._calcular_energia(counts, problem)
            return energia
        
        # Otimizar
        inicio = time.time()
        opt_result = minimize(
            objective,
            params_init,
            method=self.config['optimization']['optimizer'],
            options={'maxiter': self.config['optimization']['maxiter']}
        )
        tempo_exec = time.time() - inicio
        
        # Calcular approximation ratio
        energia_final = opt_result.fun
        approx_ratio = abs(energia_final / problem.optimal_value) if problem.optimal_value != 0 else 0.0
        
        logger.info(f"      Energia: {energia_final:.4f}")
        logger.info(f"      Approx Ratio: {approx_ratio:.4f}")
        logger.info(f"      Tempo: {tempo_exec:.2f}s")
        
        return {
            'seed': seed,
            'tipo_ruido': tipo_ruido if nivel_ruido > 0 else 'sem_ruido',
            'nivel_ruido': nivel_ruido,
            'energia_final': energia_final,
            'approx_ratio': approx_ratio,
            'iteracoes': getattr(opt_result, 'nit', getattr(opt_result, 'nfev', 0)),
            'tempo_execucao': tempo_exec,
            'convergiu': opt_result.success
        }
    
    def _calcular_energia(self, counts, problem):
        """Calcula energia esperada do Hamiltoniano."""
        energia = 0.0
        total_shots = sum(counts.values())
        
        for bitstring, count in counts.items():
            # Converter bitstring para cut
            cut = [int(b) for b in bitstring[::-1]]  # Reverter (Qiskit little-endian)
            
            # Calcular contribuição
            energia_config = 0.0
            for (i, j), coef in problem.to_hamiltonian().items():
                # Z_i Z_j = (-1)^(s_i ⊕ s_j)
                zi_zj = 1 if cut[i] == cut[j] else -1
                energia_config += coef * zi_zj
            
            prob = count / total_shots
            energia += prob * energia_config
        
        return energia
    
    def _salvar_resultados(self):
        """Salva resultados em CSV."""
        df = pd.DataFrame(self.results)
        
        arquivo = self.output_dir / 'resultados.csv'
        df.to_csv(arquivo, index=False)
        
        logger.info(f"  Resultados salvos: {arquivo}")
        
        # Análise sumária
        resumo = df.groupby(['tipo_ruido', 'nivel_ruido']).agg({
            'energia_final': ['mean', 'std'],
            'approx_ratio': ['mean', 'std']
        })
        
        arquivo_resumo = self.output_dir / 'resumo.csv'
        resumo.to_csv(arquivo_resumo)
        
        logger.info(f"  Resumo salvo: {arquivo_resumo}")
    
    def _gerar_manifesto(self, problem, builder, tempo_total):
        """Gera manifesto de execução (reprodutibilidade)."""
        manifesto = {
            'run_id': self.run_id,
            'timestamp': datetime.now().isoformat(),
            'config_file': str(self.config_path),
            'tempo_total_s': tempo_total,
            
            'problem': {
                'type': self.config['problem']['type'],
                'n_nodes': problem.n_nodes,
                'n_edges': len(problem.edges),
                'optimal_value': problem.optimal_value
            },
            
            'model': {
                'algorithm': 'QAOA',
                'p_layers': builder.p_layers,
                'n_parameters': builder.n_params
            },
            
            'circuit': builder.get_circuit_stats(),
            
            'experiments': {
                'total': len(self.results),
                'seeds': self.config['reproducibility']['seeds'],
                'noise_levels': self.config['noise']['params']['p']
            },
            
            'environment': {
                'python_version': sys.version,
                'qiskit_version': self._get_qiskit_version()
            }
        }
        
        arquivo = self.output_dir / 'manifesto.json'
        with open(arquivo, 'w') as f:
            json.dump(manifesto, f, indent=2)
        
        logger.info(f"  Manifesto salvo: {arquivo}")
    
    def _get_qiskit_version(self):
        """Obtém versão do Qiskit."""
        try:
            import qiskit
            return qiskit.__version__
        except:
            return "unknown"


def main():
    """Função principal."""
    import argparse
    
    parser = argparse.ArgumentParser(description='QAOA Framework - Execução de Experimentos')
    parser.add_argument(
        '--config',
        type=str,
        default='qaoa_framework/configs/experiment_qaoa.yaml',
        help='Caminho para arquivo de configuração YAML'
    )
    
    args = parser.parse_args()
    
    # Executar experimento
    experiment = QAOAExperiment(args.config)
    experiment.executar()


if __name__ == "__main__":
    main()
