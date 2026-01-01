#!/usr/bin/env python3
"""
Script para executar QAOA 100 qubits completo com ruído benéfico
Sem necessidade de entrada interativa
"""

import os
import sys
import time
import logging
from pathlib import Path

# Importar framework QAOA
try:
    from framework_qaoa_100qubits import (
        ConfigQAOA,
        ConstrutorCircuitoQAOA,
        OtimizadorQAOA,
        AnalisadorHiperparametrosQAOA,
        VisualizadorQAOA,
        experimento_completo_ruido_benefico
    )
    print("✅ Framework QAOA importado com sucesso")
except ImportError as e:
    print(f"❌ Erro ao importar framework: {e}")
    sys.exit(1)

# Configurar logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-8s | %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    """Executar experimento completo QAOA 100 qubits."""
    
    print("\n" + "="*80)
    print("🚀 EXPERIMENTO COMPLETO: QAOA 100 QUBITS COM RUÍDO BENÉFICO")
    print("="*80)
    print("\nConfiguração do Experimento:")
    print("-" * 80)
    print("  • Qubits: 100")
    print("  • P-layers: 3")
    print("  • Problema: MaxCut (grafos aleatórios)")
    print("  • Tipos de ruído: Depolarizing, Phase Damping, Amplitude Damping")
    print("  • Níveis de ruído: 0.0%, 0.05%, 0.1%, 0.15%, 0.2%")
    print("  • Otimizador: COBYLA")
    print("  • Iterações máximas: 100")
    print("  • Shots: 1024 medições por iteração")
    print("  • Seeds: 5 repetições (42-46)")
    print("  • Tempo estimado: 30-45 minutos")
    print("-" * 80)
    
    try:
        print("\n📊 Iniciando experimento completo...\n")
        start_time = time.time()
        
        # Executar experimento completo
        resultado = experimento_completo_ruido_benefico(
            n_qubits=100,
            p_layers=3,
            problema='maxcut',
            pasta_saida='resultados_qaoa_100qubits_completo'
        )
        
        elapsed = time.time() - start_time
        
        print("\n" + "="*80)
        print("✅ EXPERIMENTO CONCLUÍDO COM SUCESSO!")
        print("="*80)
        print(f"\n⏱️  Tempo total de execução: {elapsed/60:.1f} minutos ({elapsed:.1f} segundos)")
        
        if resultado:
            print("\n📈 Resultados disponíveis em:")
            print(f"   - Pasta: {Path('resultados_qaoa_100qubits_completo').absolute()}")
            print("   - Arquivos CSV com dados detalhados")
            print("   - Gráficos em HTML (Plotly)")
            print("   - Análises estatísticas")
        
        print("\n" + "="*80 + "\n")
        return 0
        
    except KeyboardInterrupt:
        print("\n⚠️  Experimento cancelado pelo usuário")
        return 1
    except Exception as e:
        print(f"\n❌ Erro durante experimento: {e}")
        logger.exception("Erro detalhado:")
        return 1


if __name__ == "__main__":
    sys.exit(main())
