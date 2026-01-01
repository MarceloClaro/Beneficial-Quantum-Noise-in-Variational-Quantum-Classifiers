#!/usr/bin/env python3
"""
Script de execução QAOA 100 qubits com análise de ruído benéfico
Executa o experimento completo conforme solicitado
"""

import sys
import os

# Adicionar diretório atual ao path
sys.path.insert(0, os.getcwd())

# Importar o script principal
import executar_qaoa_100qubits

# Executar
if __name__ == "__main__":
    print("\n" + "="*80)
    print("🚀 EXPERIMENTO COMPLETO: QAOA 100 QUBITS COM ANÁLISE DE RUÍDO BENÉFICO")
    print("="*80 + "\n")
    
    try:
        # O script será executado quando importado
        print("✅ Experimento iniciado com sucesso!")
    except Exception as e:
        print(f"❌ Erro ao executar: {e}")
        sys.exit(1)
