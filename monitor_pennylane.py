"""
monitor_pennylane.py
Monitor em tempo real da execução do framework PennyLane
Rastreia progresso, tempo estimado e artefatos gerados
"""

import os
import json
import time
import logging
from pathlib import Path
from datetime import datetime, timedelta

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-8s | %(message)s'
)
logger = logging.getLogger("monitor_pennylane")


def monitorar_progresso():
    """Monitorar execução do framework PennyLane em tempo real"""
    
    logger.info("=" * 80)
    logger.info("MONITOR: Framework PennyLane Investigativo")
    logger.info("=" * 80)
    
    resultado_dir = None
    inicio = time.time()
    
    while True:
        try:
            # Procurar por diretório de resultados
            dirs = list(Path(".").glob("resultados_*"))
            
            if dirs and not resultado_dir:
                resultado_dir = sorted(dirs)[-1]
                logger.info(f"\n✓ Diretório de resultados encontrado: {resultado_dir}")
            
            if resultado_dir:
                resultado_dir = Path(resultado_dir)
                
                # Contar arquivos de circuito
                circuitos = list(resultado_dir.glob("circuito_*.png"))
                barren_plateaus = list(resultado_dir.glob("barren3d_*.png"))
                logs = list(resultado_dir.glob("*.log"))
                jsons = list(resultado_dir.glob("*.json"))
                
                # Verificar log principal
                log_files = list(resultado_dir.glob("execution_log_*.log"))
                if log_files:
                    log_file = log_files[0]
                    with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
                        linhas = f.readlines()
                    
                    # Extrair último trial
                    for linha in reversed(linhas[-50:]):
                        if "[" in linha and "/" in linha:
                            logger.info(f"Progresso: {linha.strip()}")
                            break
                
                tempo_decorrido = time.time() - inicio
                horas = tempo_decorrido / 3600
                
                logger.info(f"\n📊 Status (tempo decorrido: {horas:.1f}h):")
                logger.info(f"  • Circuitos: {len(circuitos)}")
                logger.info(f"  • Barren Plateaus: {len(barren_plateaus)}")
                logger.info(f"  • Logs: {len(logs)}")
                logger.info(f"  • JSONs: {len(jsons)}")
                
                # Estimativa de conclusão
                if len(circuitos) > 0:
                    taxa = len(circuitos) / tempo_decorrido
                    tempo_total_estimado = 600 / taxa
                    tempo_restante = tempo_total_estimado - tempo_decorrido

                    logger.info("\nEstimativa:")
                    logger.info(f"  • Taxa: {taxa:.4f} trials/sec")
                    logger.info(f"  • Tempo total: {tempo_total_estimado/3600:.1f}h")
                    logger.info(f"  • Restante: {tempo_restante/3600:.1f}h")
                    conclusao = datetime.now() + timedelta(seconds=tempo_restante)
                    logger.info(f"  • Conclusão: {conclusao.strftime('%H:%M:%S')}")
            
            else:
                logger.info("⏳ Aguardando início da execução...")
            
            logger.info("-" * 80)
            time.sleep(30)  # Verificar a cada 30 segundos
            
        except KeyboardInterrupt:
            logger.info("\n✓ Monitoramento interrompido")
            break
        except Exception as e:
            logger.warning(f"Erro ao monitorar: {e}")
            time.sleep(60)


if __name__ == "__main__":
    monitorar_progresso()
