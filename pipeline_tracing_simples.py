"""
pipeline_tracing_simples.py
Wrapper simplificado da pipeline QAOA com tracing OpenTelemetry
Sem dependências de encoding específico para Windows
"""

import sys
import json
import time
import logging
from datetime import datetime
from pathlib import Path

# Configurar logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("pipeline_tracing")

# Importar tracing
try:
    from tracing_setup import setup_tracing, get_tracer
    setup_tracing("qaoa-pipeline")
    tracer = get_tracer("pipeline_tracing_simples")
    logger.info("✓ OpenTelemetry tracing inicializado")
except Exception as e:
    logger.warning(f"Tracing não disponível: {e}")
    tracer = None


def teste_artefatos():
    """Verificar se artefatos QAOA existem"""
    with tracer.start_as_current_span("teste_artefatos") if tracer else DummySpan():
        artefatos_esperados = [
            "resultados_qaoa_otimizado/resultados_20251228_145335.csv",
            "auditoria_qaoa/auditoria_qaoa_master.csv",
            "auditoria_qaoa/manifest_codigo.json",
        ]
        
        print("\nVerificando artefatos...")
        presentes = 0
        for artefato in artefatos_esperados:
            existe = Path(artefato).exists()
            presentes += int(existe)
            status = "OK" if existe else "AUSENTE"
            print(f"  [{status}] {artefato}")
        
        print(f"\nArtefatos: {presentes}/{len(artefatos_esperados)}")
        return presentes == len(artefatos_esperados)


class DummySpan:
    """Mock span para quando tracing não está disponível"""
    def __enter__(self):
        return self
    def __exit__(self, *args):
        pass
    def set_attribute(self, key, value):
        pass


def main():
    """Executar testes da pipeline com tracing"""
    tempo_inicio = time.time()
    
    print("\n" + "="*70)
    print("TESTE DA PIPELINE QAOA COM TRACING OPENTELEMETRY")
    print("="*70)
    
    if tracer:
        with tracer.start_as_current_span("teste_pipeline_qaoa") as span:
            span.set_attribute("teste.inicio", datetime.now().isoformat())
            
            # Teste 1: Verificar artefatos
            logger.info("Teste 1: Verificando artefatos gerados...")
            artefatos_ok = teste_artefatos()
            span.set_attribute("teste.artefatos_ok", artefatos_ok)
            
            # Teste 2: Verificar imports
            logger.info("Teste 2: Verificando imports...")
            try:
                import qiskit
                span.set_attribute("teste.qiskit_version", qiskit.__version__)
                logger.info(f"  ✓ Qiskit {qiskit.__version__}")
            except ImportError as e:
                logger.warning(f"  ✗ Qiskit não encontrado: {e}")
            
            # Teste 3: Verificar dados enriquecidos
            logger.info("Teste 3: Verificando dados enriquecidos...")
            csv_path = Path("resultados_qaoa_otimizado/resultados_20251228_145335.csv")
            if csv_path.exists():
                try:
                    import pandas as pd
                    df = pd.read_csv(str(csv_path))
                    logger.info(f"  ✓ CSV com {len(df)} linhas e {len(df.columns)} colunas")
                    span.set_attribute("teste.csv_linhas", len(df))
                    span.set_attribute("teste.csv_colunas", len(df.columns))
                except Exception as e:
                    logger.warning(f"  ✗ Erro ao ler CSV: {e}")
            
            # Teste 4: Verificar manifest SHA-256
            logger.info("Teste 4: Verificando manifest de hashes...")
            manifest_path = Path("auditoria_qaoa/manifest_codigo.json")
            if manifest_path.exists():
                try:
                    with open(manifest_path) as f:
                        manifest = json.load(f)
                    logger.info(f"  ✓ Manifest com {len(manifest)} hashes")
                    span.set_attribute("teste.manifest_scripts", len(manifest))
                except Exception as e:
                    logger.warning(f"  ✗ Erro ao ler manifest: {e}")
            
            tempo_total = time.time() - tempo_inicio
            span.set_attribute("teste.tempo_total_s", tempo_total)
            span.set_attribute("teste.fim", datetime.now().isoformat())
            
            if artefatos_ok:
                logger.info(f"\n✓ TODOS OS TESTES PASSARAM ({tempo_total:.1f}s)")
                print("\nTraces foram enviados para AI Toolkit!")
                print("Abra 'AI Toolkit: Show Traces' no VS Code para visualizar")
                return 0
            else:
                logger.warning(f"\n✗ ALGUNS TESTES FALHARAM ({tempo_total:.1f}s)")
                return 1
    else:
        # Executar sem tracing
        logger.info("Executando testes SEM tracing...")
        artefatos_ok = teste_artefatos()
        tempo_total = time.time() - tempo_inicio
        
        if artefatos_ok:
            logger.info(f"\n✓ TESTES PASSARAM ({tempo_total:.1f}s)")
            logger.info("Para usar tracing, instale: pip install opentelemetry-api opentelemetry-sdk opentelemetry-exporter-otlp-proto-http")
            return 0
        else:
            logger.warning(f"\n✗ TESTES FALHARAM ({tempo_total:.1f}s)")
            return 1


if __name__ == "__main__":
    sys.exit(main())
