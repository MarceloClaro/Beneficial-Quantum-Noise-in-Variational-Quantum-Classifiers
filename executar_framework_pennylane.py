"""
executar_framework_pennylane.py
Executar Framework Investigativo Completo (PennyLane) com Tracing OpenTelemetry
Versão otimizada para Windows com suporte a tracing
"""

import subprocess
import sys
import os
import json
from pathlib import Path
from datetime import datetime

# Configurar encoding UTF-8 para Windows
os.environ['PYTHONIOENCODING'] = 'utf-8'

# Tentar importar tracing
try:
    from tracing_setup import setup_tracing, get_tracer
    setup_tracing("framework_pennylane")
    tracer = get_tracer("executar_framework_pennylane")
    print("[OK] OpenTelemetry tracing inicializado")
except Exception as e:
    print(f"[AVISO] Tracing não disponível: {e}")
    tracer = None


class DummySpan:
    """Mock span para quando tracing não está disponível"""
    def __enter__(self):
        return self
    def __exit__(self, *args):
        pass
    def set_attribute(self, key, value):
        pass


def executar_framework_pennylane():
    """Executar framework investigativo completo com PennyLane"""

    with tracer.start_as_current_span("framework_pennylane_execution") if tracer else DummySpan() as span:
        if span:
            span.set_attribute("framework", "PennyLane")
            span.set_attribute("tipo", "investigativo_completo")
            span.set_attribute("inicio", datetime.now().isoformat())

        print("\n" + "="*80)
        print("FRAMEWORK INVESTIGATIVO COMPLETO - PennyLane")
        print("="*80)
        print("\nExecutando: framework_investigativo_completo.py")
        print("Componentes:")
        print("  - TREX (Tensor Rank-Entropy Extensibility)")
        print("  - AUEC (Adaptive Unified Error Correction)")
        print("  - PennyLane device: default.mixed (noise support)")
        print("  - Rastreamento com OpenTelemetry")
        print()

        try:
            # Executar script principal com utf-8 encoding
            resultado = subprocess.run(
                [sys.executable, "-c", """
import os
os.environ['PYTHONIOENCODING'] = 'utf-8'

# Importar e executar framework
from framework_investigativo_completo import main
main()
"""],
                capture_output=True,
                text=True,
                timeout=3600,
                encoding='utf-8',
                errors='replace'
            )

            # Processa output
            if resultado.stdout:
                print(resultado.stdout)

            if resultado.stderr:
                # Filtrar erros não críticos
                linhas_erro = resultado.stderr.split('\n')
                erros_criticos = [l for l in linhas_erro if any(
                    kw in l.lower() for kw in
                    ['error', 'traceback', 'exception', 'failed']
                )]

                if erros_criticos:
                    print("\n⚠️  ERROS ENCONTRADOS:")
                    for erro in erros_criticos[:10]:
                        print(f"   {erro}")

            if resultado.returncode == 0:
                print("\n✅ Framework PennyLane executado com SUCESSO")
                if span:
                    span.set_attribute("status", "sucesso")
                    span.set_attribute("exit_code", 0)
                return 0
            else:
                print(f"\n❌ Framework encerrou com código {resultado.returncode}")
                if span:
                    span.set_attribute("status", "falha")
                    span.set_attribute("exit_code", resultado.returncode)
                return 1

        except subprocess.TimeoutExpired:
            print("\n⏱️  TIMEOUT: Framework levou mais de 1 hora para executar")
            if span:
                span.set_attribute("status", "timeout")
            return 1

        except Exception as e:
            print(f"\n❌ Erro ao executar framework: {e}")
            if span:
                span.set_attribute("status", "erro")
                span.set_attribute("erro", str(e))
            return 1

        finally:
            if span:
                span.set_attribute("fim", datetime.now().isoformat())


def main():
    """Função principal"""
    print("\n🔬 EXECUTOR FRAMEWORK PENNYLANE COM TRACING")
    print("="*80)

    # Verificar dependencies
    print("\n[1/2] Verificando dependências...")
    deps_ok = True

    try:
        import pennylane
        print(f"  ✓ PennyLane {pennylane.__version__}")
    except ImportError:
        print("  ✗ PennyLane não encontrado")
        deps_ok = False

    try:
        import numpy
        print(f"  ✓ NumPy {numpy.__version__}")
    except ImportError:
        print("  ✗ NumPy não encontrado")
        deps_ok = False

    try:
        import scipy
        print(f"  ✓ SciPy {scipy.__version__}")
    except ImportError:
        print("  ✗ SciPy não encontrado")
        deps_ok = False

    if not deps_ok:
        print("\n❌ Dependências faltando!")
        return 1

    # Executar framework
    print("\n[2/2] Executando framework PennyLane...")
    retorno = executar_framework_pennylane()

    print("\n" + "="*80)
    if retorno == 0:
        print("✅ EXECUÇÃO COMPLETA")
    else:
        print("❌ EXECUÇÃO COM ERROS")

    return retorno


if __name__ == "__main__":
    sys.exit(main())
