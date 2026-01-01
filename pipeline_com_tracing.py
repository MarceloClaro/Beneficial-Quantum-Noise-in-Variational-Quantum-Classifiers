"""
pipeline_com_tracing.py
Wrapper da pipeline QAOA com OpenTelemetry Tracing integrado
Rastreia cada etapa: experimento → enriquecimento → consolidação → validação
"""

import subprocess
import sys
import json
import time
from datetime import datetime
from pathlib import Path
from tracing_setup import setup_tracing, get_tracer

# Configurar tracing global
setup_tracing("qaoa-pipeline-auditoria")
tracer = get_tracer("pipeline_com_tracing")


def executar_etapa(script: str, descricao: str, timeout: int = 600):
    """
    Executar etapa da pipeline com rastreamento

    Args:
        script: Nome do script Python
        descricao: Descrição da etapa
        timeout: Tempo máximo em segundos

    Returns:
        tuple: (sucesso: bool, tempo_s: float, saída: str)
    """
    with tracer.start_as_current_span(f"etapa_{script.split('.')[0]}") as span:
        span.set_attribute("etapa.script", script)
        span.set_attribute("etapa.descricao", descricao)
        span.set_attribute("etapa.timestamp", datetime.now().isoformat())

        print(f"\n📊 [{descricao}]")
        print(f"{'='*70}")

        tempo_inicio = time.time()
        try:
            # Executar script
            resultado = subprocess.run(
                [sys.executable, script],
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd="."
            )

            tempo_decorrido = time.time() - tempo_inicio

            # Rastrear resultado
            span.set_attribute("etapa.status", "sucesso" if resultado.returncode == 0 else "falha")
            span.set_attribute("etapa.tempo_s", tempo_decorrido)
            span.set_attribute("etapa.exit_code", resultado.returncode)

            if resultado.returncode == 0:
                print(f"✅ {descricao} — SUCESSO ({tempo_decorrido:.1f}s)")
                return (True, tempo_decorrido, resultado.stdout)
            else:
                print(f"❌ {descricao} — FALHA ({tempo_decorrido:.1f}s)")
                print(f"Erro: {resultado.stderr[:200]}")
                span.set_attribute("etapa.erro", resultado.stderr[:500])
                return (False, tempo_decorrido, resultado.stderr)

        except subprocess.TimeoutExpired:
            tempo_decorrido = time.time() - tempo_inicio
            print(f"⏱️  {descricao} — TIMEOUT ({tempo_decorrido:.1f}s > {timeout}s)")
            span.set_attribute("etapa.status", "timeout")
            span.set_attribute("etapa.tempo_s", tempo_decorrido)
            return (False, tempo_decorrido, "Timeout")

        except Exception as e:
            tempo_decorrido = time.time() - tempo_inicio
            print(f"❌ {descricao} — EXCEÇÃO: {e}")
            span.set_attribute("etapa.status", "exceção")
            span.set_attribute("etapa.erro", str(e))
            span.set_attribute("etapa.tempo_s", tempo_decorrido)
            return (False, tempo_decorrido, str(e))



def verificar_artefatos():
    """Verificar se todos os artefatos foram gerados"""
    with tracer.start_as_current_span("verificacao_artefatos") as span:
        arquivos_esperados = [
            "resultados_qaoa_otimizado/resultados_20251228_145335.csv",
            "auditoria_qaoa/auditoria_qaoa_master.csv",
            "auditoria_qaoa/auditoria_qaoa_energia.png",
            "auditoria_qaoa/manifest_codigo.json",
        ]

        presentes = 0
        for arquivo in arquivos_esperados:
            existe = Path(arquivo).exists()
            presentes += int(existe)
            span.set_attribute(f"artefato.{arquivo}", existe)
            status = "✅" if existe else "❌"
            print(f"  {status} {arquivo}")

        total = len(arquivos_esperados)
        print(f"\n📊 Artefatos: {presentes}/{total}")
        span.set_attribute("artefatos.presentes", presentes)
        span.set_attribute("artefatos.total", total)

        return presentes == total


def main():
    """Executar pipeline QAOA com tracing completo"""
    tempo_total_inicio = time.time()

    with tracer.start_as_current_span("pipeline_qaoa_completa") as span_principal:
        span_principal.set_attribute("pipeline.inicio", datetime.now().isoformat())

        print("\n" + "🚀 PIPELINE QAOA COM TRACING".center(70))
        print("="*70)

        resultados = []

        # Etapa 1: Experimento
        sucesso1, tempo1, saida1 = executar_etapa(
            "experimento_qaoa_otimizado.py",
            "Etapa 1: Experimento QAOA (6 qubits × 4 ruídos)"
        )
        resultados.append(("Experimento", sucesso1, tempo1))

        if not sucesso1:
            print("\n❌ Pipeline parada na Etapa 1")
            span_principal.set_attribute("pipeline.status", "falha_etapa_1")
            return 1

        # Etapa 2: Enriquecimento
        sucesso2, tempo2, saida2 = executar_etapa(
            "enriquecer_resultados_qaoa.py",
            "Etapa 2: Enriquecimento (Metadados + Métricas)"
        )
        resultados.append(("Enriquecimento", sucesso2, tempo2))

        if not sucesso2:
            print("\n❌ Pipeline parada na Etapa 2")
            span_principal.set_attribute("pipeline.status", "falha_etapa_2")
            return 1

        # Etapa 3: Consolidação
        sucesso3, tempo3, saida3 = executar_etapa(
            "auditoria_qaoa_resultados.py",
            "Etapa 3: Consolidação (Master CSV + Gráficos)"
        )
        resultados.append(("Consolidação", sucesso3, tempo3))

        if not sucesso3:
            print("\n❌ Pipeline parada na Etapa 3")
            span_principal.set_attribute("pipeline.status", "falha_etapa_3")
            return 1

        # Etapa 4: Validação
        sucesso4, tempo4, saida4 = executar_etapa(
            "validar_auditoria_qaoa.py",
            "Etapa 4: Validação QUALIS A1"
        )
        resultados.append(("Validação", sucesso4, tempo4))

        # Etapa 5: Rastreabilidade
        sucesso5, tempo5, saida5 = executar_etapa(
            "calculador_hashes_qaoa.py",
            "Etapa 5: Rastreabilidade (SHA-256)"
        )
        resultados.append(("Rastreabilidade", sucesso5, tempo5))

        # Verificar artefatos
        print("\n" + "📋 Verificando Artefatos".center(70))
        artefatos_ok = verificar_artefatos()

        # Resumo com tracing
        tempo_total = time.time() - tempo_total_inicio

        print("\n" + "="*70)
        print("RESUMO COM TRACING".center(70))
        print("="*70)

        for etapa, sucesso, tempo in resultados:
            status = "✅ PASS" if sucesso else "❌ FAIL"
            print(f"{etapa:.<50} {tempo:>6.1f}s {status}")

        print(f"\n{'Tempo total':.<50} {tempo_total:>6.1f}s")

        # Atributos finais para span
        span_principal.set_attribute("pipeline.fim", datetime.now().isoformat())
        span_principal.set_attribute("pipeline.tempo_total_s", tempo_total)
        span_principal.set_attribute("pipeline.artefatos_ok", artefatos_ok)
        span_principal.set_attribute("pipeline.etapas_pass", sum(1 for _, s, _ in resultados if s))

        if all(s for _, s, _ in resultados) and artefatos_ok:
            print("\n✅ PIPELINE CONCLUÍDA COM SUCESSO!")
            print("💡 Traces enviados para AI Toolkit")
            span_principal.set_attribute("pipeline.status", "sucesso")
            return 0
        else:
            print("\n❌ PIPELINE COM FALHAS")
            span_principal.set_attribute("pipeline.status", "falha")
            return 1


if __name__ == "__main__":
    sys.exit(main())
