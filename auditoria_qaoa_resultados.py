#!/usr/bin/env python3
"""
Auditoria de Resultados QAOA
Consolida CSV/JSON de diferentes execuções e gera gráficos para QUALIS A1.
Inclui hashes SHA-256 para rastreabilidade de versão de código.
"""

import os
import sys
import json
import glob
import platform
from pathlib import Path
from datetime import datetime

from typing import List, Dict, Any

import pandas as pd

# Importar calculador de hashes (se disponível)
try:
    from calculador_hashes_qaoa import CalculadorHashesQAOA
    HASHES_AVAILABLE = True
except ImportError:
    HASHES_AVAILABLE = False

# Visualização
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    import plotly.express as px
    PLOTLY_OK = True
except Exception:
    px = None
    PLOTLY_OK = False

BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / "auditoria_qaoa"
OUTPUT_DIR.mkdir(exist_ok=True)

RESULT_DIRS = [
    BASE_DIR / "resultados_qaoa_experimento_completo",
    BASE_DIR / "resultados_qaoa_otimizado",
]


def coletar_jsons(diretorio: Path) -> List[Path]:
    return sorted(diretorio.glob("resumo_*.json"))


def coletar_csvs(diretorio: Path) -> List[Path]:
    return sorted(diretorio.glob("resultados_*.csv"))


def extrair_registros_json(json_path: Path) -> List[Dict[str, Any]]:
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    registros = []
    timestamp = data.get("timestamp")
    configuracao = data.get("configuracao", {})
    for exp in data.get("experimentos", []):
        registro = {
            "timestamp": timestamp,
            "fonte": str(json_path.parent.name),
            "arquivo": json_path.name,
            "Experimento": exp.get("nome") or exp.get("Experimento"),
            "Tipo Ruído": exp.get("Tipo Ruído") or None,
            "Nível Ruído": exp.get("Nível Ruído") or None,
            "Energia Final": exp.get("energia") or exp.get("Energia Final"),
            "Iterações": exp.get("Iterações"),
            "Tempo (s)": exp.get("tempo") or exp.get("Tempo (s)"),
            "Status": exp.get("status"),
            "Erro": exp.get("erro"),
            # Configuracao
            "n_qubits": configuracao.get("n_qubits"),
            "p_layers": configuracao.get("p_layers"),
            "shots": configuracao.get("shots"),
            "max_iter": configuracao.get("max_iter"),
            "seed": configuracao.get("seed"),
        }
        registros.append(registro)
    return registros


def extrair_registros_csv(csv_path: Path) -> List[Dict[str, Any]]:
    df = pd.read_csv(csv_path)
    registros = []
    for _, row in df.iterrows():
        registros.append({
            "timestamp": None,
            "fonte": str(csv_path.parent.name),
            "arquivo": csv_path.name,
            "Experimento": row.get("nome") or row.get("Experimento"),
            "Tipo Ruído": row.get("Tipo Ruído") if "Tipo Ruído" in df.columns else None,
            "Nível Ruído": row.get("Nível Ruído") if "Nível Ruído" in df.columns else None,
            "Energia Final": row.get("Energia Final") if "Energia Final" in df.columns else None,
            "Iterações": row.get("Iterações") if "Iterações" in df.columns else None,
            "Tempo (s)": row.get("Tempo (s)") if "Tempo (s)" in df.columns else None,
            "Status": row.get("status") or row.get("Status"),
            "Erro": row.get("erro") or row.get("Erro"),
            "n_qubits": None,
            "p_layers": None,
            "shots": None,
            "max_iter": None,
            "seed": None,
            # Campos enriquecidos, se presentes
            "Run Timestamp": row.get("Run Timestamp") if "Run Timestamp" in df.columns else None,
            "Run ID": row.get("Run ID") if "Run ID" in df.columns else None,
            "Qubits": row.get("Qubits") if "Qubits" in df.columns else None,
            "P-layers": row.get("P-layers") if "P-layers" in df.columns else None,
            "Shots": row.get("Shots") if "Shots" in df.columns else None,
            "Max Iter": row.get("Max Iter") if "Max Iter" in df.columns else None,
            "Seed": row.get("Seed") if "Seed" in df.columns else None,
            "Energia Normalizada (%)": (
                row.get("Energia Normalizada (%)")
                if "Energia Normalizada (%)" in df.columns else None
            ),
            "Melhora vs Sem Ruído (%)": (
                row.get("Melhora vs Sem Ruído (%)")
                if "Melhora vs Sem Ruído (%)" in df.columns else None
            ),
            "Classificação": row.get("Classificação") if "Classificação" in df.columns else None,
            "AUE": row.get("AUE") if "AUE" in df.columns else None,
            "TREX": row.get("TREX") if "TREX" in df.columns else None,
            "Tempo total (s)": row.get("Tempo total (s)") if "Tempo total (s)" in df.columns else None,
            "Número de Experimentos": (
                row.get("Número de Experimentos")
                if "Número de Experimentos" in df.columns else None
            ),
            "Observações": row.get("Observações") if "Observações" in df.columns else None,
        })
    return registros


def gerar_graficos(df: pd.DataFrame) -> None:
    # Gráfico 1: Energia por experimento (quando disponível)
    df_energy = df.dropna(subset=["Energia Final"]).copy()
    if not df_energy.empty:
        plt.figure(figsize=(10, 5))
        plt.bar(df_energy["Experimento"], df_energy["Energia Final"], color="#3b82f6")
        plt.xticks(rotation=30, ha="right")
        plt.ylabel("Energia (MaxCut)")
        plt.title("Energia por Experimento (QAOA)")
        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / "auditoria_qaoa_energia.png", dpi=300)
        plt.close()
        if PLOTLY_OK and px is not None:
            fig = px.bar(df_energy, x="Experimento", y="Energia Final", color="Tipo Ruído",
                         title="Energia por Experimento (QAOA)")
            fig.write_html(str(OUTPUT_DIR / "auditoria_qaoa_energia.html"))

    # Gráfico 2: Tempo por experimento
    df_time = df.dropna(subset=["Tempo (s)"]).copy()
    if not df_time.empty:
        plt.figure(figsize=(10, 5))
        plt.bar(df_time["Experimento"], df_time["Tempo (s)"], color="#10b981")
        plt.xticks(rotation=30, ha="right")
        plt.ylabel("Tempo (s)")
        plt.title("Tempo de Execução por Experimento (QAOA)")
        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / "auditoria_qaoa_tempo.png", dpi=300)
        plt.close()
        if PLOTLY_OK and px is not None:
            fig = px.bar(df_time, x="Experimento", y="Tempo (s)", color="Tipo Ruído",
                         title="Tempo de Execução por Experimento (QAOA)")
            fig.write_html(str(OUTPUT_DIR / "auditoria_qaoa_tempo.html"))


def salvar_ambiente() -> None:
    info = {
        "timestamp": datetime.now().isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {}
    }
    try:
        import qiskit
        import qiskit_aer
        import numpy
        import pandas
        info["packages"] = {
            "qiskit": getattr(qiskit, "__version__", None),
            "qiskit_aer": getattr(qiskit_aer, "__version__", None),
            "numpy": getattr(numpy, "__version__", None),
            "pandas": getattr(pandas, "__version__", None),
        }
    except Exception:
        pass
    with open(OUTPUT_DIR / "ambiente_execucao.json", "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2)


def main():
    todos_registros: List[Dict[str, Any]] = []

    for d in RESULT_DIRS:
        if not d.exists():
            continue
        # JSONs
        for j in coletar_jsons(d):
            try:
                todos_registros.extend(extrair_registros_json(j))
            except Exception as e:
                print(f"Erro lendo {j}: {e}")
        # CSVs
        for c in coletar_csvs(d):
            try:
                todos_registros.extend(extrair_registros_csv(c))
            except Exception as e:
                print(f"Erro lendo {c}: {e}")

    # DataFrame consolidado
    df = pd.DataFrame(todos_registros)
    if df.empty:
        print("Nenhum registro encontrado.")
        return 1

    # Ordenar e salvar
    df_sorted = df.sort_values(by=["timestamp", "Experimento"], na_position="last")
    master_csv = OUTPUT_DIR / "auditoria_qaoa_master.csv"
    df_sorted.to_csv(master_csv, index=False)
    print(f"✅ Consolidado salvo em: {master_csv}")

    # Gráficos
    gerar_graficos(df_sorted)
    print(f"✅ Gráficos salvos em: {OUTPUT_DIR}")

    # Ambiente
    salvar_ambiente()
    print(f"✅ Ambiente de execução documentado em: {OUTPUT_DIR / 'ambiente_execucao.json'}")

    # Hashes de código (rastreabilidade)
    if HASHES_AVAILABLE:
        calc = CalculadorHashesQAOA()
        manifest = calc.gerar_manifest()
        manifest_path = OUTPUT_DIR / "manifest_codigo.json"
        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2)
        print(f"✅ Hashes SHA-256 salvos em: {manifest_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
