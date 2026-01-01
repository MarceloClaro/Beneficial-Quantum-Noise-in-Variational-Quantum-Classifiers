#!/usr/bin/env python3
"""
Enriquece CSV de resultados QAOA com metadados do run e métricas derivadas.
- Lê o JSON de resumo e o CSV do mesmo timestamp
- Acrescenta colunas: Run Timestamp, Run ID, Qubits, P-layers, Shots, Max Iter, Seed
- Acrescenta métricas: Energia Normalizada (%), Melhora vs Sem Ruído (%), Classificação, AUE, TREX
- Acrescenta uma linha de resumo com tempo total e nº de experimentos
- Inclui hashes SHA-256 do manifest de código, se disponível
"""
import re
import sys
import json
from pathlib import Path
from datetime import datetime
import pandas as pd

try:
    from calculador_hashes_qaoa import CalculadorHashesQAOA
    HASHES_AVAILABLE = True
except ImportError:
    HASHES_AVAILABLE = False

BASE_DIR = Path(__file__).parent
DIR_RESULTADOS = BASE_DIR / "resultados_qaoa_otimizado"

# Utilidades
TS_PATTERN = re.compile(r"(\d{8}_\d{6})")

def encontrar_pares_csv_json(diretorio: Path):
    pares = []
    csvs = list(diretorio.glob("resultados_*.csv"))
    jsons = list(diretorio.glob("resumo_*.json"))
    mapa_json = {TS_PATTERN.search(p.name).group(1): p for p in jsons if TS_PATTERN.search(p.name)}
    for c in csvs:
        m = TS_PATTERN.search(c.name)
        if not m:
            continue
        ts = m.group(1)
        j = mapa_json.get(ts)
        if j and j.exists():
            pares.append((c, j, ts))
    return sorted(pares)


def enriquecer_arquivo(csv_path: Path, json_path: Path, run_ts_str: str):
    # Carregar dados
    df = pd.read_csv(csv_path)
    with open(json_path, "r", encoding="utf-8") as f:
        resumo = json.load(f)

    # Metadados
    run_timestamp_iso = resumo.get("timestamp")
    cfg = resumo.get("configuracao", {})
    tempo_total = resumo.get("tempo_total_segundos")
    num_exps = resumo.get("num_experimentos")

    # Run ID a partir do timestamp do arquivo
    run_id = f"run_{run_ts_str}"

    # Energia de referência (Sem Ruído)
    energia_sem = None
    if not df.empty:
        # O CSV gerado tem coluna "Experimento" e "Energia Final"
        sem_ruido_row = df.loc[df["Experimento"].str.contains("Sem Ru", na=False)]
        if not sem_ruido_row.empty:
            energia_sem = float(sem_ruido_row.iloc[0]["Energia Final"]) if "Energia Final" in df.columns else None

    energia_max = float(df["Energia Final"].max()) if "Energia Final" in df.columns and not df.empty else None

    # Funções auxiliares
    def energia_normalizada(val):
        if energia_max and energia_max > 0:
            return (val / energia_max) * 100.0
        return None

    def melhora_vs_sem(val):
        if energia_sem and energia_sem > 0:
            return ((val - energia_sem) / energia_sem) * 100.0
        return 0.0 if energia_sem else None

    def classificacao(delta):
        if delta is None:
            return None
        if delta > 1.0:
            return "Ruído benéfico"
        if delta < -1.0:
            return "prejudicial"
        return "negligenciável"

    def aue_score(val_norm, delta):
        # Heurística simples: pondera ganho pelo nível de energia normalizada
        if val_norm is None or delta is None:
            return None
        return max(0.0, delta) * (val_norm / 100.0)

    def trex_index(val, tempo):
        # Eficiência energia/tempo
        if val is None or tempo is None or tempo == 0:
            return None
        return val / tempo

    # Acrescentar colunas
    df["Run Timestamp"] = run_timestamp_iso
    df["Run ID"] = run_id
    df["Qubits"] = cfg.get("n_qubits")
    df["P-layers"] = cfg.get("p_layers")
    df["Shots"] = cfg.get("shots")
    df["Max Iter"] = cfg.get("max_iter")
    df["Seed"] = cfg.get("seed")

    # Métricas derivadas
    if "Energia Final" in df.columns:
        df["Energia Normalizada (%)"] = df["Energia Final"].apply(energia_normalizada)
        df["Melhora vs Sem Ruído (%)"] = df["Energia Final"].apply(melhora_vs_sem)
        df["Classificação"] = df["Melhora vs Sem Ruído (%)"].apply(classificacao)
        # AUE: max(0, delta) * (energia_normalizada/100)
        aue_delta = df["Melhora vs Sem Ruído (%)"].fillna(0).clip(lower=0)
        aue_norm = df["Energia Normalizada (%)"].fillna(0) / 100.0
        df["AUE"] = aue_delta * aue_norm
        # TREX: Energia Final / Tempo (s)
        if "Tempo (s)" in df.columns:
            tempo = df["Tempo (s)"].replace(0, pd.NA)
            df["TREX"] = df["Energia Final"] / tempo

    # Linha de resumo
    resumo_row = {
        "Experimento": "Resumo do Run",
        "Tipo Ruído": None,
        "Nível Ruído": None,
        "Energia Final": None,
        "Iterações": None,
        "Tempo (s)": None,
        "Status": None,
        "Run Timestamp": run_timestamp_iso,
        "Run ID": run_id,
        "Qubits": cfg.get("n_qubits"),
        "P-layers": cfg.get("p_layers"),
        "Shots": cfg.get("shots"),
        "Max Iter": cfg.get("max_iter"),
        "Seed": cfg.get("seed"),
        "Energia Normalizada (%)": None,
        "Melhora vs Sem Ruído (%)": None,
        "Classificação": None,
        "AUE": None,
        "TREX": None,
        "Tempo total (s)": tempo_total,
        "Número de Experimentos": num_exps,
        "Observações": "Resultados enriquecidos com metadados e métricas derivadas",
        "Hashes SHA-256 (manifest)": None  # Preenchido se disponível
    }
    # Adicionar hashes de código se disponível
    if HASHES_AVAILABLE:
        calc = CalculadorHashesQAOA()
        manifest = calc.gerar_manifest()
        # Serializar manifest como string JSON compacta para incluir na célula
        resumo_row["Hashes SHA-256 (manifest)"] = json.dumps(manifest)
    # Inserir colunas de resumo no DataFrame (garantindo que existam)
    for k in ["Tempo total (s)", "Número de Experimentos", "Observações", "Hashes SHA-256 (manifest)"]:
        if k not in df.columns:
            df[k] = None
    df = pd.concat([df, pd.DataFrame([resumo_row])], ignore_index=True)

    # Salvar sobrescrevendo o CSV original
    df.to_csv(csv_path, index=False)
    print(f"✅ CSV enriquecido: {csv_path}")


def main():
    pares = encontrar_pares_csv_json(DIR_RESULTADOS)
    if not pares:
        print("Nenhum par CSV/JSON encontrado em resultados_qaoa_otimizado.")
        return 1

    # Enriquecer todos os pares encontrados; no seu caso, inclui o timestamp 20251228_145335
    for csv_path, json_path, ts in pares:
        try:
            enriquecer_arquivo(csv_path, json_path, ts)
        except Exception as e:
            print(f"Erro enriquecendo {csv_path.name}: {e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
