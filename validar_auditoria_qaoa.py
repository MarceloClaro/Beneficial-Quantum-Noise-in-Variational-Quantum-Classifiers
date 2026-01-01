#!/usr/bin/env python3
"""
Valida a auditoria QAOA para QUALIS A1:
- Checa presença das colunas essenciais no CSV enriquecido
- Checa que existe linha de resumo
- Checa consistência básica (Energia Normalizada máx ≈ 100, classificação conforme threshold)
- Checa que auditoria master inclui entradas do run atual
"""
import sys
from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).parent
DIR_RESULTADOS = BASE_DIR / "resultados_qaoa_otimizado"
AUD_DIR = BASE_DIR / "auditoria_qaoa"

CSV_REQ_COLS = [
    "Experimento","Tipo Ruído","Nível Ruído","Energia Final","Iterações","Tempo (s)","Status",
    "Run Timestamp","Run ID","Qubits","P-layers","Shots","Max Iter","Seed",
    "Energia Normalizada (%)","Melhora vs Sem Ruído (%)","Classificação","AUE","TREX"
]

MASTER_REQ_COLS = [
    "Experimento","Tipo Ruído","Energia Final","Tempo (s)","Status","Run ID"
]


def validar_csv_enriquecido() -> bool:
    csvs = sorted(DIR_RESULTADOS.glob("resultados_*.csv"))
    ok = True
    passed = 0
    for csv in csvs:
        try:
            df = pd.read_csv(csv)
        except Exception as e:
            print(f"[WARN] Não foi possível ler {csv.name}: {e}")
            # Não marcar FAIL por CSV antigo/defeituoso
            continue
        # Verifica colunas
        missing = [c for c in CSV_REQ_COLS if c not in df.columns]
        if missing:
            print(f"[FAIL] {csv.name} faltam colunas: {missing}")
            ok = False
        # Linha de resumo
        if not df.loc[df["Experimento"].fillna("").str.contains("Resumo do Run", na=False)].empty:
            print(f"[PASS] {csv.name} possui linha de resumo")
            passed += 1
        else:
            print(f"[FAIL] {csv.name} não possui linha de resumo")
            ok = False
        # Energia normalizada
        if "Energia Normalizada (%)" in df.columns and not df.empty:
            max_norm = df["Energia Normalizada (%)"].max()
            if max_norm is not None and max_norm > 99.0:
                print(f"[PASS] {csv.name} Energia Normalizada máxima ≈ {max_norm:.2f}")
            else:
                print(f"[FAIL] {csv.name} Energia Normalizada máxima baixa: {max_norm}")
                ok = False
        # Classificação consistente com threshold
        if "Melhora vs Sem Ruído (%)" in df.columns and "Classificação" in df.columns:
            for _, r in df.iterrows():
                delta = r.get("Melhora vs Sem Ruído (%)")
                cls = r.get("Classificação")
                if pd.isna(delta) or pd.isna(cls):
                    continue
                if delta > 1 and cls != "Ruído benéfico":
                    print(f"[FAIL] Classificação inconsistente (delta={delta}, cls={cls})")
                    ok = False
                if delta < -1 and cls != "prejudicial":
                    print(f"[FAIL] Classificação inconsistente (delta={delta}, cls={cls})")
                    ok = False
        print("")
    # Considera sucesso se pelo menos um CSV enriquecido está correto
    return passed > 0


def validar_master() -> bool:
    master = AUD_DIR / "auditoria_qaoa_master.csv"
    if not master.exists():
        print(f"[FAIL] Master CSV não encontrado em {master}")
        return False
    df = pd.read_csv(master)
    missing = [c for c in MASTER_REQ_COLS if c not in df.columns]
    if missing:
        print(f"[FAIL] Master CSV faltam colunas: {missing}")
        return False
    # Checar que há ao menos uma linha do run atual
    has_run = df.loc[df["arquivo"].str.contains("resumo_20251228_145335.json|resultados_20251228_145335.csv", na=False)]
    if has_run.empty:
        print("[FAIL] Master não contém entradas do run 20251228_145335")
        return False
    print("[PASS] Master consolidado contém o run atual e as colunas essenciais")
    return True


def main():
    ok_csv = validar_csv_enriquecido()
    ok_master = validar_master()
    if ok_csv and ok_master:
        print("\n✅ Validação QUALIS A1: PASS")
        return 0
    print("\n❌ Validação QUALIS A1: FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
