#!/usr/bin/env python
"""
TEST_PIPELINE_VERIFICACAO.py
Teste simples e robusto da pipeline de auditoria QAOA
Verifica se todos os artefatos existem e validações passam
"""

import json
import sys
from pathlib import Path
from datetime import datetime

def check_file_exists(path: str, description: str) -> bool:
    """Verificar se arquivo existe"""
    p = Path(path)
    if p.exists():
        size_kb = p.stat().st_size / 1024
        print(f"  ✅ {description}: {path} ({size_kb:.1f}KB)")
        return True
    else:
        print(f"  ❌ {description}: NÃO ENCONTRADO")
        return False

def check_csv_has_columns(path: str, required_columns: list) -> bool:
    """Verificar se CSV tem colunas necessárias"""
    try:
        import pandas as pd
        df = pd.read_csv(path, nrows=1)
        missing = [col for col in required_columns if col not in df.columns]
        if not missing:
            print(f"  ✅ CSV tem todas as {len(required_columns)} colunas requeridas")
            return True
        else:
            print(f"  ❌ CSV está faltando colunas: {missing}")
            return False
    except Exception as e:
        print(f"  ❌ Erro lendo CSV: {e}")
        return False

def check_json_structure(path: str) -> bool:
    """Verificar estrutura JSON"""
    try:
        with open(path) as f:
            data = json.load(f)
        # Aceita ambas estruturas: {"scripts": {...}} ou {...} diretamente
        scripts_dict = data.get("scripts", data)
        if len(scripts_dict) >= 5:
            print(f"  ✅ JSON manifest com {len(scripts_dict)} scripts")
            return True
        else:
            print("  ❌ JSON inválido ou incompleto")
            return False
    except Exception as e:
        print(f"  ❌ Erro lendo JSON: {e}")
        return False

def main():
    """Verificar pipeline de testes"""
    print("\n" + "✅ VERIFICAÇÃO PIPELINE QAOA".center(70))
    print("="*70)
    print(f"Data: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    results = []

    # Checklist 1: Arquivos de Resultado
    print("\n" + "1️⃣  ARQUIVOS DE RESULTADO".center(70))
    csv_files = list(Path("resultados_qaoa_otimizado").glob("*.csv"))
    if csv_files:
        latest_csv = max(csv_files, key=lambda p: p.stat().st_mtime)
        check1 = check_file_exists(str(latest_csv), "CSV mais recente")
        results.append(("CSV Resultado", check1))
    else:
        print("  ❌ Nenhum CSV encontrado em resultados_qaoa_otimizado/")
        results.append(("CSV Resultado", False))

    # Checklist 2: Arquivos de Auditoria
    print("\n" + "2️⃣  ARQUIVOS DE AUDITORIA".center(70))
    check2a = check_file_exists("auditoria_qaoa/auditoria_qaoa_master.csv", "Master CSV")
    check2b = check_file_exists("auditoria_qaoa/auditoria_qaoa_energia.png", "Gráfico Energia")
    check2c = check_file_exists("auditoria_qaoa/auditoria_qaoa_tempo.png", "Gráfico Tempo")
    check2d = check_file_exists("auditoria_qaoa/ambiente_execucao.json", "Ambiente JSON")
    results.append(("Auditoria Consolidada", all([check2a, check2b, check2c, check2d])))

    # Checklist 3: Rastreabilidade
    print("\n" + "3️⃣  RASTREABILIDADE DE CÓDIGO".center(70))
    check3 = check_file_exists("auditoria_qaoa/manifest_codigo.json", "Manifest de hashes")
    if check3:
        check3_valid = check_json_structure("auditoria_qaoa/manifest_codigo.json")
        results.append(("Manifest SHA-256", check3 and check3_valid))
    else:
        results.append(("Manifest SHA-256", False))

    # Checklist 4: Documentação
    print("\n" + "4️⃣  DOCUMENTAÇÃO".center(70))
    check4a = check_file_exists("auditoria_qaoa/README_AUDITORIA_QAOA.md", "README Auditoria")
    check4b = check_file_exists("PIPELINE_AUTOMATIZADA_QAOA.md", "Pipeline Guide")
    check4c = check_file_exists("VISAO_GERAL_PIPELINE.md", "Visão Geral")
    results.append(("Documentação", all([check4a, check4b, check4c])))

    # Checklist 5: Colunas Enriquecidas
    print("\n" + "5️⃣  ENRIQUECIMENTO DO CSV".center(70))
    if csv_files:
        required_cols = [
            "Experimento", "Energia Final", "Tempo (s)",
            "Run ID", "Run Timestamp", "Qubits", "P-layers",
            "Energia Normalizada (%)", "Melhora vs Sem Ruído (%)",
            "Classificação", "AUE", "TREX"
        ]
        check5 = check_csv_has_columns(str(latest_csv), required_cols)
        results.append(("Colunas Enriquecidas", check5))
    else:
        results.append(("Colunas Enriquecidas", False))

    # Resumo
    print("\n" + "="*70)
    print("RESUMO DA VERIFICAÇÃO".center(70))
    print("="*70)

    for check, success in results:
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"{check:.<50} {status}")

    total_pass = sum(1 for _, s in results if s)
    total = len(results)

    print(f"\n{'Total':.<50} {total_pass}/{total}")

    if total_pass == total:
        print("\n✅ PIPELINE FUNCIONANDO CORRETAMENTE!")
        print("Todos os artefatos foram gerados com sucesso.")
        print("\nPróximos passos:")
        print("  1. Revisar auditoria_qaoa/auditoria_qaoa_master.csv")
        print("  2. Abrir auditoria_qaoa/auditoria_qaoa_energia.html (interativo)")
        print("  3. Verificar manifest_codigo.json para hashes SHA-256")
        print("  4. Consultar PIPELINE_AUTOMATIZADA_QAOA.md para detalhes")
        return 0
    else:
        print(f"\n⚠️  PIPELINE INCOMPLETA ({total - total_pass}/{total} falhas)")
        print("Etapas faltantes:")
        for check, success in results:
            if not success:
                print(f"  • {check}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
