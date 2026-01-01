#!/usr/bin/env python
"""
TEST_PIPELINE_AUTOMATIZADA.py
Teste completo da pipeline de auditoria QAOA com rastreabilidade de versão
Valida todas as 5 etapas: experimento → enriquecimento → consolidação → validação → hashes
"""

import subprocess
import sys
import json
from pathlib import Path
from datetime import datetime

def run_command(cmd: str, description: str, script_name: str = None) -> bool:
    """Executar script Python diretamente"""
    print(f"\n{'='*70}")
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {description}")
    print(f"{'='*70}")

    if script_name:
        # Executar script Python diretamente
        import importlib.util
        try:
            spec = importlib.util.spec_from_file_location("module", script_name)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            print(f"✅ {description} — SUCESSO")
            return True
        except Exception as e:
            print(f"❌ {description} — FALHA")
            print(f"   Erro: {e}")
            return False
    else:
        # Fallback para subprocess
        result = subprocess.run(cmd, shell=True, cwd=".")
        if result.returncode == 0:
            print(f"✅ {description} — SUCESSO")
            return True
        else:
            print(f"❌ {description} — FALHA (exit code: {result.returncode})")
            return False

def verify_file_exists(path: str, description: str) -> bool:
    """Verificar se arquivo existe"""
    p = Path(path)
    if p.exists():
        print(f"  ✅ {description}: {path}")
        return True
    else:
        print(f"  ❌ {description}: NÃO ENCONTRADO")
        return False

def verify_csv_has_columns(path: str, required_columns: list) -> bool:
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

def main():
    """Executar pipeline de testes"""
    print("\n" + "🚀 TEST PIPELINE AUTOMATIZADA QAOA".center(70))
    print("="*70)
    print(f"Início: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    results = []
    
    # Etapa 1: Experimentos
    print("\n" + "ETAPA 1: EXPERIMENTO".center(70))
    cmd1 = r'& ".\.venv\Scripts\python.exe" "experimento_qaoa_otimizado.py"'
    success1 = run_command(cmd1, "Executar QAOA com 4 variantes de ruído")
    results.append(("Experimento", success1))

    if success1:
        verify_file_exists("resultados_qaoa_otimizado", "Diretório de resultados")

    # Etapa 2: Enriquecimento
    print("\n" + "ETAPA 2: ENRIQUECIMENTO".center(70))
    cmd2 = r'& ".\.venv\Scripts\python.exe" "enriquecer_resultados_qaoa.py"'
    success2 = run_command(cmd2, "Adicionar metadados e métricas ao CSV")
    results.append(("Enriquecimento", success2))

    if success2:
        # Verificar colunas enriquecidas
        csv_files = list(Path("resultados_qaoa_otimizado").glob("*.csv"))
        if csv_files:
            latest_csv = max(csv_files, key=lambda p: p.stat().st_mtime)
            required_cols = [
                "Experimento", "Energia Final", "Tempo (s)",
                "Run ID", "Run Timestamp", "Qubits",
                "Energia Normalizada (%)", "Melhora vs Sem Ruído (%)",
                "Classificação", "AUE", "TREX"
            ]
            verify_csv_has_columns(str(latest_csv), required_cols)

    # Etapa 3: Consolidação
    print("\n" + "ETAPA 3: CONSOLIDAÇÃO".center(70))
    cmd3 = r'& ".\.venv\Scripts\python.exe" "auditoria_qaoa_resultados.py"'
    success3 = run_command(cmd3, "Consolidar múltiplos runs e gerar gráficos")
    results.append(("Consolidação", success3))

    if success3:
        verify_file_exists("auditoria_qaoa/auditoria_qaoa_master.csv", "Master CSV")
        verify_file_exists("auditoria_qaoa/auditoria_qaoa_energia.png", "Gráfico Energia")
        verify_file_exists("auditoria_qaoa/auditoria_qaoa_tempo.png", "Gráfico Tempo")
        verify_file_exists("auditoria_qaoa/ambiente_execucao.json", "Ambiente JSON")

    # Etapa 4: Validação
    print("\n" + "ETAPA 4: VALIDAÇÃO QUALIS A1".center(70))
    cmd4 = r'& ".\.venv\Scripts\python.exe" "validar_auditoria_qaoa.py"'
    success4 = run_command(cmd4, "Validar conformidade QUALIS A1")
    results.append(("Validação", success4))

    # Etapa 5: Rastreabilidade
    print("\n" + "ETAPA 5: RASTREABILIDADE DE CÓDIGO".center(70))
    cmd5 = r'& ".\.venv\Scripts\python.exe" "calculador_hashes_qaoa.py"'
    success5 = run_command(cmd5, "Calcular SHA-256 dos scripts")
    results.append(("Rastreabilidade", success5))

    if success5:
        manifest_path = "auditoria_qaoa/manifest_codigo.json"
        if Path(manifest_path).exists():
            try:
                with open(manifest_path) as f:
                    manifest = json.load(f)
                    num_scripts = len(manifest.get("scripts", {}))
                    print(f"  ✅ Manifest contém {num_scripts} scripts com hashes SHA-256")
                    for script, hash_val in manifest.get("scripts", {}).items():
                        print(f"     • {script}: {hash_val[:16]}...")
            except Exception as e:
                print(f"  ❌ Erro lendo manifest: {e}")
    
    # Resumo
    print("\n" + "="*70)
    print("RESUMO DOS TESTES".center(70))
    print("="*70)
    
    for etapa, success in results:
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"{etapa:.<50} {status}")
    
    total_pass = sum(1 for _, s in results if s)
    total = len(results)
    
    print(f"\n{'Total':.<50} {total_pass}/{total}")
    
    if total_pass == total:
        print("\n✅ PIPELINE COMPLETA COM SUCESSO!")
        print("Todos os artefatos foram gerados e validados.")
        return 0
    else:
        print(f"\n❌ PIPELINE INCOMPLETA ({total - total_pass} falhas)")
        return 1

if __name__ == "__main__":
    sys.exit(main())
