# ============================================================================
# audit.py
# Auditoria de integridade com SHA-256 para rastreabilidade Qualis A1
# ============================================================================
import hashlib
import pathlib
import json
import os
import datetime
from typing import Dict, Tuple

def hash_file(filepath: str) -> str:
    """Calcula SHA-256 de um arquivo."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def hash_all_files(root_dir: str = ".",
                   outfile: str = "checksums.sha256",
                   exclude_patterns: list = None) -> str:
    """
    Calcula SHA-256 de todos os arquivos em um diretório.

    Args:
        root_dir: Diretório raiz
        outfile: Arquivo de saída com checksums
        exclude_patterns: Padrões a excluir (ex: ['*.log', 'checksums*', '__pycache__'])

    Returns:
        Caminho do arquivo com checksums
    """
    if exclude_patterns is None:
        exclude_patterns = ['*.log', 'checksums*', '__pycache__', '.git', '*.pyc']

    checksums = {}
    root_path = pathlib.Path(root_dir)

    print(f"\n📊 AUDITANDO INTEGRIDADE ({root_dir})...")
    count = 0

    for filepath in sorted(root_path.rglob("*")):
        if not filepath.is_file():
            continue

        # Verificar exclusões
        should_exclude = False
        for pattern in exclude_patterns:
            if pattern.startswith('.'):
                # padrão literal (ex: '.git')
                if pattern in filepath.parts:
                    should_exclude = True
                    break
            elif '*' in pattern:
                # wildcard (ex: '*.log')
                from fnmatch import fnmatch
                if fnmatch(filepath.name, pattern):
                    should_exclude = True
                    break

        if should_exclude:
            continue

        try:
            file_hash = hash_file(str(filepath))
            rel_path = filepath.relative_to(root_path)
            checksums[str(rel_path)] = file_hash
            count += 1
            if count % 10 == 0:
                print(f"   {count} arquivos processados...")
        except Exception as e:
            print(f"   ⚠️  Erro ao processar {filepath}: {e}")

    # Salvar checksums
    with open(outfile, 'w') as f:
        for filepath in sorted(checksums.keys()):
            f.write(f"{checksums[filepath]}  {filepath}\n")

    print(f"✅ {count} arquivos auditados")
    print(f"   Checksums salvos em: {outfile}\n")

    return outfile


def verify_checksums(checksum_file: str) -> Tuple[bool, Dict]:
    """
    Verifica integridade dos arquivos contra checksums salvos.

    Returns:
        Tupla (all_valid: bool, results: dict com detalhes)
    """
    results = {
        "total": 0,
        "valid": 0,
        "missing": 0,
        "modified": [],
        "timestamp": datetime.datetime.utcnow().isoformat() + "Z"
    }

    if not os.path.exists(checksum_file):
        print(f"❌ Arquivo de checksums não encontrado: {checksum_file}")
        return False, results

    print(f"\n🔍 VERIFICANDO INTEGRIDADE ({checksum_file})...")

    with open(checksum_file, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if not line.strip():
            continue

        parts = line.strip().split("  ")
        if len(parts) != 2:
            continue

        expected_hash, filepath = parts
        results["total"] += 1

        if not os.path.exists(filepath):
            results["missing"] += 1
            print(f"   ❌ FALTANDO: {filepath}")
        else:
            actual_hash = hash_file(filepath)
            if actual_hash == expected_hash:
                results["valid"] += 1
            else:
                results["modified"].append(filepath)
                print(f"   ⚠️  MODIFICADO: {filepath}")

    all_valid = len(results["modified"]) == 0 and results["missing"] == 0

    print(f"\n✅ Válidos: {results['valid']}/{results['total']}")
    if results["missing"] > 0:
        print(f"❌ Faltando: {results['missing']}")
    if results["modified"]:
        print(f"⚠️  Modificados: {len(results['modified'])}")

    return all_valid, results


def audit_report(results_dir: str = ".") -> str:
    """Gera relatório de auditoria JSON."""
    hash_all_files(results_dir, f"{results_dir}/checksums.sha256")
    all_valid, results = verify_checksums(f"{results_dir}/checksums.sha256")

    audit_data = {
        "audit_timestamp": datetime.datetime.utcnow().isoformat() + "Z",
        "directory": results_dir,
        "integrity_status": "VALID" if all_valid else "COMPROMISED",
        "statistics": results,
        "recommendation": "SAFE FOR PUBLICATION" if all_valid else "INVESTIGATE MODIFICATIONS"
    }

    report_file = f"{results_dir}/audit_report.json"
    with open(report_file, 'w') as f:
        json.dump(audit_data, f, indent=2)

    print(f"\n📋 Relatório de auditoria: {report_file}")
    return report_file


if __name__ == "__main__":
    # Exemplo de uso
    hash_all_files(".", "checksums_test.sha256")
    verify_checksums("checksums_test.sha256")
