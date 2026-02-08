#!/usr/bin/env python3
"""Validador de artefatos para banca CNPq/Qualis A1."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

REQUIRED_FILES = [
    "environment_snapshot.txt",
    "cv_comparison.csv",
    "02_validation_report.json",
    "02_validation_report.md",
    "checksums_final.sha256",
]

REQUIRED_PREFIXES = [
    "01_protocolo_pre_registrado_",
]

REQUIRED_GATES = [
    "gate_pvalue",
    "gate_permutation",
    "gate_ci_above_random",
    "gate_stability_cv_lt_10pct",
]


def find_prefixed_file(results_dir: Path, prefix: str) -> bool:
    return any(p.name.startswith(prefix) for p in results_dir.iterdir() if p.is_file())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", required=True, help="Diretório de resultados gerado pelo pipeline")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.exists() or not results_dir.is_dir():
        print(f"❌ Diretório inválido: {results_dir}")
        return 2

    ok = True

    for fname in REQUIRED_FILES:
        fpath = results_dir / fname
        if fpath.exists():
            print(f"✅ arquivo presente: {fname}")
        else:
            print(f"❌ arquivo ausente: {fname}")
            ok = False

    for prefix in REQUIRED_PREFIXES:
        if find_prefixed_file(results_dir, prefix):
            print(f"✅ arquivo com prefixo presente: {prefix}*")
        else:
            print(f"❌ arquivo com prefixo ausente: {prefix}*")
            ok = False

    report_path = results_dir / "02_validation_report.json"
    if report_path.exists():
        report = json.loads(report_path.read_text(encoding="utf-8"))
        gates = report.get("quality_gates", {})
        for gate in REQUIRED_GATES:
            gate_val = gates.get(gate, False)
            if bool(gate_val):
                print(f"✅ gate aprovado: {gate}=true")
            else:
                print(f"❌ gate não aprovado: {gate}={gate_val}")
                ok = False

    if ok:
        print("\n✅ Validação concluída: pacote apto para banca/auditoria.")
        return 0

    print("\n❌ Validação falhou: revisar itens acima.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
