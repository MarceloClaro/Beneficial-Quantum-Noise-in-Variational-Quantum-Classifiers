#!/usr/bin/env python3
"""
Calculador de hashes SHA-256 para rastreabilidade de código.
Gera um manifest de versão dos scripts principais para cada run.
"""
import hashlib
import json
from pathlib import Path
from typing import Dict


class CalculadorHashesQAOA:
    """Calcula e gerencia hashes SHA-256 de scripts QAOA."""

    SCRIPTS_PRINCIPAIS = [
        "experimento_qaoa_otimizado.py",
        "enriquecer_resultados_qaoa.py",
        "auditoria_qaoa_resultados.py",
        "validar_auditoria_qaoa.py",
        "framework_qaoa_100qubits.py",
    ]

    def __init__(self, base_dir: Path = None):
        self.base_dir = base_dir or Path(__file__).parent

    def calcular_hash_arquivo(self, caminho: Path) -> str:
        """Calcula SHA-256 de um arquivo."""
        sha256 = hashlib.sha256()
        try:
            with open(caminho, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    sha256.update(chunk)
            return sha256.hexdigest()
        except FileNotFoundError:
            return None

    def gerar_manifest(self) -> Dict[str, str]:
        """Gera manifest com hashes de todos os scripts principais."""
        manifest = {}
        for script in self.SCRIPTS_PRINCIPAIS:
            caminho = self.base_dir / script
            hash_val = self.calcular_hash_arquivo(caminho)
            if hash_val:
                manifest[script] = hash_val
        return manifest

    def salvar_manifest(self, caminho_output: Path, manifest: Dict[str, str]) -> None:
        """Salva manifest em JSON."""
        with open(caminho_output, "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2)


def main():
    calc = CalculadorHashesQAOA()
    manifest = calc.gerar_manifest()

    print("📋 Manifest de Hashes SHA-256 (Rastreabilidade de Código)")
    print("=" * 70)
    for script, hash_val in manifest.items():
        if hash_val:
            print(f"{script:45} | {hash_val}")
        else:
            print(f"{script:45} | [NÃO ENCONTRADO]")

    # Salvar manifest
    manifest_path = Path(__file__).parent / "auditoria_qaoa" / "manifest_codigo.json"
    manifest_path.parent.mkdir(exist_ok=True)
    calc.salvar_manifest(manifest_path, manifest)
    print(f"\n✅ Manifest salvo em: {manifest_path}")


if __name__ == "__main__":
    main()
