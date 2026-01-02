#!/usr/bin/env python3
# =============================================================================
# DEEPCHEM INSTALLATION SCRIPT
# =============================================================================
"""
Script para instalação e verificação do DeepChem com suas dependências.

DeepChem é uma biblioteca para aprendizado profundo em química e biologia molecular.
Suporta datasets moleculares como HIV, Malária, TB para descoberta de medicamentos.

INSTALAÇÃO:
    python install_deepchem.py

REQUISITOS:
    - Python 3.8+
    - pip
    - Opcional: conda (recomendado para RDKit)

DATASETS SUPORTADOS:
    - HIV: Atividade antiviral (41,127 compostos)
    - Malária: Atividade antimalárica
    - TB: Atividade contra tuberculose
    - E mais 20+ datasets moleculares

REFERÊNCIAS:
    - DeepChem: https://deepchem.io/
    - Documentação: https://deepchem.readthedocs.io/
"""

import os
import sys
import subprocess
import importlib
from pathlib import Path
from typing import List, Tuple, Optional


class DeepChemInstaller:
    """Instalador e verificador para DeepChem."""
    
    def __init__(self, verbose: bool = True):
        self.verbose = verbose
        self.python_executable = sys.executable
        
    def log(self, message: str):
        """Imprime mensagem se verbose."""
        if self.verbose:
            print(message)
    
    def check_python_version(self) -> bool:
        """Verifica se Python é 3.8+."""
        self.log("\n[1/7] Verificando versão do Python...")
        
        version = sys.version_info
        if version.major == 3 and version.minor >= 8:
            self.log(f"✓ Python {version.major}.{version.minor}.{version.micro} OK")
            return True
        else:
            self.log(f"✗ Python {version.major}.{version.minor}.{version.micro} não suportado")
            self.log("  DeepChem requer Python 3.8 ou superior")
            return False
    
    def check_pip(self) -> bool:
        """Verifica se pip está disponível."""
        self.log("\n[2/7] Verificando pip...")
        
        try:
            result = subprocess.run(
                [self.python_executable, "-m", "pip", "--version"],
                capture_output=True,
                text=True,
                check=True
            )
            self.log(f"✓ {result.stdout.strip()}")
            return True
        except Exception as e:
            self.log(f"✗ pip não encontrado: {e}")
            return False
    
    def upgrade_pip(self) -> bool:
        """Atualiza pip para última versão."""
        self.log("\n[3/7] Atualizando pip...")
        
        try:
            subprocess.run(
                [self.python_executable, "-m", "pip", "install", "--upgrade", "pip"],
                check=True,
                capture_output=True
            )
            self.log("✓ pip atualizado")
            return True
        except Exception as e:
            self.log(f"⚠ Não foi possível atualizar pip: {e}")
            return False
    
    def install_dependencies(self) -> bool:
        """Instala dependências básicas do DeepChem."""
        self.log("\n[4/7] Instalando dependências básicas...")
        
        dependencies = [
            "numpy",
            "pandas",
            "scikit-learn",
            "scipy",
            "joblib",
        ]
        
        try:
            for dep in dependencies:
                self.log(f"  Instalando {dep}...")
                subprocess.run(
                    [self.python_executable, "-m", "pip", "install", dep],
                    check=True,
                    capture_output=True
                )
            self.log("✓ Dependências instaladas")
            return True
        except Exception as e:
            self.log(f"✗ Erro instalando dependências: {e}")
            return False
    
    def install_deepchem(self) -> bool:
        """Instala DeepChem via pip."""
        self.log("\n[5/7] Instalando DeepChem...")
        self.log("  (Isso pode demorar alguns minutos...)")
        
        try:
            # Instala DeepChem
            result = subprocess.run(
                [self.python_executable, "-m", "pip", "install", "deepchem"],
                check=True,
                capture_output=True,
                text=True
            )
            self.log("✓ DeepChem instalado")
            return True
        except subprocess.CalledProcessError as e:
            self.log(f"✗ Erro instalando DeepChem:")
            self.log(f"  {e.stderr}")
            return False
    
    def verify_deepchem(self) -> bool:
        """Verifica se DeepChem foi instalado corretamente."""
        self.log("\n[6/7] Verificando instalação do DeepChem...")
        
        try:
            import deepchem as dc
            self.log(f"✓ DeepChem {dc.__version__} importado com sucesso")
            
            # Testa imports principais
            from deepchem.molnet import load_hiv
            self.log("✓ MoleculeNet disponível")
            
            return True
        except ImportError as e:
            self.log(f"✗ Erro ao importar DeepChem: {e}")
            return False
    
    def test_datasets(self) -> Tuple[bool, List[str]]:
        """Testa carregamento dos principais datasets."""
        self.log("\n[7/7] Testando datasets moleculares...")
        
        successful_datasets = []
        
        try:
            import deepchem as dc
            
            # Teste 1: HIV dataset
            self.log("\n  Testando HIV dataset...")
            try:
                tasks, datasets, transformers = dc.molnet.load_hiv(
                    featurizer='ECFP',
                    splitter='scaffold',
                    reload=False
                )
                train, valid, test = datasets
                self.log(f"    ✓ HIV: {len(train)} amostras de treino")
                successful_datasets.append("HIV")
            except Exception as e:
                self.log(f"    ✗ HIV falhou: {e}")
            
            # Teste 2: BACE (proxy para malária - dataset similar)
            self.log("\n  Testando BACE dataset (similar a Malária)...")
            try:
                tasks, datasets, transformers = dc.molnet.load_bace_classification(
                    featurizer='ECFP',
                    splitter='scaffold',
                    reload=False
                )
                train, valid, test = datasets
                self.log(f"    ✓ BACE: {len(train)} amostras de treino")
                successful_datasets.append("BACE")
            except Exception as e:
                self.log(f"    ✗ BACE falhou: {e}")
            
            # Teste 3: Tox21 (inclui TB-related toxicity)
            self.log("\n  Testando Tox21 dataset...")
            try:
                tasks, datasets, transformers = dc.molnet.load_tox21(
                    featurizer='ECFP',
                    splitter='scaffold',
                    reload=False
                )
                train, valid, test = datasets
                self.log(f"    ✓ Tox21: {len(train)} amostras de treino")
                successful_datasets.append("Tox21")
            except Exception as e:
                self.log(f"    ✗ Tox21 falhou: {e}")
            
            return len(successful_datasets) > 0, successful_datasets
            
        except Exception as e:
            self.log(f"✗ Erro testando datasets: {e}")
            return False, []
    
    def install_rdkit_info(self):
        """Mostra informações sobre instalação do RDKit (opcional)."""
        self.log("\n" + "=" * 80)
        self.log("INSTALAÇÃO OPCIONAL: RDKit")
        self.log("=" * 80)
        self.log("""
RDKit é uma biblioteca química opcional que melhora o DeepChem.
Ela fornece funcionalidades avançadas de química computacional.

INSTALAÇÃO RECOMENDADA (via conda):
    conda install -c conda-forge rdkit

INSTALAÇÃO VIA PIP (pode não funcionar em todos os sistemas):
    pip install rdkit-pypi

NOTA: RDKit NÃO é obrigatório para usar os datasets moleculares básicos.
      O DeepChem funciona sem RDKit usando featurizers alternativos.
""")
    
    def run_full_installation(self) -> bool:
        """Executa instalação completa."""
        self.log("=" * 80)
        self.log("INSTALADOR DEEPCHEM - Framework Quantum Advanced V8")
        self.log("=" * 80)
        
        # Checagens
        if not self.check_python_version():
            return False
        
        if not self.check_pip():
            return False
        
        # Instalação
        self.upgrade_pip()
        
        if not self.install_dependencies():
            return False
        
        if not self.install_deepchem():
            return False
        
        if not self.verify_deepchem():
            return False
        
        # Testes
        success, datasets = self.test_datasets()
        
        if success:
            self.log("\n" + "=" * 80)
            self.log("✓ INSTALAÇÃO CONCLUÍDA COM SUCESSO!")
            self.log("=" * 80)
            self.log(f"\nDatasets disponíveis: {', '.join(datasets)}")
            self.log("\nPróximos passos:")
            self.log("  1. Execute: python run_framework_quantum_advanced_v8.py --dataset hiv")
            self.log("  2. Ou use: python exemplos_framework_quantum_v8.py")
            self.log("  3. Veja FRAMEWORK_QUANTUM_ADVANCED_V8_README.md para mais informações")
            
            self.install_rdkit_info()
            return True
        else:
            self.log("\n" + "=" * 80)
            self.log("⚠ INSTALAÇÃO CONCLUÍDA COM AVISOS")
            self.log("=" * 80)
            self.log("\nDeepChem está instalado mas alguns datasets falharam.")
            self.log("Você ainda pode usar o framework com datasets padrão (iris, wine, etc.)")
            return False


def main():
    """Função principal."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Instalador DeepChem para Framework Quantum Advanced V8"
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Modo silencioso (menos output)"
    )
    parser.add_argument(
        "--test-only",
        action="store_true",
        help="Apenas testa se DeepChem já está instalado"
    )
    
    args = parser.parse_args()
    
    installer = DeepChemInstaller(verbose=not args.quiet)
    
    if args.test_only:
        print("Testando instalação existente do DeepChem...\n")
        if installer.verify_deepchem():
            success, datasets = installer.test_datasets()
            if success:
                print(f"\n✓ DeepChem está instalado e funcionando")
                print(f"  Datasets disponíveis: {', '.join(datasets)}")
                return 0
            else:
                print("\n⚠ DeepChem instalado mas datasets não carregam")
                return 1
        else:
            print("\n✗ DeepChem não está instalado")
            print("  Execute: python install_deepchem.py")
            return 1
    else:
        success = installer.run_full_installation()
        return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
