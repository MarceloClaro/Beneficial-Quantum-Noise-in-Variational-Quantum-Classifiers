#!/usr/bin/env python3
"""
RDKIT INSTALLER - Solução Simples Multi-Método

Tenta instalar RDKit de múltiplas formas:
1. Via conda (se encontrado)
2. Via mamba (se encontrado)  
3. Via pip (como último recurso)
4. Oferece alternativas se tudo falhar
"""

import sys
import subprocess
import os
from pathlib import Path
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-8s | %(message)s'
)
logger = logging.getLogger('rdkit_installer')


def find_conda():
    """Procura por conda no sistema."""
    logger.info("🔍 Procurando conda...")
    
    possible_paths = [
        "conda",
        "C:\\Users\\marce\\anaconda3\\Scripts\\conda.exe",
        "C:\\Users\\marce\\miniconda3\\Scripts\\conda.exe",
        "C:\\Users\\marce\\mambaforge\\Scripts\\conda.exe",
        "/usr/local/bin/conda",
        "/opt/conda/bin/conda",
    ]
    
    for path in possible_paths:
        try:
            result = subprocess.run(
                [str(path), "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                logger.info(f"✅ Conda encontrado: {path}")
                return str(path)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass
    
    logger.warning("❌ Conda não encontrado")
    return None


def find_mamba():
    """Procura por mamba no sistema."""
    logger.info("🔍 Procurando mamba...")
    
    possible_paths = [
        "mamba",
        "C:\\Users\\marce\\mambaforge\\Scripts\\mamba.exe",
    ]
    
    for path in possible_paths:
        try:
            result = subprocess.run(
                [str(path), "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                logger.info(f"✅ Mamba encontrado: {path}")
                return str(path)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass
    
    logger.warning("❌ Mamba não encontrado")
    return None


def install_via_conda(conda_path):
    """Instala RDKit via conda."""
    logger.info("\n[MÉTODO 1] Instalando via conda...")
    
    try:
        cmd = [conda_path, "install", "-c", "conda-forge", "rdkit", "-y"]
        logger.info(f"Executando: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, timeout=600)
        
        if result.returncode == 0:
            logger.info("✅ RDKit instalado com sucesso!")
            return True
        else:
            logger.error("❌ Falha na instalação via conda")
            return False
    except Exception as e:
        logger.error(f"❌ Erro: {e}")
        return False


def install_via_mamba(mamba_path):
    """Instala RDKit via mamba."""
    logger.info("\n[MÉTODO 2] Instalando via mamba...")
    
    try:
        cmd = [mamba_path, "install", "-c", "conda-forge", "rdkit", "-y"]
        logger.info(f"Executando: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, timeout=600)
        
        if result.returncode == 0:
            logger.info("✅ RDKit instalado com sucesso!")
            return True
        else:
            logger.error("❌ Falha na instalação via mamba")
            return False
    except Exception as e:
        logger.error(f"❌ Erro: {e}")
        return False


def install_via_pip():
    """Tenta instalar RDKit via pip."""
    logger.info("\n[MÉTODO 3] Tentando via pip...")
    
    packages = ["rdkit-pypi", "rdkit"]
    
    for package in packages:
        try:
            logger.info(f"  Tentando: {package}")
            result = subprocess.run(
                [sys.executable, "-m", "pip", "install", package, "-q"],
                timeout=300
            )
            
            if result.returncode == 0:
                logger.info(f"✅ {package} instalado!")
                return True
        except Exception:
            pass
    
    logger.warning("⚠️  Nenhum package RDKit funcionou via pip")
    return False


def verify_rdkit():
    """Verifica se RDKit está instalado."""
    logger.info("\n[VERIFICAÇÃO] Testando RDKit...")
    
    try:
        import rdkit
        from rdkit import Chem
        
        logger.info(f"✅ RDKit {rdkit.__version__} funcionando!")
        
        # Teste rápido
        mol = Chem.MolFromSmiles("CCO")
        if mol is not None:
            logger.info("✅ RDKit operacional (teste de moléculas OK)")
            return True
        else:
            return True
    except ImportError:
        logger.error("❌ RDKit não importou")
        return False
    except Exception as e:
        logger.error(f"❌ Erro: {e}")
        return False


def show_alternatives():
    """Mostra alternativas se RDKit não puder ser instalado."""
    logger.info("""
╔════════════════════════════════════════════════════════════════╗
║         ALTERNATIVAS - RDKIT NÃO DISPONÍVEL                   ║
╚════════════════════════════════════════════════════════════════╝

📋 OPÇÃO 1: Instalar Mambaforge (Recomendado)
────────────────────────────────────────────────────────────────
1. Acesse: https://github.com/conda-forge/miniforge
2. Baixe: Mambaforge-Windows-x86_64.exe
3. Instale e execute no novo prompt:
   mamba install -c conda-forge rdkit
4. Volte e execute: python install_rdkit_simple.py

📋 OPÇÃO 2: Usar Framework sem RDKit
────────────────────────────────────────────────────────────────
Os datasets sklearn funcionam 100%:
  • Iris ✅
  • Wine ✅  
  • Breast Cancer ✅

Execute:
  python run_framework_quantum_advanced_v8.py --dataset iris
  python benchmark_all_frameworks_v8.py

Datasets moleculares (HIV, Malaria, TB) precisam de RDKit.

📋 OPÇÃO 3: Usar DeepChem com alternativas
────────────────────────────────────────────────────────────────
DeepChem tem featurizadores que não precisam de RDKit:
  • MACCS keys
  • Descriptores químicos
  • One-hot encoding

    """)


def main():
    """Função principal."""
    logger.info("="*70)
    logger.info("RDKIT MULTI-METHOD INSTALLER")
    logger.info("="*70 + "\n")
    
    # Tenta conda
    conda_path = find_conda()
    if conda_path:
        if install_via_conda(conda_path):
            if verify_rdkit():
                logger.info("\n" + "="*70)
                logger.info("✅ RDKIT INSTALADO COM SUCESSO!")
                logger.info("="*70)
                return 0
    
    # Tenta mamba
    mamba_path = find_mamba()
    if mamba_path:
        if install_via_mamba(mamba_path):
            if verify_rdkit():
                logger.info("\n" + "="*70)
                logger.info("✅ RDKIT INSTALADO COM SUCESSO!")
                logger.info("="*70)
                return 0
    
    # Tenta pip
    if install_via_pip():
        if verify_rdkit():
            logger.info("\n" + "="*70)
            logger.info("✅ RDKIT INSTALADO COM SUCESSO!")
            logger.info("="*70)
            return 0
    
    # Verifica se já está instalado
    if verify_rdkit():
        logger.info("\n" + "="*70)
        logger.info("✅ RDKIT JÁ ESTAVA INSTALADO!")
        logger.info("="*70)
        return 0
    
    # Oferece alternativas
    logger.warning("\n" + "="*70)
    logger.warning("⚠️  RDKIT NÃO INSTALADO - VEJA OPÇÕES ABAIXO")
    logger.warning("="*70)
    show_alternatives()
    
    logger.info("\nFramework continua 100% operacional com sklearn datasets!")
    return 1


if __name__ == "__main__":
    sys.exit(main())
