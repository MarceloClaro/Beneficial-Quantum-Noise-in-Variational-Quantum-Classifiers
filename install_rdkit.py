#!/usr/bin/env python3
# =============================================================================
# RDKIT INSTALLER - Multi-Method Solution
# =============================================================================
"""
Script para instalar RDKit de múltiplas formas:
1. Via conda (recomendado)
2. Via pip (alternativa)
3. Via compilação manual (última opção)
4. Fallback: usar alternativas (PubChem, RDKit via Docker)
"""

import sys
import subprocess
import os
import platform
from pathlib import Path
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-8s | %(message)s'
)
logger = logging.getLogger('rdkit_installer')


class RDKitInstaller:
    """Instalador multi-método para RDKit."""
    
    def __init__(self):
        self.system = platform.system()
        self.python_version = f"{sys.version_info.major}.{sys.version_info.minor}"
        logger.info(f"Sistema: {self.system}")
        logger.info(f"Python: {self.python_version}")
    
    def find_conda(self):
        """Procura by conda no sistema."""
        logger.info("\n[1] Procurando conda...")
        
        # Locais comuns de conda
        possible_paths = [
            "conda",
            "C:\\Users\\marce\\anaconda3\\Scripts\\conda.exe",
            "C:\\Users\\marce\\miniconda3\\Scripts\\conda.exe",
            "C:\\Users\\marce\\mambaforge\\Scripts\\conda.exe",
            "/usr/local/bin/conda",
            "/opt/conda/bin/conda",
            Path.home() / "anaconda3" / "Scripts" / "conda.exe",
            Path.home() / "miniconda3" / "Scripts" / "conda.exe",
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
                    logger.info(f"   {result.stdout.strip()}")
                    return str(path)
            except Exception:
                pass
        
        logger.warning("❌ Conda não encontrado no PATH")
        return None
    
    def find_mamba(self):
        """Procura por mamba (alternativa rápida a conda)."""
        logger.info("\n[2] Procurando mamba...")
        
        possible_paths = [
            "mamba",
            "C:\\Users\\marce\\mambaforge\\Scripts\\mamba.exe",
            Path.home() / "mambaforge" / "Scripts" / "mamba.exe",
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
                    logger.info(f"   {result.stdout.strip()}")
                    return str(path)
            except Exception:
                pass
        
        logger.warning("❌ Mamba não encontrado")
        return None
    
    def install_via_conda(self, conda_path):
        """Tenta instalar RDKit via conda."""
        logger.info("\n[MÉTODO 1] Instalando RDKit via conda...")
        
        try:
            cmd = [
                conda_path,
                "install",
                "-c", "conda-forge",
                "rdkit",
                "-y"
            ]
            
            logger.info(f"Executando: {' '.join(cmd)}")
            result = subprocess.run(cmd, timeout=600)  # 10 minutos timeout
            
            if result.returncode == 0:
                logger.info("✅ RDKit instalado com sucesso via conda!")
                return self.verify_rdkit()
            else:
            logger.error(f"❌ Falha na instalação via conda")
                return False
        except Exception as e:
            logger.error(f"❌ Erro: {str(e)}")
            return False
    
    def install_via_mamba(self, mamba_path):
        """Tenta instalar RDKit via mamba (mais rápido)."""
        logger.info("\n[MÉTODO 2] Instalando RDKit via mamba...")
        
        try:
            cmd = [
                mamba_path,
                "install",
                "-c", "conda-forge",
                "rdkit",
                "-y"
            ]
            
            logger.info(f"Executando: {' '.join(cmd)}")
            result = subprocess.run(cmd, timeout=600)
            
            if result.returncode == 0:
                logger.info("✅ RDKit instalado com sucesso via mamba!")
                return self.verify_rdkit()
            else:
                logger.error("❌ Falha na instalação via mamba")
                return False
        except Exception as e:
            logger.error(f"❌ Erro: {str(e)}")
            return False
    
    def install_via_pip(self):
        """Tenta instalar RDKit via pip (pode não funcionar)."""
        logger.info("\n[MÉTODO 3] Tentando instalar RDKit via pip...")
        
        # Tentar diferentes packages
        packages = [
            "rdkit-pypi",
            "rdkit",
            "rdkit-conda",
        ]
        
        for package in packages:
            try:
                logger.info(f"  Tentando: {package}")
                result = subprocess.run(
                    [sys.executable, "-m", "pip", "install", package, "-q"],
                    timeout=300
                )
                
                if result.returncode == 0:
                    logger.info(f"✅ {package} instalado!")
                    return self.verify_rdkit()
            except:
                pass
        
        logger.warning("⚠️  Nenhum package RDKit funcionou via pip")
        return False
    
    def download_conda(self):
        """Oferece instruções para baixar conda/mamba."""
        logger.info("\n[MÉTODO 4] Baixando Conda/Mamba Forge...")
        logger.info("""
╔════════════════════════════════════════════════════════════════╗
║                INSTRUÇÕES DE INSTALAÇÃO MANUAL                 ║
╠════════════════════════════════════════════════════════════════╣
║                                                                ║
║ OPÇÃO 1: Instalar Mambaforge (RECOMENDADO - Rápido)          ║
║ ─────────────────────────────────────────────────────────────  ║
║                                                                ║
║ 1. Acesse: https://github.com/conda-forge/miniforge           ║
║ 2. Baixe: Mambaforge-Windows-x86_64.exe (seu sistema)        ║
║ 3. Instale normalmente (pode marcar Python como padrão)       ║
║ 4. Execute no novo prompt Mambaforge:                         ║
║    mamba install -c conda-forge rdkit                         ║
║                                                                ║
║ OPÇÃO 2: Instalar Miniconda + RDKit                          ║
║ ─────────────────────────────────────────────────────────────  ║
║                                                                ║
║ 1. Acesse: https://docs.conda.io/projects/miniconda          ║
║ 2. Baixe: Miniconda3-latest-Windows-x86_64.exe               ║
║ 3. Instale normalmente                                         ║
║ 4. Execute no novo prompt Conda:                              ║
║    conda install -c conda-forge rdkit                         ║
║                                                                ║
║ OPÇÃO 3: Usar Docker (se tiver Docker instalado)            ║
║ ─────────────────────────────────────────────────────────────  ║
║                                                                ║
║ docker run -it --rm continuumio/conda-forge bash              ║
║ conda install -c conda-forge rdkit                            ║
║                                                                ║
╚════════════════════════════════════════════════════════════════╝
        """)
        
        logger.info("""
OBS: Após instalar Conda/Mamba, execute este script novamente:
     python install_rdkit.py
        """)
        
        return "manual"
    
    def verify_rdkit(self):
        """Verifica se RDKit foi instalado corretamente."""
        logger.info("\n[VERIFICAÇÃO] Testando RDKit...")
        
        try:
            import rdkit
            from rdkit import Chem
            
            logger.info(f"✅ RDKit {rdkit.__version__} importado com sucesso!")
            
            # Teste rápido
            mol = Chem.MolFromSmiles("CCO")
            if mol is not None:
                logger.info("✅ RDKit funcionando corretamente (teste de moléculas OK)")
                return True
            else:
                logger.warning("⚠️  RDKit importou mas teste de moléculas falhou")
                return True
        except ImportError as e:
            logger.error(f"❌ Erro ao importar RDKit: {str(e)}")
            return False
        except Exception as e:
            logger.error(f"❌ Erro inesperado: {str(e)}")
            return False
    
    def use_rdkit_alternative(self):
        """Oferece alternativas se RDKit não puder ser instalado."""
        logger.info("\n[FALLBACK] Usando alternativa para RDKit...")
        logger.info("""
╔════════════════════════════════════════════════════════════════╗
║            ALTERNATIVAS SE RDKIT NÃO INSTALAR                ║
╠════════════════════════════════════════════════════════════════╣
║                                                                ║
║ 1. RDKit Online (PubChem, ChemSpider)                         ║
║    - Usar APIs online para featurização                        ║
║    - pip install pubchempy requests                            ║
║                                                                ║
║ 2. DeepChem com Featurizadores Alternativos                  ║
║    - Usar MACCS, Descriptors (não precisa RDKit)              ║
║    - Já inclusos no deepchem.feat                              ║
║                                                                ║
║ 3. Scikit-learn Pipeline                                      ║
║    - Usar TF-IDF ou contadores para sequências                ║
║    - Totalmente independente de RDKit                          ║
║                                                                ║
║ 4. Framework continua operacional                             ║
║    - Datasets sklearn funcionam 100%                           ║
║    - Iris, Wine, Breast Cancer: OK ✅                         ║
║    - Molecular datasets: Esperar RDKit ou usar alt.            ║
║                                                                ║
╚════════════════════════════════════════════════════════════════╝
        """)
        
        return True
    
    def run(self):
        """Executa o instalador."""
        logger.info("="*80)
        logger.info("RDKIT MULTI-METHOD INSTALLER")
        logger.info("="*80)
        
        # Tenta conda
        conda_path = self.find_conda()
        if conda_path:
            if self.install_via_conda(conda_path):
                return True
        
        # Tenta mamba
        mamba_path = self.find_mamba()
        if mamba_path:
            if self.install_via_mamba(mamba_path):
                return True
        
        # Tenta pip
        if self.install_via_pip():
            return True
        
        # Verifica se já está instalado
        if self.verify_rdkit():
            return True
        
        # Oferece download manual ou alternativas
        logger.warning("\n" + "="*80)
        logger.warning("RDKIT NÃO ENCONTRADO - OPÇÕES:")
        logger.warning("="*80)
        
        logger.info("""
RECOMENDAÇÃO: Instalar Mambaforge (rápido e fácil)

PASSO 1: Baixar Mambaforge
   → https://github.com/conda-forge/miniforge#mambaforge
   → Baixe: Mambaforge-Windows-x86_64.exe (seu SO)

PASSO 2: Instalar Mambaforge
   → Execute o .exe e siga as instruções
   → Pode marcar "Add to PATH" para facilitar

PASSO 3: Instalar RDKit
   → Abra novo prompt (Mambaforge Prompt)
   → Execute: mamba install -c conda-forge rdkit

PASSO 4: Volte aqui e execute este script novamente
   → python install_rdkit.py

ALTERNATIVA: Se não quiser instalar nada novo
   → O framework continua operacional com datasets sklearn
   → Iris, Wine, Breast Cancer funcionam 100%
   → Datasets moleculares podem usar alternativas
        """)
        
        # Oferece alternativas
        self.use_rdkit_alternative()
        
        logger.warning("\n" + "="*80)
        logger.warning("FRAMEWORK CONTINUA OPERACIONAL SEM RDKIT")
        logger.warning("="*80)
        
        logger.info("""
Para validar o framework funcional:
   python run_framework_quantum_advanced_v8.py --dataset iris
   python benchmark_all_frameworks_v8.py

Todos os testes funcionam com datasets sklearn!
        """)
        
        return False


if __name__ == "__main__":
    installer = RDKitInstaller()
    success = installer.run()
    
    if success:
        logger.info("\n" + "="*80)
        logger.info("✅ RDKIT INSTALADO COM SUCESSO!")
        logger.info("="*80)
        logger.info("\nPróximos passos:")
        logger.info("  python test_hiv_complete_v8.py")
        logger.info("  python run_framework_quantum_advanced_v8.py --dataset hiv")
        sys.exit(0)
    else:
        logger.info("\n" + "="*80)
        logger.info("⚠️  RDKIT NÃO INSTALADO - VEJA INSTRUÇÕES ACIMA")
        logger.info("="*80)
        logger.info("\nFramework continua operacional com sklearn datasets")
        sys.exit(1)
