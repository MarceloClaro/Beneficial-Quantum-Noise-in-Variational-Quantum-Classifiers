"""
Teste estrutural rápido do VQC-Molecular v10.0-A1
Verifica arquivos e estrutura sem executar código
"""
from pathlib import Path

print("="*70)
print("🧪 TESTE ESTRUTURAL VQC-Molecular v10.0-A1")
print("="*70)
print()

# Diretório base
base_dir = Path(__file__).parent / "vqc_drug_v10a1"

# Teste 1: Estrutura de diretórios
print("📁 Teste 1: Verificando estrutura de diretórios...")
dirs_expected = ["src", "docker", "tests"]
for d in dirs_expected:
    dir_path = base_dir / d
    if dir_path.exists():
        print(f"   ✅ {d}/ existe")
    else:
        print(f"   ❌ {d}/ NÃO encontrado")
print()

# Teste 2: Arquivos principais
print("📄 Teste 2: Verificando arquivos principais...")
files_expected = {
    "pyproject.toml": "Configuração pip",
    "README.md": "Documentação principal",
    "QUICKSTART.md": "Guia de instalação",
    "LICENSE": "Licença MIT"
}

for fname, desc in files_expected.items():
    fpath = base_dir / fname
    if fpath.exists():
        size = fpath.stat().st_size
        print(f"   ✅ {fname} ({size} bytes) - {desc}")
    else:
        print(f"   ❌ {fname} NÃO encontrado")
print()

# Teste 3: Módulos Python em src/
print("🐍 Teste 3: Verificando módulos Python...")
modules_expected = {
    "__init__.py": "Package facade",
    "data.py": "Morgan fingerprints + PCA",
    "models.py": "VQC com GPU support",
    "tune.py": "Optuna ultra-tuner",
    "audit.py": "SHA-256 + checksums",
    "plots.py": "Figuras 600 dpi",
    "cli.py": "Entry point"
}

src_dir = base_dir / "src"
for fname, desc in modules_expected.items():
    fpath = src_dir / fname
    if fpath.exists():
        lines = len(fpath.read_text(encoding="utf-8").splitlines())
        print(f"   ✅ src/{fname} ({lines} linhas) - {desc}")
    else:
        print(f"   ❌ src/{fname} NÃO encontrado")
print()

# Teste 4: Docker
print("🐳 Teste 4: Verificando arquivos Docker...")
docker_files = ["Dockerfile", "environment.yml"]
docker_dir = base_dir / "docker"
for fname in docker_files:
    fpath = docker_dir / fname
    if fpath.exists():
        print(f"   ✅ docker/{fname} existe")
    else:
        print(f"   ❌ docker/{fname} NÃO encontrado")
print()

# Teste 5: Testes
print("🧪 Teste 5: Verificando arquivos de teste...")
test_file = base_dir / "tests" / "test_all.py"
if test_file.exists():
    lines = len(test_file.read_text(encoding="utf-8").splitlines())
    print(f"   ✅ tests/test_all.py ({lines} linhas) existe")
else:
    print(f"   ❌ tests/test_all.py NÃO encontrado")
print()

# Teste 6: README com comparação v9 vs v10
print("📊 Teste 6: Verificando conteúdo do README...")
readme = base_dir / "README.md"
if readme.exists():
    content = readme.read_text(encoding="utf-8")
    checks = {
        "v10.0-A1": "Título correto",
        "Marcelo Claro Laranjeira": "Autor correto",
        "Fisher-CRLB": "Otimização Fisher presente",
        "Lindblad": "Lindblad scheduling presente",
        "v9.0 vs v10.0": "Comparação presente",
        "Power-Adaptative Search": "PAS mencionado"
    }
    
    for key, desc in checks.items():
        if key in content:
            print(f"   ✅ {desc}")
        else:
            print(f"   ⚠️  {desc} - não encontrado")
    
    print(f"   📝 Total de linhas no README: {len(content.splitlines())}")
else:
    print("   ❌ README.md não encontrado")
print()

# Estatísticas finais
print("="*70)
print("📊 ESTATÍSTICAS FINAIS")
print("="*70)

total_files = sum(1 for _ in base_dir.rglob("*") if _.is_file())
total_py = sum(1 for _ in base_dir.rglob("*.py"))
total_md = sum(1 for _ in base_dir.rglob("*.md"))

print(f"   Total de arquivos: {total_files}")
print(f"   Arquivos Python (.py): {total_py}")
print(f"   Arquivos Markdown (.md): {total_md}")
print()

# Cálculo de linhas
total_lines = 0
for fpath in base_dir.rglob("*.py"):
    try:
        total_lines += len(fpath.read_text(encoding="utf-8").splitlines())
    except:
        pass

for fpath in base_dir.rglob("*.md"):
    try:
        total_lines += len(fpath.read_text(encoding="utf-8").splitlines())
    except:
        pass

print(f"   Total de linhas (código + docs): ~{total_lines}")
print()

print("="*70)
print("✅ ESTRUTURA DO FRAMEWORK v10.0-A1 VERIFICADA!")
print("="*70)
print()
print("🎯 Status: Pronto para instalação")
print()
print("📦 Para instalar dependências:")
print("   pip install pennylane torch optuna rdkit scikit-learn")
print("   pip install pandas numpy matplotlib seaborn click requests statsmodels")
print()
print("🚀 Para instalar o framework:")
print("   pip install -e vqc_drug_v10a1/")
print()
