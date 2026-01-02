#!/usr/bin/env python3
"""
Conversor de arquivos Markdown da pasta fase4_secoes para DOCX e PDF.

Este script converte automaticamente todos os arquivos .md da pasta
artigo_cientifico/fase4_secoes para os formatos DOCX e PDF usando pandoc.
"""

import os
import sys
import subprocess
import tempfile
import re
from pathlib import Path
from typing import List, Tuple


# Constantes para limpeza de conteúdo
MAX_ERROR_MSG_LENGTH = 200

# Padrões regex para correção de fórmulas matemáticas problemáticas
# Converte notação dimensional $[\mathcal{...}] = $ para formato LaTeX correto
PATTERN_DIMENSIONAL_MATHCAL = r'\$\[\\mathcal\{([^}]+)\}([^\$]*?)\]\s*=\s*\$'
REPLACEMENT_DIMENSIONAL_MATHCAL = r'$\\text{dim}[\\mathcal{\1}\2] =$'

# Converte notação dimensional genérica $[...] = $ para formato LaTeX correto  
PATTERN_DIMENSIONAL_GENERIC = r'\$\[([^\]]+)\]\s*=\s*\$'
REPLACEMENT_DIMENSIONAL_GENERIC = r'$\\text{dim}[\1] =$'


def find_md_files(directory: str) -> List[Path]:
    """
    Encontra todos os arquivos Markdown em um diretório.
    
    Args:
        directory: Caminho do diretório para buscar arquivos .md
        
    Returns:
        Lista de objetos Path com os arquivos .md encontrados
    """
    dir_path = Path(directory)
    if not dir_path.exists():
        print(f"❌ Erro: Diretório não encontrado: {directory}")
        sys.exit(1)
    
    md_files = list(dir_path.glob("*.md"))
    return sorted(md_files)


def convert_to_docx(md_file: Path, output_dir: Path = None) -> bool:
    """
    Converte um arquivo Markdown para DOCX usando pandoc.
    
    Args:
        md_file: Caminho do arquivo Markdown
        output_dir: Diretório de saída (usa o mesmo do arquivo .md se None)
        
    Returns:
        True se a conversão foi bem-sucedida, False caso contrário
    """
    if output_dir is None:
        output_dir = md_file.parent
    
    output_file = output_dir / f"{md_file.stem}.docx"
    
    try:
        # Comando pandoc para converter MD -> DOCX
        cmd = [
            "pandoc",
            str(md_file),
            "-o", str(output_file),
            "--from", "markdown",
            "--to", "docx",
            "--standalone"
        ]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        print(f"✅ DOCX criado: {output_file.name}")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Erro ao converter {md_file.name} para DOCX:")
        print(f"   {e.stderr}")
        return False
    except Exception as e:
        print(f"❌ Erro inesperado ao converter {md_file.name}: {e}")
        return False


def convert_to_pdf(md_file: Path, output_dir: Path = None) -> bool:
    """
    Converte um arquivo Markdown para PDF usando pandoc.
    
    Args:
        md_file: Caminho do arquivo Markdown
        output_dir: Diretório de saída (usa o mesmo do arquivo .md se None)
        
    Returns:
        True se a conversão foi bem-sucedida, False caso contrário
    """
    if output_dir is None:
        output_dir = md_file.parent
    
    output_file = output_dir / f"{md_file.stem}.pdf"
    
    # Ler e sanitizar conteúdo do arquivo
    with open(md_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Substituir caracteres problemáticos para LaTeX
    # Preservar matemática, mas limpar checkboxes e emojis
    content_clean = content.replace('- [x]', '- [DONE]')
    content_clean = content_clean.replace('- [ ]', '- [TODO]')
    # Remover emojis comuns
    content_clean = re.sub(r'[✅✓✔☑✗✘☒❌⚠️]', '', content_clean)
    
    # Fix problematic math patterns with brackets
    # Pattern like $[content] = $ needs to be wrapped properly
    # Replace $[...] pattern when it appears as dimensional analysis
    content_clean = re.sub(PATTERN_DIMENSIONAL_MATHCAL, 
                           REPLACEMENT_DIMENSIONAL_MATHCAL, content_clean)
    content_clean = re.sub(PATTERN_DIMENSIONAL_GENERIC, 
                           REPLACEMENT_DIMENSIONAL_GENERIC, content_clean)
    
    try:
        # Usar arquivo temporário para conversão
        with tempfile.NamedTemporaryFile(mode='w', suffix='.md', delete=False, encoding='utf-8') as tmp:
            tmp.write(content_clean)
            tmp_path = tmp.name
        
        try:
            # Comando pandoc para converter MD -> PDF
            # Usando xelatex com configurações para suportar Unicode e matemática
            cmd = [
                "pandoc",
                tmp_path,
                "-o", str(output_file),
                "--from", "markdown+tex_math_dollars",
                "--to", "pdf",
                "--pdf-engine", "xelatex",
                "--standalone",
                "-V", "geometry:margin=2cm",
                "-V", "mainfont=DejaVu Serif",
                "-V", "sansfont=DejaVu Sans",
                "-V", "monofont=DejaVu Sans Mono"
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            print(f"✅ PDF criado: {output_file.name}")
            return True
            
        finally:
            # Limpar arquivo temporário
            os.unlink(tmp_path)
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Erro ao converter {md_file.name} para PDF:")
        error_msg = e.stderr[:MAX_ERROR_MSG_LENGTH] if e.stderr else "Erro desconhecido"
        print(f"   {error_msg}")
        return False
    except Exception as e:
        print(f"❌ Erro inesperado ao converter {md_file.name}: {e}")
        return False


def convert_all_files(directory: str) -> Tuple[int, int, int]:
    """
    Converte todos os arquivos .md de um diretório para DOCX e PDF.
    
    Args:
        directory: Caminho do diretório com arquivos .md
        
    Returns:
        Tupla com (total de arquivos, conversões DOCX bem-sucedidas, conversões PDF bem-sucedidas)
    """
    md_files = find_md_files(directory)
    
    if not md_files:
        print(f"⚠️ Nenhum arquivo .md encontrado em {directory}")
        return 0, 0, 0
    
    print(f"\n📁 Encontrados {len(md_files)} arquivos Markdown em {directory}")
    print("=" * 70)
    
    docx_success = 0
    pdf_success = 0
    
    for i, md_file in enumerate(md_files, 1):
        print(f"\n[{i}/{len(md_files)}] Convertendo: {md_file.name}")
        print("-" * 70)
        
        # Converter para DOCX
        if convert_to_docx(md_file):
            docx_success += 1
        
        # Converter para PDF
        if convert_to_pdf(md_file):
            pdf_success += 1
    
    return len(md_files), docx_success, pdf_success


def main():
    """Função principal do script."""
    
    # Definir o diretório padrão
    base_dir = Path(__file__).resolve().parent.parent
    default_dir = base_dir / "artigo_cientifico" / "fase4_secoes"
    
    # Permitir especificar diretório via argumento
    if len(sys.argv) > 1:
        target_dir = sys.argv[1]
    else:
        target_dir = str(default_dir)
    
    print("\n" + "=" * 70)
    print("🔄 CONVERSOR DE MARKDOWN PARA DOCX E PDF")
    print("=" * 70)
    print(f"📂 Diretório: {target_dir}")
    
    # Verificar se pandoc está instalado
    try:
        result = subprocess.run(
            ["pandoc", "--version"],
            capture_output=True,
            text=True,
            check=True
        )
        pandoc_version = result.stdout.split('\n')[0]
        print(f"✅ {pandoc_version}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("❌ Erro: Pandoc não está instalado!")
        print("   Instale com: sudo apt-get install pandoc")
        sys.exit(1)
    
    # Converter todos os arquivos
    total, docx_ok, pdf_ok = convert_all_files(target_dir)
    
    # Resumo final
    print("\n" + "=" * 70)
    print("📊 RESUMO DA CONVERSÃO")
    print("=" * 70)
    print(f"📄 Total de arquivos: {total}")
    print(f"✅ DOCX criados: {docx_ok}/{total}")
    print(f"✅ PDF criados: {pdf_ok}/{total}")
    
    if docx_ok == total and pdf_ok == total:
        print("\n🎉 Todas as conversões foram concluídas com sucesso!")
        return 0
    else:
        print("\n⚠️ Algumas conversões falharam. Verifique os erros acima.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
