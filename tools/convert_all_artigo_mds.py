#!/usr/bin/env python3
"""
Conversor de TODOS os arquivos Markdown da pasta artigo_cientifico para DOCX e PDF.

Este script converte automaticamente todos os arquivos .md encontrados em
artigo_cientifico e suas subpastas para os formatos DOCX e PDF usando pandoc.
"""

import os
import sys
import subprocess
import tempfile
import re
from pathlib import Path
from typing import List, Tuple, Dict


# Constantes para limpeza de conteúdo
MAX_ERROR_MSG_LENGTH = 200

# Padrões regex para correção de fórmulas matemáticas problemáticas
# Converte notação dimensional $[\mathcal{...}] = $ para formato LaTeX correto
PATTERN_DIMENSIONAL_MATHCAL = r'\$\[\\mathcal\{([^}]+)\}([^\$]*?)\]\s*=\s*\$'
REPLACEMENT_DIMENSIONAL_MATHCAL = r'$\\text{dim}[\\mathcal{\1}\2] =$'

# Converte notação dimensional genérica $[...] = $ para formato LaTeX correto  
PATTERN_DIMENSIONAL_GENERIC = r'\$\[([^\]]+)\]\s*=\s*\$'
REPLACEMENT_DIMENSIONAL_GENERIC = r'$\\text{dim}[\1] =$'


def find_md_files_recursive(directory: str) -> List[Path]:
    """
    Encontra todos os arquivos Markdown em um diretório e suas subpastas.
    
    Args:
        directory: Caminho do diretório raiz para buscar arquivos .md
        
    Returns:
        Lista de objetos Path com os arquivos .md encontrados
    """
    dir_path = Path(directory)
    if not dir_path.exists():
        print(f"❌ Erro: Diretório não encontrado: {directory}")
        sys.exit(1)
    
    # Usar rglob para busca recursiva
    md_files = list(dir_path.rglob("*.md"))
    return sorted(md_files)


def convert_to_docx(md_file: Path) -> bool:
    """
    Converte um arquivo Markdown para DOCX usando pandoc.
    
    Args:
        md_file: Caminho do arquivo Markdown
        
    Returns:
        True se a conversão foi bem-sucedida, False caso contrário
    """
    output_file = md_file.parent / f"{md_file.stem}.docx"
    
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
        
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"      ❌ Erro DOCX: {e.stderr[:MAX_ERROR_MSG_LENGTH]}")
        return False
    except Exception as e:
        print(f"      ❌ Erro inesperado DOCX: {str(e)[:MAX_ERROR_MSG_LENGTH]}")
        return False


def convert_to_pdf(md_file: Path) -> bool:
    """
    Converte um arquivo Markdown para PDF usando pandoc.
    
    Args:
        md_file: Caminho do arquivo Markdown
        
    Returns:
        True se a conversão foi bem-sucedida, False caso contrário
    """
    output_file = md_file.parent / f"{md_file.stem}.pdf"
    
    # Ler e sanitizar conteúdo do arquivo
    with open(md_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Substituir caracteres problemáticos para LaTeX
    # Preservar matemática, mas limpar checkboxes e emojis
    content_clean = content.replace('- [x]', '- [DONE]')
    content_clean = content_clean.replace('- [ ]', '- [TODO]')
    # Remover emojis comuns
    content_clean = re.sub(r'[✅✓✔☑✗✘☒❌⚠️📊📝📚🎯🔬🔗📎ℹ️📄🔄]', '', content_clean)
    
    # Fix problematic math patterns with brackets
    content_clean = re.sub(PATTERN_DIMENSIONAL_MATHCAL, 
                           REPLACEMENT_DIMENSIONAL_MATHCAL, content_clean)
    content_clean = re.sub(PATTERN_DIMENSIONAL_GENERIC, 
                           REPLACEMENT_DIMENSIONAL_GENERIC, content_clean)
    
    # Fix double backslashes in inline math (in markdown tables)
    # Replace \\| with \| (single backslash) in inline math
    content_clean = re.sub(r'\$([^\$]*?)\\\\([|])', r'$\1\\\2', content_clean)
    
    # Fix \\| patterns in general (norm notation)
    content_clean = re.sub(r'\\\\\\\\', r'\\\\', content_clean)  # Reduce quadruple to double
    
    # Ensure \text commands are properly closed (common LaTeX error)
    content_clean = re.sub(r'\\text\{([^}]*?)$', r'\\text{\1}', content_clean, flags=re.MULTILINE)
    
    try:
        # Usar arquivo temporário para conversão
        with tempfile.NamedTemporaryFile(mode='w', suffix='.md', delete=False, encoding='utf-8') as tmp:
            tmp.write(content_clean)
            tmp_path = tmp.name
        
        try:
            # Comando pandoc para converter MD -> PDF
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
            
            return True
            
        finally:
            # Limpar arquivo temporário
            os.unlink(tmp_path)
        
    except subprocess.CalledProcessError as e:
        print(f"      ❌ Erro PDF: {e.stderr[:MAX_ERROR_MSG_LENGTH] if e.stderr else 'Erro desconhecido'}")
        return False
    except Exception as e:
        print(f"      ❌ Erro inesperado PDF: {str(e)[:MAX_ERROR_MSG_LENGTH]}")
        return False


def convert_all_files(directory: str) -> Tuple[Dict[str, int], List[Path]]:
    """
    Converte todos os arquivos .md de um diretório (recursivamente) para DOCX e PDF.
    
    Args:
        directory: Caminho do diretório raiz com arquivos .md
        
    Returns:
        Tupla com (dict de estatísticas, lista de arquivos com falha)
    """
    md_files = find_md_files_recursive(directory)
    
    if not md_files:
        print(f"⚠️ Nenhum arquivo .md encontrado em {directory}")
        return {"total": 0, "docx_success": 0, "pdf_success": 0}, []
    
    print(f"\n📁 Encontrados {len(md_files)} arquivos Markdown em {directory} e subpastas")
    print("=" * 80)
    
    stats = {
        "total": len(md_files),
        "docx_success": 0,
        "pdf_success": 0,
        "both_success": 0
    }
    
    # Agrupar arquivos por diretório para melhor visualização
    files_by_dir = {}
    for md_file in md_files:
        rel_path = md_file.relative_to(directory)
        dir_name = str(rel_path.parent) if rel_path.parent != Path('.') else "raiz"
        if dir_name not in files_by_dir:
            files_by_dir[dir_name] = []
        files_by_dir[dir_name].append(md_file)
    
    failed_files = []
    
    for dir_name in sorted(files_by_dir.keys()):
        files = files_by_dir[dir_name]
        print(f"\n📂 {dir_name} ({len(files)} arquivo(s))")
        print("-" * 80)
        
        for md_file in files:
            print(f"   📄 {md_file.name}")
            
            # Converter para DOCX
            docx_ok = convert_to_docx(md_file)
            if docx_ok:
                stats["docx_success"] += 1
                print(f"      ✅ DOCX: {md_file.stem}.docx")
            else:
                failed_files.append(md_file)
            
            # Converter para PDF
            pdf_ok = convert_to_pdf(md_file)
            if pdf_ok:
                stats["pdf_success"] += 1
                print(f"      ✅ PDF: {md_file.stem}.pdf")
            elif not docx_ok:  # Só adicionar uma vez à lista de falhas
                pass
            else:
                failed_files.append(md_file)
            
            # Contar sucessos completos
            if docx_ok and pdf_ok:
                stats["both_success"] += 1
    
    return stats, failed_files


def main():
    """Função principal do script."""
    
    # Definir o diretório padrão
    base_dir = Path(__file__).resolve().parent.parent
    default_dir = base_dir / "artigo_cientifico"
    
    # Permitir especificar diretório via argumento
    if len(sys.argv) > 1:
        target_dir = sys.argv[1]
    else:
        target_dir = str(default_dir)
    
    print("\n" + "=" * 80)
    print("🔄 CONVERSOR COMPLETO DE MARKDOWN PARA DOCX E PDF")
    print("=" * 80)
    print(f"📂 Diretório: {target_dir}")
    print("🔍 Modo: Busca recursiva em todas as subpastas")
    
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
        print("   Instale com: sudo apt-get install pandoc texlive-xetex texlive-fonts-recommended")
        sys.exit(1)
    
    # Converter todos os arquivos
    stats, failed_files = convert_all_files(target_dir)
    
    # Resumo final
    print("\n" + "=" * 80)
    print("📊 RESUMO DA CONVERSÃO")
    print("=" * 80)
    print(f"📄 Total de arquivos: {stats['total']}")
    print(f"✅ DOCX criados: {stats['docx_success']}/{stats['total']}")
    print(f"✅ PDF criados: {stats['pdf_success']}/{stats['total']}")
    print(f"🎯 Ambos criados: {stats['both_success']}/{stats['total']}")
    
    if failed_files:
        print(f"\n⚠️ Arquivos com falha ({len(failed_files)}):")
        for f in failed_files[:10]:  # Mostrar apenas os primeiros 10
            print(f"   - {f.relative_to(target_dir)}")
        if len(failed_files) > 10:
            print(f"   ... e mais {len(failed_files) - 10} arquivo(s)")
    
    if stats['both_success'] == stats['total']:
        print("\n🎉 Todas as conversões foram concluídas com sucesso!")
        return 0
    else:
        print("\n⚠️ Algumas conversões falharam. Verifique os erros acima.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
