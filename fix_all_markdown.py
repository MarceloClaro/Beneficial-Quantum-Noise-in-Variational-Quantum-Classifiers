#!/usr/bin/env python3
"""
Script para corrigir automaticamente erros de markdownlint em todos os arquivos .md do projeto
"""
import re
import sys
from pathlib import Path
from typing import List

def fix_markdown_content(content: str) -> str:
    """
    Aplica correções de markdown no conteúdo.

    Args:
        content: Conteúdo do arquivo markdown

    Returns:
        Conteúdo corrigido
    """
    # 1. MD009: Remover trailing spaces (exceto quebras de linha duplas)
    lines = content.split('\n')
    fixed_lines = []
    for line in lines:
        # Remove espaços no final, exceto se forem 2 ou mais (quebra de linha MD)
        if line.endswith('  '):
            fixed_lines.append(line)  # Mantém quebras de linha MD
        else:
            fixed_lines.append(line.rstrip())
    content = '\n'.join(fixed_lines)

    # 2. MD040: Adicionar linguagem aos blocos de código sem especificação
    content = re.sub(r'```\n((?:(?!```).)*)```', r'```text\n\1```', content, flags=re.DOTALL)

    # 3. MD036: Converter ênfase em headings (linhas que são só **texto**)
    lines = content.split('\n')
    fixed_lines = []
    for i, line in enumerate(lines):
        if re.match(r'^\*\*[^*]+\*\*$', line.strip()) and len(line.strip()) < 100:
            if i + 1 < len(lines) and lines[i + 1].strip().startswith('-'):
                text = re.sub(r'^\*\*([^*]+)\*\*$', r'\1', line.strip())
                fixed_lines.append(f'#### {text}')
                continue
        fixed_lines.append(line)
    content = '\n'.join(fixed_lines)

    # 4. MD032: Adicionar linhas em branco antes/depois de listas
    lines = content.split('\n')
    fixed_lines = []
    in_list = False

    for i, line in enumerate(lines):
        is_list_item = line.strip().startswith(('-', '*', '+')) or re.match(r'^\d+\.', line.strip())
        prev_line = lines[i-1] if i > 0 else ''

        if is_list_item and not in_list:
            if prev_line.strip() and not prev_line.strip().startswith('#'):
                fixed_lines.append('')
            in_list = True
        elif not is_list_item and in_list:
            if i > 0 and lines[i-1].strip():
                fixed_lines.insert(len(fixed_lines), '')
            in_list = False

        fixed_lines.append(line)

    content = '\n'.join(fixed_lines)

    # 5. MD022: Adicionar linhas em branco antes/depois de headings
    lines = content.split('\n')
    fixed_lines = []

    for i, line in enumerate(lines):
        is_heading = line.strip().startswith('#')
        prev_line = lines[i-1] if i > 0 else ''
        next_line = lines[i+1] if i < len(lines) - 1 else ''

        if is_heading and i > 0 and prev_line.strip() and not prev_line.strip().startswith('#'):
            fixed_lines.append('')

        fixed_lines.append(line)

        if is_heading and next_line and next_line.strip() and not next_line.strip().startswith('#'):
            if i < len(lines) - 1 and not lines[i+1].strip():
                continue

    content = '\n'.join(fixed_lines)

    # 6. MD031: Linhas em branco ao redor de blocos de código
    content = re.sub(r'([^\n])\n```', r'\1\n\n```', content)
    content = re.sub(r'```\n([^\n])', r'```\n\n\1', content)

    # 7. MD034: Converter URLs nuas em links
    def wrap_bare_url(match):
        url = match.group(0)
        if url.startswith('<') and url.endswith('>'):
            return url
        return f'<{url}>'

    content = re.sub(r'(?<![\[(<])(https?://[^\s\)<>]+)(?![\])>])', wrap_bare_url, content)

    return content


def process_markdown_file(file_path: Path) -> bool:
    """
    Processa um arquivo markdown aplicando correções.

    Args:
        file_path: Caminho para o arquivo markdown

    Returns:
        True se o arquivo foi modificado, False caso contrário
    """
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            original_content = f.read()

        fixed_content = fix_markdown_content(original_content)

        if fixed_content != original_content:
            with open(file_path, 'w', encoding='utf-8', newline='\n') as f:
                f.write(fixed_content)
            return True
        return False
    except Exception as e:
        print(f"❌ Erro ao processar {file_path}: {e}")
        return False


def find_all_markdown_files(root_dir: Path) -> List[Path]:
    """
    Encontra todos os arquivos .md no diretório e subdiretórios.

    Args:
        root_dir: Diretório raiz para busca

    Returns:
        Lista de caminhos para arquivos .md
    """
    return list(root_dir.rglob('*.md'))


def main():
    """Função principal."""
    # Determina o diretório raiz
    if len(sys.argv) > 1:
        root_dir = Path(sys.argv[1])
    else:
        root_dir = Path(__file__).parent

    print(f"\n🔍 Buscando arquivos .md em: {root_dir}")

    # Encontra todos os arquivos markdown
    md_files = find_all_markdown_files(root_dir)
    total_files = len(md_files)

    print(f"📄 Encontrados {total_files} arquivos markdown\n")

    # Processa cada arquivo
    modified_count = 0
    unchanged_count = 0

    for i, md_file in enumerate(md_files, 1):
        relative_path = md_file.relative_to(root_dir)
        print(f"[{i}/{total_files}] Processando: {relative_path}...", end=' ')

        if process_markdown_file(md_file):
            print("✅ Corrigido")
            modified_count += 1
        else:
            print("⏭️ Sem mudanças")
            unchanged_count += 1

    # Resumo final
    print(f"\n{'='*60}")
    print("📊 RESUMO DA CORREÇÃO")
    print(f"{'='*60}")
    print(f"Total de arquivos processados: {total_files}")
    print(f"✅ Arquivos corrigidos: {modified_count}")
    print(f"⏭️ Arquivos sem mudanças: {unchanged_count}")
    print(f"{'='*60}\n")


if __name__ == '__main__':
    main()
