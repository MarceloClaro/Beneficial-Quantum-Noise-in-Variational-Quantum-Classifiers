#!/usr/bin/env python3
"""
Script para corrigir automaticamente erros comuns de markdownlint no README.md
"""
import re
from pathlib import Path

def fix_markdown_lint(file_path: str) -> None:
    """
    Corrige erros comuns de markdownlint.

    Args:
        file_path: Caminho para o arquivo README.md
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    original_content = content

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
    # Detecta blocos que começam com ``` seguido de newline
    content = re.sub(r'```\n((?:(?!```).)*)```', r'```text\n\1```', content, flags=re.DOTALL)

    # 3. MD036: Converter ênfase em headings (linhas que são só **texto**)
    lines = content.split('\n')
    fixed_lines = []
    for i, line in enumerate(lines):
        # Se linha é apenas **texto** ou *texto* e não tem nada antes/depois
        if re.match(r'^\*\*[^*]+\*\*$', line.strip()) and len(line.strip()) < 100:
            # Verifica contexto - se parece ser um heading (lista abaixo)
            if i + 1 < len(lines) and lines[i + 1].strip().startswith('-'):
                # Converte para heading nível 4
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
        next_line = lines[i+1] if i < len(lines) - 1 else ''

        # Se começamos uma lista e linha anterior não é vazia e não é heading
        if is_list_item and not in_list:
            if prev_line.strip() and not prev_line.strip().startswith('#'):
                fixed_lines.append('')
            in_list = True
        # Se terminamos uma lista
        elif not is_list_item and in_list:
            if i > 0 and lines[i-1].strip():  # Linha anterior era lista
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

        # Adiciona linha antes do heading se necessário
        if is_heading and i > 0 and prev_line.strip() and not prev_line.strip().startswith('#'):
            fixed_lines.append('')

        fixed_lines.append(line)

        # Adiciona linha depois do heading se necessário
        if is_heading and next_line and next_line.strip() and not next_line.strip().startswith('#'):
            # Não adiciona se próxima linha já é vazia ou é outro heading
            if i < len(lines) - 1 and not lines[i+1].strip():
                continue

    content = '\n'.join(fixed_lines)

    # 6. MD031: Linhas em branco ao redor de blocos de código
    content = re.sub(r'([^\n])\n```', r'\1\n\n```', content)
    content = re.sub(r'```\n([^\n])', r'```\n\n\1', content)

    # 7. MD034: Converter URLs nuas em links
    # Detecta URLs que não estão em links ou código
    def wrap_bare_url(match):
        url = match.group(0)
        if url.startswith('<') and url.endswith('>'):
            return url
        return f'<{url}>'

    # URLs em texto normal (não em links ou código)
    content = re.sub(r'(?<![\[(<])(https?://[^\s\)<>]+)(?![\])>])', wrap_bare_url, content)

    # Salvar se houver mudanças
    if content != original_content:
        with open(file_path, 'w', encoding='utf-8', newline='\n') as f:
            f.write(content)
        print(f"✅ Arquivo corrigido: {file_path}")
        print("📊 Total de mudanças aplicadas")
    else:
        print("ℹ️ Nenhuma mudança necessária")


if __name__ == '__main__':
    readme_path = Path(__file__).parent / 'README.md'
    fix_markdown_lint(str(readme_path))
