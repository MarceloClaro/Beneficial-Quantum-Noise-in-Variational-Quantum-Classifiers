#!/usr/bin/env python3
"""
Conversor de Markdown para PDF com suporte completo a LaTeX via MathJax.

Este script converte Markdown para HTML com renderização MathJax e abre
no navegador para impressão em PDF (Ctrl+P → Save as PDF).

Esta é a solução mais confiável para documentos científicos com equações matemáticas.
"""

import sys
import os
from pathlib import Path
from typing import Optional

def convert_md_to_html_with_mathjax(md_file: str, output_html: Optional[str] = None) -> str:
    """
    Converte Markdown para HTML com suporte MathJax para LaTeX.

    Args:
        md_file: Caminho para o arquivo Markdown
        output_html: Caminho opcional para o HTML de saída

    Returns:
        Caminho do arquivo HTML gerado
    """

    # Verificar se o arquivo existe
    if not os.path.exists(md_file):
        print(f"❌ Erro: Arquivo não encontrado: {md_file}")
        sys.exit(1)

    # Determinar nome do HTML de saída
    if output_html is None:
        output_html = Path(md_file).stem + '.html'

    # Ler o conteúdo do Markdown
    with open(md_file, 'r', encoding='utf-8') as f:
        md_content = f.read()

    # Converter Markdown para HTML
    try:
        import markdown2
        # Usar extras para melhor formatação
        html_body = markdown2.markdown(
            md_content,
            extras=[
                'tables',
                'fenced-code-blocks',
                'strike',
                'cuddled-lists',
                'header-ids',
                'footnotes'
            ]
        )
    except ImportError:
        print("❌ Erro: markdown2 não está instalado")
        print("   Execute: pip install markdown2")
        sys.exit(1)

    # Template HTML completo com MathJax
    html_template = """<!DOCTYPE html>
<html lang="pt-BR">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>

    <!-- MathJax para renderização de LaTeX -->
    <script>
        MathJax = {{
            tex: {{
                inlineMath: [['$', '$'], ['\\(', '\\)']],
                displayMath: [['$$', '$$'], ['\\[', '\\]']],
                processEscapes: true,
                processEnvironments: true
            }},
            options: {{
                skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre']
            }},
            svg: {{
                fontCache: 'global'
            }}
        }};
    </script>
    <script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js" id="MathJax-script" async></script>

    <!-- Estilos para impressão -->
    <style>
        /* === Estilos para Tela === */
        @media screen {{
            body {{
                font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
                line-height: 1.6;
                color: #2c3e50;
                max-width: 900px;
                margin: 0 auto;
                padding: 40px 20px;
                background-color: #f5f7fa;
            }}

            .container {{
                background: white;
                padding: 60px;
                border-radius: 8px;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            }}
        }}

        /* === Estilos para Impressão === */
        @media print {{
            @page {{
                size: A4;
                margin: 2cm;
            }}

            body {{
                font-family: 'Times New Roman', Times, serif;
                line-height: 1.5;
                color: #000;
                background: white;
            }}

            .container {{
                padding: 0;
                box-shadow: none;
            }}

            /* Evitar quebras de página indesejadas */
            h1, h2, h3, h4, h5, h6 {{
                page-break-after: avoid;
                page-break-inside: avoid;
            }}

            table, figure, pre {{
                page-break-inside: avoid;
            }}

            /* Forçar quebra antes de seções principais */
            h2 {{
                page-break-before: always;
            }}

            h2:first-of-type {{
                page-break-before: avoid;
            }}
        }}

        /* === Tipografia === */
        h1 {{
            font-size: 2.5em;
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 15px;
            margin-top: 0;
            margin-bottom: 30px;
        }}

        h2 {{
            font-size: 2em;
            color: #34495e;
            border-bottom: 2px solid #95a5a6;
            padding-bottom: 10px;
            margin-top: 40px;
            margin-bottom: 20px;
        }}

        h3 {{
            font-size: 1.5em;
            color: #7f8c8d;
            margin-top: 30px;
            margin-bottom: 15px;
        }}

        h4 {{
            font-size: 1.2em;
            color: #95a5a6;
            margin-top: 20px;
            margin-bottom: 10px;
        }}

        p {{
            margin: 15px 0;
            text-align: justify;
        }}

        strong, b {{
            color: #2c3e50;
            font-weight: 600;
        }}

        em, i {{
            color: #34495e;
        }}

        /* === Listas === */
        ul, ol {{
            margin: 15px 0;
            padding-left: 30px;
        }}

        li {{
            margin: 8px 0;
        }}

        /* === Tabelas === */
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 25px 0;
            font-size: 0.95em;
        }}

        table thead {{
            background-color: #3498db;
            color: white;
        }}

        table th, table td {{
            padding: 12px 15px;
            border: 1px solid #ddd;
            text-align: left;
        }}

        table tbody tr:nth-child(even) {{
            background-color: #f8f9fa;
        }}

        table tbody tr:hover {{
            background-color: #e8f4f8;
        }}

        @media print {{
            table thead {{
                background-color: #2980b9 !important;
                -webkit-print-color-adjust: exact;
                print-color-adjust: exact;
            }}

            table tbody tr:nth-child(even) {{
                background-color: #f0f0f0 !important;
                -webkit-print-color-adjust: exact;
                print-color-adjust: exact;
            }}
        }}

        /* === Código === */
        code {{
            background-color: #f4f4f4;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: 'Courier New', monospace;
            font-size: 0.9em;
            color: #c7254e;
        }}

        pre {{
            background-color: #f8f9fa;
            border-left: 4px solid #3498db;
            padding: 15px;
            overflow-x: auto;
            border-radius: 4px;
            margin: 20px 0;
        }}

        pre code {{
            background: none;
            color: #2c3e50;
            padding: 0;
        }}

        /* === Equações Matemáticas === */
        .MathJax {{
            font-size: 1.1em !important;
        }}

        mjx-container[display="true"] {{
            margin: 20px 0 !important;
        }}

        /* === Separadores === */
        hr {{
            border: none;
            border-top: 2px solid #ecf0f1;
            margin: 40px 0;
        }}

        /* === Links === */
        a {{
            color: #3498db;
            text-decoration: none;
        }}

        a:hover {{
            text-decoration: underline;
        }}

        @media print {{
            a {{
                color: #000;
                text-decoration: none;
            }}

            a[href]:after {{
                content: " (" attr(href) ")";
                font-size: 0.8em;
                color: #666;
            }}
        }}

        /* === Blocos de destaque === */
        blockquote {{
            border-left: 4px solid #3498db;
            padding-left: 20px;
            margin: 20px 0;
            color: #555;
            font-style: italic;
            background-color: #f8f9fa;
            padding: 15px 20px;
        }}

        /* === Checkmarks e status === */
        .status-icon {{
            display: inline-block;
            width: 1.2em;
            text-align: center;
        }}
    </style>
</head>
<body>
    <div class="container">
        {body}
    </div>

    <!-- Instruções de impressão -->
    <script>
        window.onload = function() {{
            console.log('✅ Documento carregado com sucesso!');
            console.log('📄 Para salvar como PDF:');
            console.log('   1. Pressione Ctrl+P (Windows) ou Cmd+P (Mac)');
            console.log('   2. Selecione "Salvar como PDF" ou "Microsoft Print to PDF"');
            console.log('   3. Ajuste margens e layout conforme necessário');
            console.log('   4. Clique em "Salvar"');
        }};
    </script>
</body>
</html>
"""

    # Gerar HTML final
    title = Path(md_file).stem.replace('_', ' ').title()
    html_final = html_template.format(title=title, body=html_body)

    # Salvar HTML
    with open(output_html, 'w', encoding='utf-8') as f:
        f.write(html_final)

    print(f"✅ HTML gerado com sucesso: {output_html}")
    print("📐 MathJax habilitado para renderização de LaTeX")
    print("\n📄 Para gerar PDF:")
    print("   1. O arquivo será aberto no navegador")
    print("   2. Aguarde a renderização das equações (2-3 segundos)")
    print("   3. Pressione Ctrl+P (Windows) ou Cmd+P (Mac)")
    print("   4. Selecione 'Salvar como PDF'")
    print("   5. Clique em 'Salvar'")

    return output_html


def main():
    """Função principal do script."""

    if len(sys.argv) < 2:
        print("Uso: python md_to_pdf_mathjax.py <arquivo.md> [saida.html]")
        print("\nExemplo:")
        print("  python md_to_pdf_mathjax.py OBJETIVOS_PROJETO.md")
        sys.exit(1)

    md_file = sys.argv[1]
    output_html = sys.argv[2] if len(sys.argv) > 2 else None

    # Converter para HTML
    html_file = convert_md_to_html_with_mathjax(md_file, output_html)

    # Abrir no navegador
    try:
        import webbrowser
        abs_path = os.path.abspath(html_file)
        webbrowser.open(f'file://{abs_path}')
        print(f"\n🌐 Abrindo no navegador: {abs_path}")
    except Exception as e:
        print(f"\n⚠️ Não foi possível abrir automaticamente: {e}")
        print(f"   Abra manualmente: {os.path.abspath(html_file)}")


if __name__ == '__main__':
    main()
