#!/usr/bin/env python3
"""
Converte Markdown para PDF usando markdown2 e weasyprint.
"""
import sys
import os
from pathlib import Path


def convert_md_to_pdf(md_file: str, output_pdf: str = None):
    """
    Converte um arquivo Markdown para PDF usando markdown2 e pdfkit.
    
    Args:
        md_file: Caminho para o arquivo Markdown de entrada
        output_pdf: Caminho opcional para o arquivo PDF de saída.
                   Se None, usa o mesmo nome do arquivo de entrada com extensão .pdf
    
    Returns:
        None
        
    Raises:
        SystemExit: Se o arquivo não for encontrado ou markdown2 não estiver instalado
        
    Example:
        >>> convert_md_to_pdf('README.md', 'README.pdf')
        >>> convert_md_to_pdf('ARTICLE.md')  # Gera ARTICLE.pdf
    """

    # Verificar se o arquivo existe
    if not os.path.exists(md_file):
        print(f"❌ Erro: Arquivo não encontrado: {md_file}")
        sys.exit(1)

    # Determinar nome do PDF de saída
    if output_pdf is None:
        output_pdf = Path(md_file).stem + '.pdf'

    # Ler o conteúdo do Markdown
    with open(md_file, 'r', encoding='utf-8') as f:
        md_content = f.read()

    # Converter Markdown para HTML
    try:
        import markdown2
        html_content = markdown2.markdown(md_content, extras=['tables', 'fenced-code-blocks'])
    except ImportError:
        print("❌ Erro: markdown2 não está instalado. Execute: pip install markdown2")
        sys.exit(1)

    # Adicionar CSS para melhor formatação
    html_with_style = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <style>
            @page {{
                size: A4;
                margin: 2cm;
            }}
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                line-height: 1.6;
                color: #333;
                max-width: 210mm;
                margin: 0 auto;
                padding: 20px;
            }}
            h1 {{
                color: #2c3e50;
                border-bottom: 3px solid #3498db;
                padding-bottom: 10px;
                margin-top: 30px;
            }}
            h2 {{
                color: #34495e;
                border-bottom: 2px solid #95a5a6;
                padding-bottom: 8px;
                margin-top: 25px;
            }}
            h3 {{
                color: #7f8c8d;
                margin-top: 20px;
            }}
            h4 {{
                color: #95a5a6;
                margin-top: 15px;
            }}
            table {{
                border-collapse: collapse;
                width: 100%;
                margin: 20px 0;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }}
            th {{
                background-color: #3498db;
                color: white;
                padding: 12px;
                text-align: left;
                border: 1px solid #2980b9;
            }}
            td {{
                padding: 10px;
                border: 1px solid #ddd;
            }}
            tr:nth-child(even) {{
                background-color: #f8f9fa;
            }}
            code {{
                background-color: #f4f4f4;
                padding: 2px 6px;
                border-radius: 3px;
                font-family: 'Courier New', monospace;
                font-size: 0.9em;
            }}
            pre {{
                background-color: #f4f4f4;
                padding: 15px;
                border-radius: 5px;
                overflow-x: auto;
                border-left: 4px solid #3498db;
            }}
            strong {{
                color: #2c3e50;
            }}
            ul, ol {{
                margin: 15px 0;
                padding-left: 30px;
            }}
            li {{
                margin: 8px 0;
            }}
            hr {{
                border: none;
                border-top: 2px solid #ecf0f1;
                margin: 30px 0;
            }}
            blockquote {{
                border-left: 4px solid #3498db;
                padding-left: 20px;
                margin: 20px 0;
                color: #555;
                font-style: italic;
            }}
        </style>
    </head>
    <body>
        {html_content}
    </body>
    </html>
    """

    # Usar reportlab para converter HTML simples para PDF
    try:
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.units import cm
        from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
        from reportlab.lib import colors
        from reportlab.lib.enums import TA_LEFT, TA_CENTER
        import re

        # Criar documento PDF
        doc = SimpleDocTemplate(
            output_pdf,
            pagesize=A4,
            rightMargin=2*cm,
            leftMargin=2*cm,
            topMargin=2*cm,
            bottomMargin=2*cm
        )

        # Estilos
        styles = getSampleStyleSheet()
        styles.add(ParagraphStyle(name='CustomTitle', parent=styles['Heading1'], fontSize=18, textColor=colors.HexColor('#2c3e50'), spaceAfter=12))
        styles.add(ParagraphStyle(name='CustomHeading2', parent=styles['Heading2'], fontSize=14, textColor=colors.HexColor('#34495e'), spaceAfter=10))
        styles.add(ParagraphStyle(name='CustomHeading3', parent=styles['Heading3'], fontSize=12, textColor=colors.HexColor('#7f8c8d'), spaceAfter=8))

        story = []

        # Processar linhas do markdown
        lines = md_content.split('\n')
        i = 0
        while i < len(lines):
            line = lines[i].strip()

            if not line:
                story.append(Spacer(1, 0.3*cm))
            elif line.startswith('# '):
                story.append(Paragraph(line[2:], styles['CustomTitle']))
            elif line.startswith('## '):
                story.append(Spacer(1, 0.3*cm))
                story.append(Paragraph(line[3:], styles['CustomHeading2']))
            elif line.startswith('### '):
                story.append(Paragraph(line[4:], styles['CustomHeading3']))
            elif line.startswith('#### '):
                story.append(Paragraph('<b>' + line[5:] + '</b>', styles['Normal']))
            elif line.startswith('- ') or line.startswith('* '):
                story.append(Paragraph('• ' + line[2:], styles['Normal']))
            elif line.startswith('|') and '|' in line:
                # Tabela
                table_data = []
                while i < len(lines) and lines[i].strip().startswith('|'):
                    row = [cell.strip() for cell in lines[i].strip().split('|')[1:-1]]
                    if not all(set(cell) <= {'-', ' ', ':'} for cell in row):  # Pular linha separadora
                        table_data.append(row)
                    i += 1
                if table_data:
                    t = Table(table_data)
                    t.setStyle(TableStyle([
                        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#3498db')),
                        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                        ('FONTSIZE', (0, 0), (-1, 0), 10),
                        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                        ('GRID', (0, 0), (-1, -1), 1, colors.black)
                    ]))
                    story.append(Spacer(1, 0.3*cm))
                    story.append(t)
                    story.append(Spacer(1, 0.3*cm))
                continue
            elif line == '---':
                story.append(Spacer(1, 0.5*cm))
            else:
                # Texto normal com formatação básica
                text = line.replace('**', '<b>').replace('**', '</b>')
                text = text.replace('`', '<font name="Courier">')
                text = text.replace('`', '</font>')
                story.append(Paragraph(text, styles['Normal']))

            i += 1

        # Gerar PDF
        doc.build(story)
        print(f"✅ PDF gerado com sucesso: {output_pdf}")
        print(f"📄 Tamanho: {os.path.getsize(output_pdf) / 1024:.1f} KB")
        return True

    except Exception as e:
        print(f"❌ Erro ao gerar PDF: {str(e)}")

        # Salvar HTML como fallback
        html_output = Path(md_file).stem + '.html'
        with open(html_output, 'w', encoding='utf-8') as f:
            f.write(html_with_style)
        print(f"\n📄 HTML salvo como alternativa: {html_output}")
        print("   Abra no navegador e use Ctrl+P para salvar como PDF")
        return False


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso: python md_to_pdf.py <arquivo.md> [saida.pdf]")
        sys.exit(1)

    md_file = sys.argv[1]
    output_pdf = sys.argv[2] if len(sys.argv) > 2 else None

    convert_md_to_pdf(md_file, output_pdf)
