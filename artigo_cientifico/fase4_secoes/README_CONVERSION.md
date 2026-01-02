# Conversão de Markdown para DOCX e PDF

Este diretório contém os arquivos Markdown das seções do artigo científico (Fase 4) e suas conversões para os formatos DOCX e PDF.

## Arquivos

Todos os 12 arquivos Markdown foram convertidos com sucesso:

| Arquivo Markdown | DOCX | PDF | Tamanho PDF |
|-----------------|------|-----|-------------|
| agradecimentos_referencias.md | ✅ | ✅ | 69 KB |
| conclusao_completa.md | ✅ | ✅ | 70 KB |
| contraprova_casos_limite.md | ✅ | ✅ | 70 KB |
| discussao_completa.md | ✅ | ✅ | 135 KB |
| introducao_completa.md | ✅ | ✅ | 87 KB |
| metodologia_completa.md | ✅ | ✅ | 124 KB |
| prova_teorema.md | ✅ | ✅ | 73 KB |
| resultados_completo.md | ✅ | ✅ | 110 KB |
| resumo_abstract.md | ✅ | ✅ | 43 KB |
| revisao_literatura_completa.md | ✅ | ✅ | 108 KB |
| secao_didatica_leigos.md | ✅ | ✅ | 79 KB |
| teorema_beneficio_condicionado.md | ✅ | ✅ | 91 KB |

**Total:** 12 arquivos MD → 12 DOCX + 12 PDF = 24 arquivos de saída

## Como Usar o Script de Conversão

### Pré-requisitos

O script requer:
- Python 3.x
- Pandoc (instalado via `apt-get install pandoc`)
- LaTeX (texlive-xetex) para geração de PDF

### Instalação das Dependências

```bash
# Instalar pandoc
sudo apt-get update
sudo apt-get install -y pandoc

# Instalar LaTeX para PDFs
sudo apt-get install -y texlive-xetex texlive-fonts-recommended texlive-plain-generic
```

### Execução

```bash
# Converter todos os arquivos .md da pasta fase4_secoes
python3 tools/convert_fase4_to_docx_pdf.py

# Ou especificar outro diretório
python3 tools/convert_fase4_to_docx_pdf.py /caminho/para/diretorio
```

### Características do Script

O script `tools/convert_fase4_to_docx_pdf.py`:

1. **Encontra automaticamente** todos os arquivos `.md` no diretório especificado
2. **Converte para DOCX** usando pandoc com suporte a:
   - Tabelas
   - Blocos de código
   - Formatação Markdown completa
3. **Converte para PDF** usando pandoc + XeLaTeX com:
   - Suporte completo a Unicode
   - Renderização de equações LaTeX (formato `$...$` e `$$...$$`)
   - Fontes DejaVu para melhor compatibilidade
   - Margens de 2cm
4. **Sanitiza o conteúdo** automaticamente:
   - Remove emojis problemáticos para LaTeX
   - Converte checkboxes markdown (`- [x]` → `- [DONE]`)
   - Corrige padrões de fórmulas matemáticas que causam erros
5. **Gera relatório** detalhado com estatísticas de sucesso/falha

## Formato dos Arquivos

### DOCX
- Arquivos editáveis no Microsoft Word, LibreOffice, Google Docs
- Mantém formatação básica (títulos, listas, tabelas, negrito, itálico)
- Tamanho típico: 15-35 KB

### PDF
- Formato final para leitura e distribuição
- Renderização completa de equações matemáticas via LaTeX
- Formatação profissional com fontes serifadas
- Tamanho típico: 40-140 KB dependendo do conteúdo matemático

## Observações Técnicas

### Equações Matemáticas

Os arquivos suportam LaTeX inline e display:

```markdown
Equação inline: $E = mc^2$

Equação em bloco:
$$
\hat{H} = \sum_{i} \sigma_i^z
$$
```

### Tratamento de Caracteres Especiais

O script remove automaticamente:
- Emojis (✅, ❌, ⚠️, etc.) que causam problemas no LaTeX
- Checkboxes markdown são convertidos para texto

### Engine de PDF

O script usa `xelatex` como engine de PDF porque:
- Suporte nativo a Unicode
- Melhor renderização de caracteres especiais
- Compatibilidade com fontes TrueType/OpenType

## Manutenção

Para reconverter todos os arquivos (por exemplo, após atualizar o Markdown):

```bash
# Remover conversões antigas
rm artigo_cientifico/fase4_secoes/*.docx artigo_cientifico/fase4_secoes/*.pdf

# Reconverter
python3 tools/convert_fase4_to_docx_pdf.py
```

## Solução de Problemas

### Erro: "pandoc: command not found"
```bash
sudo apt-get install pandoc
```

### Erro: "xelatex not found"
```bash
sudo apt-get install texlive-xetex texlive-fonts-recommended
```

### PDF não gera corretamente
- Verifique se há caracteres especiais não suportados no Markdown
- O script já trata os casos mais comuns, mas novos padrões podem precisar de ajustes
- Verifique os logs de erro do pandoc para detalhes

## Estrutura de Arquivos

```
artigo_cientifico/fase4_secoes/
├── *.md                    # Arquivos Markdown originais
├── *.docx                  # Conversões DOCX
├── *.pdf                   # Conversões PDF
└── README_CONVERSION.md    # Esta documentação
```

## Contato

Para problemas ou sugestões relacionadas à conversão, consulte o script em:
`tools/convert_fase4_to_docx_pdf.py`

---

**Última atualização:** 02 de janeiro de 2026  
**Versão do script:** 1.0  
**Status:** Todos os arquivos convertidos com sucesso ✅
