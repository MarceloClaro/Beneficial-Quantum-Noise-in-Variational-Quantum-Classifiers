# Instruções para Uso no Overleaf

## 📄 Arquivo Criado

**Nome:** `artigo_completo_qualis_a1.tex`
**Tamanho:** ~321 KB (318.856 caracteres)
**Conteúdo:** Artigo completo com ~21.400 palavras

---

## 🚀 Como Usar no Overleaf.com

### Passo 1: Fazer Upload

1. Acesse [Overleaf.com](https://www.overleaf.com)
2. Faça login ou crie uma conta gratuita
3. Clique em **"New Project"** → **"Upload Project"**
4. Faça upload do arquivo `artigo_completo_qualis_a1.tex`

### Passo 2: Compilar

1. O Overleaf compilará automaticamente
2. Caso ocorra erro, clique em **"Recompile"**
3. Verifique o PDF gerado no painel direito

### Passo 3: Ajustes (se necessário)

Se houver erros de compilação:

#### Problemas Comuns

**Erro de encoding:**
- Adicione no início: `\usepackage[utf8]{inputenc}`
- Já está incluído no arquivo

**Erro com matemática:**
- Verifique se pacotes `amsmath`, `amssymb` estão carregados
- Já estão incluídos

**Tabelas muito largas:**
- Use `\small` antes da tabela
- Ajuste com `\resizebox{\textwidth}{!}{...}`

**Figuras faltando:**
- As figuras precisam ser carregadas separadamente
- Comentadas com `% FIGURA AQUI` no texto

---

## 📋 Estrutura do Documento

### Seções Principais

1. **Capa** (Titlepage)
2. **Resumo/Abstract** (Português e Inglês)
3. **Sumário** (Automático)
4. **Introdução** (~3.000 palavras)
5. **Revisão da Literatura**
6. **Teorema do Benefício Condicionado** (~3.400 palavras)
   - 3 Lemas com provas
7. **Prova do Teorema** (~2.900 palavras)
   - 3 passos principais
8. **Contraprova e Casos-Limite** (~2.500 palavras)
9. **Metodologia** (~1.500 palavras)
10. **Resultados** (~3.000 palavras)
11. **Discussão** (~2.500 palavras)
12. **Seção Didática para Leigos** (~1.500 palavras)
13. **Conclusão** (~800 palavras)

### Apêndices

- **Apêndice D:** Métrica de Fubini-Study (~1.100 palavras)
- **Apêndice E:** Framework AUEC (~1.250 palavras)
- **Apêndice F:** Barren Plateaus (~1.050 palavras)
- **Apêndice G:** Validação Estatística (~1.300 palavras)
- **Apêndice I:** Lista de Símbolos (~550 palavras)
- **Apêndice J:** Checklist de Verificação (~550 palavras)

### Elementos Finais

- Agradecimentos
- Disponibilidade de Dados
- Conflito de Interesses

---

## 🎨 Personalização

### Alterar Autores

Localize na linha ~43:

```latex
{\large Equipe de Pesquisa em Computação Quântica\par}
```

Substitua por:

```latex
{\large Seu Nome\par}
{\large Afiliação Institucional\par}
```

### Alterar Título (se necessário)

Localize:

```latex
{\Huge\bfseries Do Obstáculo à Oportunidade:\par}
```

### Adicionar Figuras

Substitua comentários `% FIGURA AQUI` por:

```latex
\begin{figure}[h]
    \centering
    \includegraphics[width=0.8\textwidth]{nome_figura.pdf}
    \caption{Legenda da figura}
    \label{fig:nome}
\end{figure}
```

### Adicionar Referências

No final, antes de `\end{document}`, adicione:

```latex
\bibliographystyle{plain}
\bibliography{referencias}
```

E crie arquivo `referencias.bib` no Overleaf.

---

## ✅ Checklist de Qualidade

Após compilar, verifique:

- [ ] PDF gerado sem erros
- [ ] Sumário funcionando (clicável)
- [ ] Todas as seções presentes
- [ ] Equações renderizadas corretamente
- [ ] Tabelas formatadas adequadamente
- [ ] Numeração de páginas correta

---

## 📊 Estatísticas do Documento

| Métrica | Valor |
|---------|-------|
| **Palavras Totais** | ~21.400 |
| **Páginas Estimadas** | ~60-70 (A4, 12pt) |
| **Equações Numeradas** | 127+ |
| **Tabelas** | 35+ |
| **Teoremas/Lemas** | 4 |
| **Provas Completas** | 4 |
| **Apêndices** | 6 novos |

---

## 🔧 Problemas e Soluções

### Compilação Lenta

**Causa:** Documento muito grande
**Solução:** 
- Use fast compile no Overleaf
- Comente seções desnecessárias temporariamente

### Erros de Matemática

**Sintoma:** `Missing $ inserted`
**Solução:** Verifique balanceamento de `$...$` e `\[...\]`

### Espaçamento Estranho

**Sintoma:** Grandes espaços em branco
**Solução:** Use `\raggedbottom` no preâmbulo (já incluído)

---

## 📞 Suporte

Se encontrar problemas:

1. Verifique o log de compilação (ícone de alerta no Overleaf)
2. Consulte [Overleaf Knowledge Base](https://www.overleaf.com/learn)
3. Ou ajuste manualmente as seções problemáticas

---

## 🎯 Próximos Passos Recomendados

1. **Revisar matemática:** Verificar todas as equações
2. **Adicionar figuras:** Gerar e inserir as 8+ figuras mencionadas
3. **Completar referências:** Adicionar bibliografia completa
4. **Revisão de português:** Verificar ortografia e gramática
5. **Formatação final:** Ajustar espaçamentos e quebras de página

---

## ✨ Recursos Adicionais

### Templates Overleaf Úteis

- IEEE Template: https://www.overleaf.com/latex/templates/ieee-conference-template
- Springer Template: https://www.overleaf.com/latex/templates/springer-latex-template

### Documentação LaTeX

- Wikibooks LaTeX: https://en.wikibooks.org/wiki/LaTeX
- Overleaf Tutorials: https://www.overleaf.com/learn/latex/Tutorials

---

**Versão:** 1.0  
**Data:** 02/01/2026  
**Status:** ✅ Pronto para Overleaf

**Bom trabalho! 🎉**
