# 🚀 Guia Rápido: Consultor Metodológico Qualis A1

## Início Rápido em 3 Passos

### 1️⃣ Preparar Insumos

Crie um arquivo JSON com seus dados de pesquisa:

```bash
cp exemplo_insumos_consultor.json meu_artigo.json
# Edite meu_artigo.json com seus dados
```

### 2️⃣ Executar Análise

```bash
# Usando o wrapper (recomendado)
./executar_consultor.sh --insumos meu_artigo.json --output relatorio.md

# Ou diretamente
python consultor_metodologico.py --insumos meu_artigo.json --output relatorio.md
```

### 3️⃣ Revisar Relatório

```bash
# Abra o relatório gerado
cat relatorio.md
# ou use seu editor preferido
code relatorio.md
```

---

## 💡 Casos de Uso Comuns

### Revisar Introdução Antes da Submissão

```bash
# Executar todas as tarefas (A-G)
./executar_consultor.sh --insumos introducao.json
```

### Gerar Apenas Justificativa Metodológica

```bash
# Executar apenas Tarefa A
./executar_consultor.sh --insumos dados.json --tarefas A --output justificativa.md
```

### Verificar Progressão Lógica

```bash
# Executar Tarefas D e E
./executar_consultor.sh --insumos introducao.json --tarefas D,E
```

### Criar Tabela de Definições

```bash
# Executar apenas Tarefa G
./executar_consultor.sh --insumos conceitos.json --tarefas G
```

---

## 📝 Template de Insumos Mínimo

```json
{
  "pergunta_pesquisa": "Qual é a sua pergunta de pesquisa?",
  "objetivo_geral": "Qual é o objetivo geral do estudo?",
  "objetivos_especificos": [
    "Objetivo específico 1",
    "Objetivo específico 2"
  ],
  "delimitacao_contexto": "Descreva o contexto empírico",
  "estrategia_metodologica": "Descreva a estratégia metodológica",
  "introducao_completa": "Cole o texto completo da introdução aqui...",
  "referencias_citadas": [
    "Autor1. Título. Periódico, ano.",
    "Autor2. Título. Periódico, ano."
  ],
  "conceito_central": "Nome do conceito principal",
  "trechos_definicao": []
}
```

Salve como `meu_artigo.json` e execute!

---

## 🔍 Exemplos de Saída

### Tarefa A: Justificativa Metodológica

```markdown
# Tarefa A — Justificativa Metodológica

## Versão Longa (500-900 palavras)

### 1. Alinhamento Lógico
A estratégia metodológica adotada alinha-se com a pergunta de pesquisa...

### 2. Adequação ao Fenômeno
O método escolhido é adequado porque...

[... 6 seções completas ...]

## Versão Curta (150-250 palavras)
Este estudo adota [tipo] para [objetivo]...
```

### Tarefa C: Diagnóstico de Irrelevâncias

```markdown
# Tarefa C — Diagnóstico de Irrelevâncias

### Parágrafo 3

**Trecho:** "O tema é importante..."

**Problema Identificado:** Genérico, sem substância

**Ação Recomendada:** Reescrever com especificidade

**Justificativa:** Adjetivação vazia prejudica rigor A1
```

### Tarefa E: Checklist

```markdown
# Tarefa E — Checklist de Elementos

| Elemento | Presente? | Evidência | Ajuste |
|----------|-----------|-----------|--------|
| Tema | ✅ Sim | Parágrafo 1 | Sem ajustes |
| Panorama | ⚠️ Parcial | Parágrafos 2-4 | Estruturar |
| Lacuna | ❌ Não | - | Adicionar explicitamente |
```

---

## ⚙️ Opções Avançadas

### Modo Interativo

```bash
./executar_consultor.sh --interativo
```

Você será solicitado a fornecer cada insumo via prompt.

### Executar Tarefas Específicas

```bash
# Apenas A e B
./executar_consultor.sh --insumos dados.json --tarefas A,B

# Apenas C, D e E
./executar_consultor.sh --insumos dados.json --tarefas C,D,E

# Todas (padrão)
./executar_consultor.sh --insumos dados.json --tarefas all
```

### Salvar em Arquivo Específico

```bash
./executar_consultor.sh --insumos dados.json --output ~/Desktop/analise.md
```

---

## 🐛 Solução de Problemas

### Erro: "Python não encontrado"

**Solução:**
```bash
# Ubuntu/Debian
sudo apt install python3

# macOS
brew install python3

# Verificar instalação
python3 --version
```

### Erro: "INFORMAÇÃO AUSENTE"

**Solução:**
Verifique se todos os campos obrigatórios estão preenchidos no JSON:
- `pergunta_pesquisa`
- `objetivo_geral`
- `estrategia_metodologica`

### Saída muito genérica

**Solução:**
Forneça mais detalhes nos insumos:
- Introdução completa (não apenas resumo)
- Lista completa de referências citadas
- Contexto empírico detalhado

---

## 📚 Documentação Completa

Para documentação detalhada, consulte:

- 📖 [CONSULTOR_METODOLOGICO_README.md](CONSULTOR_METODOLOGICO_README.md) - Documentação completa
- 📄 [exemplo_insumos_consultor.json](exemplo_insumos_consultor.json) - Exemplo de entrada
- 💻 [consultor_metodologico.py](consultor_metodologico.py) - Código-fonte

---

## 🎯 Próximos Passos

Depois de executar o consultor:

1. ✅ Revise o relatório gerado
2. ✅ Implemente as recomendações prioritárias
3. ✅ Reescreva parágrafos conforme sugerido (Tarefa F)
4. ✅ Complete lacunas de citação identificadas
5. ✅ Reorganize parágrafos se necessário (Tarefa D)
6. ✅ Execute novamente para validar melhorias

---

## 💬 Feedback

Encontrou problemas ou tem sugestões?

- 🐛 [Reportar bug](https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers/issues)
- 💡 [Sugerir melhoria](https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers/issues)

---

**Última atualização:** 26/12/2025
