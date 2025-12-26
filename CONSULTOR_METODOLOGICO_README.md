# Consultor Metodológico Qualis A1

## 📋 Visão Geral

Este repositório contém uma ferramenta avançada de consultoria metodológica e revisão de artigos científicos, especializada em:

- ✅ Desenho de pesquisa
- ✅ Argumentação metodológica
- ✅ Revisão de introduções acadêmicas
- ✅ Análise de clareza, coerência e progressão lógica
- ✅ Avaliação de contribuição teórica

**Orientado por padrões internacionais de publicação Qualis A1**

---

## 🎯 Funcionalidades Principais

O consultor metodológico executa **7 tarefas especializadas** (A-G):

### Tarefa A: Justificativa Metodológica Convincente
Gera justificativa metodológica de nível A1 cobrindo:
- Alinhamento lógico (problema → método → evidências)
- Adequação ao fenômeno
- Unidade de análise e contexto
- Rigor e qualidade
- Limitações e trade-offs
- Alternativas plausíveis

**Entrega:** Versão longa (500-900 palavras) + versão curta (150-250 palavras)

### Tarefa B: Contexto Específico
Explica por que o contexto empírico escolhido é eficaz:
- Pertinência empírica
- Força inferencial
- Critérios de seleção
- Condições de acesso e integridade

**Entrega:** Parágrafo publicável (120-200 palavras) + bullet-list para defesa oral

### Tarefa C: Diagnóstico de Irrelevâncias
Analisa a introdução identificando trechos irrelevantes:
- Objetivo retórico esperado
- Por que o trecho falha
- Ação recomendada (remover/condensar/mover/reescrever)

**Entrega:** Lista numerada com análise parágrafo a parágrafo

### Tarefa D: Verificação de Progressão Lógica
Verifica progressão clara entre:
1. Apresentação do tema
2. Panorama do debate
3. Lacuna/contradição
4. Problema e pergunta
5. Objetivos
6. Contribuições e estrutura

**Entrega:** Mapa da introdução + recomendações de reordenação

### Tarefa E: Checklist de Elementos Obrigatórios
Verifica presença e clareza de:
- Apresentação do tema
- Panorama (estado do debate)
- Lacuna (gap)
- Pergunta de pesquisa
- Objetivos (geral e específicos)

**Entrega:** Tabela "Elemento | Presente? | Evidência | Ajuste necessário"

### Tarefa F: Reescrita dos Primeiros Parágrafos
Reescreve os 2-4 primeiros parágrafos **SEM alterar referências**:
- Melhora coesão, progressão, definições
- Elimina adjetivação vazia
- Cria pontes para lacuna

**Entrega:** Parágrafos reescritos + lista de operações textuais realizadas

### Tarefa G: Tabela Comparativa de Definições
Cria tabela comparativa do conceito central:
- Autor(es)
- Definição/ênfase central
- Elementos constitutivos
- Implicações operacionais
- Convergências e divergências

**Entrega:** Tabela pronta para colar no artigo (markdown)

---

## 🚀 Instalação

### Pré-requisitos
- Python 3.9+
- pip

### Instalar Dependências

```bash
# Clone o repositório
git clone https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers.git
cd Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers

# Instalar (não há dependências externas além de Python padrão)
chmod +x consultor_metodologico.py
```

---

## 📖 Uso

### Modo 1: Arquivo de Insumos (Recomendado)

```bash
# Executar com arquivo JSON
python consultor_metodologico.py --insumos exemplo_insumos_consultor.json --output relatorio.md
```

### Modo 2: Interativo

```bash
# Modo interativo (solicita insumos via prompt)
python consultor_metodologico.py --interativo
```

### Modo 3: Tarefas Específicas

```bash
# Executar apenas tarefas A, B e C
python consultor_metodologico.py --insumos dados.json --tarefas A,B,C

# Executar todas as tarefas (padrão)
python consultor_metodologico.py --insumos dados.json --tarefas all
```

---

## 📝 Formato dos Insumos

Crie um arquivo JSON com os seguintes campos:

```json
{
  "pergunta_pesquisa": "Sua pergunta de pesquisa aqui...",
  
  "objetivo_geral": "Objetivo geral do estudo...",
  
  "objetivos_especificos": [
    "Objetivo específico 1",
    "Objetivo específico 2",
    "Objetivo específico 3"
  ],
  
  "delimitacao_contexto": "Local, período, população, unidade de análise...",
  
  "estrategia_metodologica": "Tipo de estudo, abordagem, método, técnicas, amostra, dados, análise...",
  
  "introducao_completa": "Texto completo da introdução...",
  
  "referencias_citadas": [
    "AUTOR1, A. Título. Periódico, v. X, n. Y, p. Z, ano.",
    "AUTOR2, B. Título. Periódico, v. X, n. Y, p. Z, ano."
  ],
  
  "conceito_central": "Nome do conceito principal",
  
  "trechos_definicao": [
    {
      "autor": "AUTOR1 (ano)",
      "definicao": "Definição do conceito...",
      "elementos": "Elementos constitutivos...",
      "implicacoes": "Implicações operacionais..."
    }
  ]
}
```

**📄 Arquivo de Exemplo:** Ver `exemplo_insumos_consultor.json`

---

## 📊 Saída Gerada

O consultor gera um relatório completo em Markdown com:

```
relatorio_metodologico.md
├── Sumário Executivo
├── Tarefa A: Justificativa Metodológica (longa + curta)
├── Tarefa B: Contexto Específico (parágrafo + bullets)
├── Tarefa C: Irrelevâncias/trechos fracos na introdução
├── Tarefa D: Progressão lógica (mapa parágrafo a parágrafo)
├── Tarefa E: Checklist dos elementos obrigatórios
├── Tarefa F: Reescrita dos primeiros parágrafos
├── Tarefa G: Tabela comparativa de definições
└── Recomendações Finais (prioridades alta/média/baixa)
```

---

## 🎓 Casos de Uso

### 1. Preparação de Artigo para Qualis A1

```bash
# Use o consultor para revisar sua introdução antes da submissão
python consultor_metodologico.py --insumos meu_artigo.json --output revisao_pre_submissao.md
```

### 2. Defesa de Tese/Dissertação

```bash
# Gere justificativa metodológica robusta para defesa
python consultor_metodologico.py --insumos tese.json --tarefas A,B
```

### 3. Revisão de Projeto de Pesquisa

```bash
# Verifique completude e progressão lógica
python consultor_metodologico.py --insumos projeto.json --tarefas D,E
```

### 4. Análise Conceitual

```bash
# Crie tabela comparativa de definições
python consultor_metodologico.py --insumos conceitos.json --tarefas G
```

---

## 🔍 Exemplos de Análise

### Exemplo 1: Diagnóstico de Irrelevâncias

**Entrada (Parágrafo):**
```
"O tema é muito importante. Vários autores estudam isso. 
É relevante investigar porque pode contribuir."
```

**Saída do Consultor:**
```
❌ Problema: Genérico, sem substância, adjetivação vazia
✅ Ação: Reescrever com especificidade
📝 Sugestão: "Este tema tem recebido atenção crescente 
(AUTOR1, 2020; AUTOR2, 2021), especialmente devido a 
[contexto específico]. Entretanto, aspectos X e Y 
permanecem subinvestigados..."
```

### Exemplo 2: Verificação de Progressão

**Entrada (Introdução):**
```
P1: Tema
P2: Objetivos  ← PROBLEMA: Objetivos antes de lacuna!
P3: Revisão
P4: Lacuna
```

**Saída do Consultor:**
```
⚠️ Salto Lógico Detectado: Objetivos (P2) aparecem 
antes da identificação da lacuna (P4).

✅ Recomendação: Reordenar para P1 → P3 → P4 → P2
Justificativa: Modelo CARS (Swales, 1990) estabelece 
progressão Território → Nicho → Ocupação.
```

### Exemplo 3: Justificativa Metodológica

**Entrada (Estratégia):**
```
"Estudo experimental com 4 datasets, 5 tipos de ruído, 
ANOVA multifatorial"
```

**Saída do Consultor (trecho):**
```
A estratégia experimental com desenho fatorial completo 
é superior a abordagens one-factor-at-a-time porque:

1. Alinhamento: Permite testar hipóteses sobre interações 
   (H₃: Ansatz × NoiseType)
2. Eficiência: Um experimento testa múltiplas hipóteses
3. Rigor: ANOVA quantifica efeitos principais E interações

Alternativas consideradas:
- Meta-análise: Descartada (estudos prévios insuficientes)
- Simulação Monte Carlo: Descartada (menor controle sobre fatores)
```

---

## ⚙️ Configuração Avançada

### Modificar Templates de Análise

Edite diretamente o arquivo `consultor_metodologico.py`:

```python
# Linha ~350: Modificar critérios de irrelevância
def _identificar_problema_paragrafo(self, paragrafo: str, objetivo: str) -> str:
    # Adicione seus critérios customizados aqui
    if "sua_palavra_chave" in paragrafo.lower():
        return "Problema customizado identificado"
```

### Adicionar Novas Tarefas

```python
def tarefa_h_nova_analise(self) -> str:
    """Tarefa H: Nova análise customizada"""
    resultado = "# Tarefa H — Nova Análise\n\n"
    # Implemente sua lógica aqui
    return resultado
```

---

## 📚 Fundamentação Teórica

Este consultor metodológico é baseado em:

### Modelo CARS (Create a Research Space)
- **Swales, J. M. (1990).** *Genre Analysis: English in Academic and Research Settings.*
- Estrutura: Território → Nicho → Ocupação

### Padrões Qualis A1
- Rigor estatístico: ANOVA, post-hoc, effect sizes
- Intervalos de confiança (95%)
- Correção para comparações múltiplas
- Análise de poder estatístico

### Framework SMART para Objetivos
- **S**pecific: Claramente definido
- **M**easurable: Métricas quantitativas
- **A**chievable: Viável
- **R**elevant: Alinhado com lacuna
- **T**ime-bound: Escopo delimitado

---

## 🛠️ Solução de Problemas

### Erro: "INFORMAÇÃO AUSENTE"

**Problema:** Campo obrigatório não fornecido no JSON

**Solução:**
```bash
# Verifique se todos os campos obrigatórios estão preenchidos:
# - pergunta_pesquisa
# - objetivo_geral
# - introducao_completa (para tarefas C-F)
```

### Erro: "Encoding UTF-8"

**Problema:** Caracteres especiais no JSON

**Solução:**
```bash
# Salve o arquivo JSON com encoding UTF-8
# No VS Code: "Save with Encoding" → UTF-8
```

### Saída Genérica

**Problema:** Análise muito genérica

**Solução:**
```bash
# Forneça mais detalhes nos insumos:
# - Introdução completa (não apenas resumo)
# - Referências citadas (lista completa)
# - Trechos de definição (com contexto)
```

---

## 📞 Suporte e Contribuições

### Reportar Bugs

Abra uma issue no GitHub com:
- Descrição do problema
- Arquivo de insumos (JSON)
- Mensagem de erro completa

### Contribuir

1. Fork o repositório
2. Crie branch (`git checkout -b feature/nova-tarefa`)
3. Commit mudanças (`git commit -am 'Adiciona Tarefa H'`)
4. Push para branch (`git push origin feature/nova-tarefa`)
5. Abra Pull Request

---

## 📄 Licença

Este projeto está licenciado sob a [MIT License](LICENSE).

---

## 🎓 Citação

Se você usar este consultor em sua pesquisa, por favor cite:

```bibtex
@software{consultor_metodologico_2025,
  author = {Claro, Marcelo},
  title = {Consultor Metodológico Qualis A1},
  year = {2025},
  url = {https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers}
}
```

---

## 🌟 Recursos Relacionados

- 📖 [PROMPT PARA AUXILIAR NO ARTIGO.md](PROMPT%20PARA%20AUXILIAR%20NO%20ARTIGO.md) - Mega-prompt completo
- 📂 [artigo_cientifico/](artigo_cientifico/) - Framework de geração de artigo
- 📊 [Exemplo de Insumos](exemplo_insumos_consultor.json) - Arquivo JSON de exemplo

---

## 📧 Contato

**Autor:** Marcelo Claro  
**Repositório:** [Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers](https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers)  
**Issues:** [GitHub Issues](https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers/issues)

---

**Última atualização:** 26 de dezembro de 2025  
**Versão:** 1.0.0
