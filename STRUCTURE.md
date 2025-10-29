# Estrutura do Projeto

```
Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers/
│
├── .git/                                   # Repositório Git
├── .gitignore                              # Arquivos ignorados pelo Git
├── .ruff.toml                              # Configuração do linter Ruff
├── .venv/                                  # Ambiente virtual Python (ignorado)
│
├── LICENSE                                 # Licença MIT
├── README.md                               # 📖 LEIA PRIMEIRO - Documentação principal
├── INSTALL.md                              # 🚀 Guia de instalação
├── STRUCTURE.md                            # 📂 Este arquivo - estrutura do projeto
├── requirements.txt                        # 📦 Dependências Python
│
├── framework_investigativo_completo.py     # 🔬 Framework principal (executável)
│
├── docs/                                   # 📚 Documentação detalhada
│   ├── AUTOMACAO_FRAMEWORK.md              # Guia de automação
│   ├── CHANGELOG_v7.2.md                   # Histórico de mudanças
│   ├── GUIA_RAPIDO_v7.2.md                 # Guia rápido de uso
│   └── RESUMO_EXECUTIVO_v7.2.md            # Resumo executivo
│
├── examples/                               # 💡 Exemplos de uso
│   └── exemplo_uso_programatico.py         # Exemplos com Python
│
└── tools/                                  # 🔧 Scripts auxiliares (obsoletos)
    ├── consolidate_results.py              # Funcionalidade integrada no framework
    └── orchestrate_framework.py            # Funcionalidade integrada no framework

# Gerado após execução:
├── resultados_YYYY-MM-DD_HH-MM-SS/         # 📊 Resultados de cada execução
    ├── README.md                           # Descrição dos resultados
    ├── metadata_orchestrator.json          # Metadados da execução
    ├── resultados_completos_artigo.csv     # Consolidado de todos os experimentos
    ├── comparacao_baselines.csv            # Comparação VQC vs. SVM/RF
    ├── figura*.html                        # Visualizações interativas (9 figuras)
    ├── experimentos_individuais/           # CSVs de cada configuração
    │   ├── exp_00001.csv
    │   ├── exp_00002.csv
    │   └── ...
    ├── circuitos/                          # Circuitos exportados (PNG/QASM)
    ├── barren_plateaus/                    # Análises de Barren Plateaus
    ├── analises_individuais/               # Análises estatísticas (ANOVA, etc.)
    └── otimizacao_bayesiana/               # Resultados Optuna (se --bayes)
```

## 🗂️ Descrição dos Arquivos

### Raiz do Projeto

- **framework_investigativo_completo.py** 
  - Framework principal
  - Execute com: `python framework_investigativo_completo.py`
  - Versão: 7.2

- **requirements.txt**
  - Lista todas as dependências Python
  - Instale com: `pip install -r requirements.txt`

- **README.md**
  - Documentação principal do projeto
  - Leia primeiro para entender o projeto

- **INSTALL.md**
  - Guia passo-a-passo de instalação
  - Troubleshooting comum

- **LICENSE**
  - Licença MIT
  - Código aberto e livre para uso

### Diretório `docs/`

- **AUTOMACAO_FRAMEWORK.md** (~350 linhas)
  - Guia completo de automação integrada (v7.2)
  - Modos de execução (grid/bayes/both)
  - Uso programático
  - Estrutura de resultados

- **CHANGELOG_v7.2.md** (~200 linhas)
  - Histórico de mudanças da v7.2
  - Novos recursos
  - Breaking changes (nenhum)
  - Testes realizados

- **GUIA_RAPIDO_v7.2.md** (~400 linhas)
  - Exemplos práticos de uso
  - Cenários comuns
  - Interpretação de resultados
  - Troubleshooting

- **RESUMO_EXECUTIVO_v7.2.md** (~300 linhas)
  - Visão executiva do projeto
  - Benefícios e impacto
  - Estatísticas de implementação

### Diretório `examples/`

- **exemplo_uso_programatico.py**
  - 5 exemplos práticos de uso via API Python
  - Consolidação de resultados
  - Análise com Pandas
  - Comparação de baselines

### Diretório `tools/` (Obsoleto)

⚠️ **Nota**: Os scripts em `tools/` são mantidos para retrocompatibilidade,
mas suas funcionalidades foram **integradas no framework principal**.

- **consolidate_results.py** → Agora: `consolidar_e_gerar_metadados()` no framework
- **orchestrate_framework.py** → Agora: execução automática no `main()`

**Recomendação**: Use o framework principal diretamente.

## 📊 Resultados Gerados

Após execução, uma pasta `resultados_YYYY-MM-DD_HH-MM-SS/` é criada com:

### Arquivos Principais

| Arquivo | Descrição |
|---------|-----------|
| `resultados_completos_artigo.csv` | Consolidado de todos os experimentos |
| `comparacao_baselines.csv` | VQC vs. SVM/RF por dataset |
| `metadata_orchestrator.json` | Inventário completo e metadados |
| `figura*.html` | 9 visualizações interativas (Plotly) |

### Subpastas

| Pasta | Conteúdo |
|-------|----------|
| `experimentos_individuais/` | CSV de cada configuração (8,280 configs × 5 seeds) |
| `circuitos/` | Circuitos quânticos exportados (PNG/QASM) |
| `barren_plateaus/` | Análises de Barren Plateaus |
| `analises_individuais/` | Análises estatísticas (ANOVA, post-hoc) |
| `otimizacao_bayesiana/` | Resultados Optuna (apenas se --bayes) |

## 🔍 Fluxo de Trabalho Típico

```
1. Instalar dependências
   → pip install -r requirements.txt

2. Executar framework
   → python framework_investigativo_completo.py

3. Aguardar conclusão
   → Logs mostram progresso em tempo real

4. Analisar resultados
   → Pasta resultados_YYYY-MM-DD_HH-MM-SS/ criada automaticamente
   → Abrir visualizações HTML no navegador
   → Carregar CSVs consolidados com Pandas

5. (Opcional) Usar API programática
   → Ver examples/exemplo_uso_programatico.py
```

## 📚 Ordem de Leitura Recomendada

Para novos usuários:

1. **README.md** - Visão geral do projeto
2. **INSTALL.md** - Configurar ambiente
3. **docs/GUIA_RAPIDO_v7.2.md** - Primeiros passos
4. **examples/exemplo_uso_programatico.py** - Exemplos práticos
5. **docs/AUTOMACAO_FRAMEWORK.md** - Recursos avançados

Para desenvolvedores/pesquisadores:

1. **framework_investigativo_completo.py** - Código-fonte principal
2. **docs/CHANGELOG_v7.2.md** - Mudanças recentes
3. **docs/RESUMO_EXECUTIVO_v7.2.md** - Visão técnica

## 🧹 Arquivos Ignorados (.gitignore)

Os seguintes arquivos/pastas são ignorados pelo Git:

- `__pycache__/` - Cache Python
- `.venv/` - Ambiente virtual
- `resultados_*/` - Resultados de execução (gerados localmente)
- `*.csv`, `*.html`, `*.png` - Artefatos de execução
- `*.pdf` - PDFs gerados
- `.vscode/`, `.idea/` - Configurações de IDE

**Motivo**: Esses arquivos são gerados automaticamente ou específicos do ambiente local.

## 🔄 Reprodutibilidade

Para reproduzir resultados em outra máquina:

1. Clone o repositório
2. Instale as dependências (`requirements.txt`)
3. Execute o framework com o mesmo modo
4. Compare `resultados_completos_artigo.csv`

**Seeds fixas**: O framework usa seeds fixas (42-46) para reprodutibilidade.

## 📝 Convenções

- **Markdown**: Documentação em `.md`
- **Python**: Código-fonte em `.py`
- **CSV**: Dados tabulares
- **JSON**: Metadados estruturados
- **HTML**: Visualizações interativas (Plotly)

## 🚀 Versionamento

- **v7.2** (atual): Automação integrada
- **v7.1**: Otimização Bayesiana
- **v7.0**: Framework base

Ver `docs/CHANGELOG_v7.2.md` para histórico completo.

## 📞 Suporte

- **Documentação**: Consulte `docs/`
- **Exemplos**: Ver `examples/`
- **Issues**: GitHub Issues (se aplicável)

---

**Última atualização**: 28 de outubro de 2025  
**Versão do framework**: 7.2
