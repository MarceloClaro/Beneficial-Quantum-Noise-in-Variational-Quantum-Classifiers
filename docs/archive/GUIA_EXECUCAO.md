# Guia de Execução do Framework Investigativo v7.2

## 🚀 Execução Rápida

### Método 1: Script Automático (Recomendado)

```bash

# Linux/macOS
./executar_framework.sh

# Windows (PowerShell)
bash executar_framework.sh

```text

O script oferece um menu interativo com as seguintes opções:

1. **Modo Rápido Bayesiano** (~15 minutos) - Ideal para testes
2. **Modo Bayesiano Completo** (~1-2 horas) - Otimização inteligente
3. **Modo Grid Search Rápido** (~5-6 horas) - Exploração básica
4. **Modo Grid Search Completo** (~15-20 horas) - Exploração exhaustiva
5. **Modo Híbrido** (~20-25 horas) - Grid + Bayesiano
6. **Modo Personalizado** - Você define os parâmetros


### Método 2: Comando Direto

#### Modo Rápido Bayesiano (Teste)

```bash
export VQC_QUICK=1
export VQC_BAYESIAN=1
python framework_investigativo_completo.py --bayes --trials 10 --dataset-bayes moons

```text

#### Modo Bayesiano Completo

```bash
export VQC_BAYESIAN=1
python framework_investigativo_completo.py --bayes --trials 200 --dataset-bayes all

```text

#### Modo Grid Search Rápido

```bash
export VQC_QUICK=1
python framework_investigativo_completo.py

```text

#### Modo Grid Search Completo

```bash
python framework_investigativo_completo.py

```text

#### Modo Híbrido (Grid + Bayesiano)

```bash
python framework_investigativo_completo.py --run-both

```text

## 📋 Pré-requisitos

- Python 3.9 ou superior
- pip (gerenciador de pacotes Python)
- 8 GB RAM mínimo (16 GB recomendado)
- Espaço em disco: ~1 GB para resultados


## 📦 Instalação de Dependências

```bash
pip install -r requirements.txt

```text

Principais pacotes instalados:

- PennyLane >= 0.30.0 (computação quântica)
- Optuna >= 3.0.0 (otimização Bayesiana)
- NumPy >= 1.23.0 (arrays e computação numérica)
- Pandas >= 2.0.0 (manipulação de dados)
- Scipy >= 1.10.0 (algoritmos científicos)
- Plotly >= 5.0.0 (visualizações interativas)
- Matplotlib >= 3.5.0 (visualizações estáticas)
- scikit-learn >= 1.3.0 (machine learning)
- statsmodels >= 0.14.0 (análises estatísticas)
- joblib >= 1.2.0 (paralelização)
- kaleido >= 0.2.1 (exportação de imagens)


## ⚙️ Parâmetros de Execução

### Variáveis de Ambiente

| Variável | Valores | Descrição |
|----------|---------|-----------|
| `VQC_QUICK` | 0 ou 1 | Modo rápido (5 épocas vs 15) |
| `VQC_BAYESIAN` | 0 ou 1 | Usar otimização Bayesiana |
| `VQC_BAYES_AFTER_GRID` | 0 ou 1 | Executar Bayesiano após Grid |
| `RESULTS_BASE_DIR` | caminho | Diretório base para resultados |

### Argumentos CLI

| Argumento | Descrição | Exemplo |
|-----------|-----------|---------|
| `--bayes` | Ativar modo Bayesiano | `--bayes` |
| `--trials N` | Número de trials Bayesianos | `--trials 100` |
| `--dataset-bayes NOME` | Dataset para Bayesiano | `--dataset-bayes moons` |
| `--epocas-bayes N` | Épocas por trial | `--epocas-bayes 20` |
| `--bayes-after-grid` | Grid seguido de Bayesiano | `--bayes-after-grid` |
| `--run-both` | Mesmo que `--bayes-after-grid` | `--run-both` |

## 📊 Resultados Gerados

Após a execução, um diretório `resultados_YYYY-MM-DD_HH-MM-SS/` é criado contendo:

### Estrutura de Diretórios

```

resultados_2025-12-23_13-39-53/
├── README.md                              # Descrição do experimento
├── metadata.json                          # Metadados raiz
├── metadata_orchestrator.json             # Metadados de consolidação
├── resultados_completos_artigo.csv        # Resultados consolidados
├── comparacao_baselines.csv               # VQC vs SVM/RF
├── analise_comparacao_inicializacoes.csv  # Análise de inicializações
├── analises_estatisticas_completo.csv     # Análises estatísticas
│
├── figura2_beneficial_noise.html          # Visualizações interativas
├── figura2b_beneficial_noise_ic95.html
├── figura3_noise_types.html
├── figura3b_noise_types_ic95.html
├── figura4_initialization.html
├── figura5_architecture_tradeoffs.html
├── figura6_effect_sizes.html
├── figura7_overfitting.html
├── figura_correlacao.html
│
├── otimizacao_bayesiana/                  # Resultados Bayesianos
│   ├── resultado_otimizacao.json
│   ├── historico_trials.csv
│   └── README_otimizacao.md
│
├── circuitos/                             # Circuitos quânticos exportados
├── barren_plateaus/                       # Análise de Barren Plateaus
├── experimentos_individuais/              # CSVs por experimento
├── analises_individuais/                  # Análises granulares
└── visualizacoes_individuais/             # Visualizações granulares

```text

### Arquivos Principais

#### 1. `resultados_completos_artigo.csv`
Dados consolidados de todos os experimentos com colunas:

- `dataset`, `arquitetura`, `estrategia_init`
- `tipo_ruido`, `nivel_ruido`
- `acuracia_treino`, `acuracia_teste`
- `gap_treino_teste`, `tempo_segundos`
- `n_parametros`, `entropia_final`, `negatividade_media`
- `seed`


#### 2. `comparacao_baselines.csv`
Comparação VQC vs. classificadores clássicos (SVM, Random Forest)

#### 3. Visualizações HTML
Figuras interativas Plotly que podem ser abertas em qualquer navegador

#### 4. `otimizacao_bayesiana/resultado_otimizacao.json`
Contém:

- Melhor acurácia encontrada
- Melhores hiperparâmetros
- Importância de cada hiperparâmetro
- Histórico completo de trials


## 🔍 Monitoramento da Execução

### Acompanhar Progresso em Tempo Real

```bash

# Ver últimas linhas do log
tail -f framework.log

# Contar experimentos concluídos
grep "✓ Acurácia" framework.log | wc -l

# Verificar erros
grep "ERROR\|Traceback" framework.log

```text

### Estimar Tempo Restante

O framework exibe progresso como: `[X/Y]` onde:

- `X` = experimentos concluídos
- `Y` = total de experimentos


Exemplo:

```

[842/8280] Dataset: moons | Acurácia: 0.8250 | Tempo: 120.5s

```text

Tempo estimado restante: `(Y - X) * tempo_médio`

## ⏱️ Tempo de Execução Estimado

| Modo | Configurações | Tempo Estimado | Uso Recomendado |
|------|--------------|----------------|-----------------|
| **Rápido Bayesiano** | 10 trials × 5 épocas | 15-30 min | Testes rápidos |
| **Bayesiano Completo** | 200 trials × 15 épocas | 1-2 horas | Otimização eficiente |
| **Grid Rápido** | 8,280 × 5 épocas | 5-6 horas | Exploração básica |
| **Grid Completo** | 8,280 × 15 épocas | 15-20 horas | Análise exhaustiva |
| **Híbrido** | Grid + Bayesiano | 20-25 horas | Máxima precisão |

## 🔧 Troubleshooting

### Erro: "PennyLane not found"

```bash
pip install pennylane>=0.30.0

```text

### Erro: "Optuna not available"

```bash
pip install optuna>=3.0.0

```text

### Memória Insuficiente
- Use o modo rápido: `export VQC_QUICK=1`
- Reduza o número de trials: `--trials 50`
- Feche outros programas


### Execução Muito Lenta
- Use modo Bayesiano em vez de Grid Search
- Reduza o número de épocas: `--epocas-bayes 5`
- Execute para um dataset específico: `--dataset-bayes moons`


### Resultados não consolidados
O framework consolida automaticamente ao final. Se houver erro, execute manualmente:

```python
from framework_investigativo_completo import consolidar_e_gerar_metadados

resultado = consolidar_e_gerar_metadados(
    pasta_resultados="resultados_2025-12-23_13-39-53",
    verbose=True
)

```text

## 📝 Exemplos de Uso

### Exemplo 1: Teste Rápido (15 minutos)

```bash
export VQC_QUICK=1
export VQC_BAYESIAN=1
python framework_investigativo_completo.py --bayes --trials 5 --dataset-bayes moons

```text

### Exemplo 2: Análise Completa de Um Dataset (2 horas)

```bash
python framework_investigativo_completo.py --bayes --trials 150 --epocas-bayes 20 --dataset-bayes iris

```text

### Exemplo 3: Exploração Completa Todos os Datasets (4 horas)

```bash
python framework_investigativo_completo.py --bayes --trials 200 --dataset-bayes all

```text

### Exemplo 4: Grid Search Completo (Artigo Científico - 20 horas)

```bash
python framework_investigativo_completo.py

```text

### Exemplo 5: Híbrido - Máxima Precisão (25 horas)

```bash
python framework_investigativo_completo.py --run-both --trials 150

```

## 🎯 Escolhendo o Modo Certo

### Para Desenvolvimento e Testes
➡️ **Modo Rápido Bayesiano** (15-30 min)

- Valida que tudo está funcionando
- Testa mudanças no código
- Explora rapidamente hiperparâmetros


### Para Artigos Científicos (Qualis A1)
➡️ **Modo Grid Completo** (15-20 horas)

- Cobertura exhaustiva do espaço de hiperparâmetros
- Análises estatísticas rigorosas
- Máxima reprodutibilidade


### Para Otimização Eficiente
➡️ **Modo Bayesiano Completo** (1-2 horas)

- 10-20x mais rápido que Grid Search
- Encontra configurações ótimas inteligentemente
- Ideal para pesquisa aplicada


### Para Máxima Precisão
➡️ **Modo Híbrido** (20-25 horas)

- Combina exploração (Grid) com refinamento (Bayesiano)
- Melhor de ambos os mundos
- Recomendado para trabalhos definitivos


## 💡 Dicas e Boas Práticas

1. **Sempre comece com modo rápido** para validar a configuração
2. **Use modo Bayesiano** para exploração inicial eficiente
3. **Reserve Grid Completo** para o experimento final do artigo
4. **Monitore o progresso** com `tail -f framework.log`
5. **Faça backup dos resultados** importantes
6. **Use nomes descritivos** ao salvar resultados importantes
7. **Verifique os metadados** para validar a execução


## 📚 Documentação Adicional

- [README.md](README.md) - Visão geral do projeto
- [INSTALL.md](INSTALL.md) - Guia de instalação detalhado
- [STRUCTURE.md](STRUCTURE.md) - Estrutura do código
- [docs/AUTOMACAO_FRAMEWORK.md](docs/AUTOMACAO_FRAMEWORK.md) - Automação v7.2
- [docs/GUIA_RAPIDO_v7.2.md](docs/GUIA_RAPIDO_v7.2.md) - Guia rápido
- [examples/exemplo_uso_programatico.py](examples/exemplo_uso_programatico.py) - Exemplos de uso


## 🆘 Suporte

Para dúvidas, problemas ou sugestões:

1. Verifique este guia primeiro
2. Consulte a documentação no diretório `docs/`
3. Abra uma issue no repositório GitHub
4. Entre em contato: marceloclaro@gmail.com


## 📄 Licença

MIT License - Copyright (c) 2025 Marcelo Claro Laranjeira

---


**⭐ Se este framework foi útil para sua pesquisa, considere citar nosso trabalho e dar uma estrela no repositório!**

