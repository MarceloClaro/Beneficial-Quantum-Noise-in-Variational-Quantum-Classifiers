# Tabela de Rastreabilidade: Código → Método

Esta tabela mapeia cada componente metodológico descrito no artigo para sua implementação exata no código-fonte, garantindo rastreabilidade total e reprodutibilidade.

## Instruções de Preenchimento

Para cada linha:
1. **Componente do Método**: Nome do elemento metodológico como descrito no artigo
2. **Arquivo/Função/Linha**: Localização exata no código (formato: `arquivo.py:função():L15-30`)
3. **Parâmetros**: Valores ou variáveis usadas (com valores padrão quando aplicável)
4. **Artefatos Gerados**: Outputs produzidos (arquivos, objetos, métricas)

## Tabela Principal

| ID | Componente do Método | Arquivo/Função/Linha | Parâmetros | Artefatos Gerados |
|----|---------------------|----------------------|------------|-------------------|
| M01 | Definição do Ansatz BasicEntangler | `framework_investigativo_completo.py:criar_circuito_ansatz():L245-260` | `n_qubits=4`, `depth=2`, `ansatz_name='BasicEntangler'` | Objeto `qml.QNode` |
| M02 | Definição do Ansatz StronglyEntangling | `framework_investigativo_completo.py:criar_circuito_ansatz():L262-278` | `n_qubits=4`, `depth=3`, `ansatz_name='StronglyEntangling'` | Objeto `qml.QNode` |
| M03 | Canal de Ruído Depolarizing | `framework_investigativo_completo.py:aplicar_ruido():L320-335` | `p=0.01`, `noise_type='Depolarizing'`, `n_qubits=4` | Operadores de Kraus $(K_0, K_1, K_2, K_3)$ |
| M04 | Canal de Ruído Amplitude Damping | `framework_investigativo_completo.py:aplicar_ruido():L337-348` | `gamma=0.05`, `noise_type='AmplitudeDamping'` | Operadores de Kraus $(K_0, K_1)$ |
| M05 | Canal de Ruído Phase Damping | `framework_investigativo_completo.py:aplicar_ruido():L350-361` | `gamma=0.03`, `noise_type='PhaseDamping'` | Operadores de Kraus $(K_0, K_1)$ |
| M06 | Schedule de Ruído Constant | `framework_investigativo_completo.py:calcular_schedule_ruido():L390-395` | `p_inicial=0.1`, `epoca=t`, `total_epocas=100`, `schedule='Constant'` | $p(t) = p_0$ |
| M07 | Schedule de Ruído Linear | `framework_investigativo_completo.py:calcular_schedule_ruido():L397-402` | `p_inicial=0.1`, `schedule='Linear'` | $p(t) = p_0 \cdot (1 - t/T)$ |
| M08 | Schedule de Ruído Cosine | `framework_investigativo_completo.py:calcular_schedule_ruido():L404-409` | `p_inicial=0.1`, `schedule='Cosine'` | $p(t) = p_0 \cdot \frac{1+\cos(\pi t/T)}{2}$ |
| M09 | Otimizador Adam | `framework_investigativo_completo.py:treinar_modelo():L450-455` | `lr=0.01`, `beta1=0.9`, `beta2=0.999`, `epsilon=1e-8` | Parâmetros $\theta$ atualizados |
| M10 | Função de Custo (Cross-Entropy) | `framework_investigativo_completo.py:calcular_custo():L480-490` | `y_true`, `y_pred` | $L = -\sum y \log(\hat{y}) + (1-y)\log(1-\hat{y})$ |
| M11 | Critério de Early Stopping | `framework_investigativo_completo.py:treinar_modelo():L520-530` | `patience=20`, `min_delta=1e-4` | Flag de convergência |
| M12 | Dataset Moons | `framework_investigativo_completo.py:carregar_dataset():L150-165` | `n_samples=200`, `noise=0.1`, `random_state=42` | $(X, y)$ com shape $(200, 2)$ |
| M13 | Dataset Circles | `framework_investigativo_completo.py:carregar_dataset():L167-182` | `n_samples=200`, `noise=0.1`, `factor=0.5` | $(X, y)$ com shape $(200, 2)$ |
| M14 | Dataset Blobs | `framework_investigativo_completo.py:carregar_dataset():L184-199` | `n_samples=200`, `centers=2`, `cluster_std=0.5` | $(X, y)$ com shape $(200, 2)$ |
| M15 | Dataset Iris | `framework_investigativo_completo.py:carregar_dataset():L201-216` | `n_samples=150`, `n_features=4`, `n_classes=3` | $(X, y)$ com shape $(150, 4)$ |
| M16 | Normalização de Features | `framework_investigativo_completo.py:preprocessar_dados():L230-238` | Método: `StandardScaler` | $X_{norm} = \frac{X - \mu}{\sigma}$ |
| M17 | Split Train/Val/Test | `framework_investigativo_completo.py:dividir_dados():L240-250` | `train_ratio=0.6`, `val_ratio=0.2`, `test_ratio=0.2`, `random_state=42` | 3 subsets |
| M18 | Métrica: Acurácia | `framework_investigativo_completo.py:calcular_metricas():L600-610` | `y_true`, `y_pred` | $\text{Acc} = \frac{TP+TN}{TP+TN+FP+FN}$ |
| M19 | Métrica: F1-Score | `framework_investigativo_completo.py:calcular_metricas():L612-622` | `y_true`, `y_pred`, `average='weighted'` | $F1 = 2 \cdot \frac{P \cdot R}{P+R}$ |
| M20 | Métrica: Precision | `framework_investigativo_completo.py:calcular_metricas():L624-634` | `y_true`, `y_pred` | $P = \frac{TP}{TP+FP}$ |
| M21 | Métrica: Recall | `framework_investigativo_completo.py:calcular_metricas():L636-646` | `y_true`, `y_pred` | $R = \frac{TP}{TP+FN}$ |
| M22 | Teste Estatístico: ANOVA | `framework_investigativo_completo.py:analise_estatistica():L700-720` | `groups`, `alpha=0.05` | $F$-statistic, $p$-value |
| M23 | Teste Post-Hoc: Tukey HSD | `framework_investigativo_completo.py:analise_estatistica():L722-740` | `groups`, `alpha=0.05` | Comparações pareadas com $p$-values |
| M24 | Correção de Bonferroni | `framework_investigativo_completo.py:corrigir_pvalues():L750-760` | `p_values`, `n_comparisons` | $p_{corrigido} = \min(p \cdot n, 1)$ |
| M25 | Tamanho de Efeito: Cohen's d | `framework_investigativo_completo.py:calcular_effect_size():L780-795` | `group1`, `group2` | $d = \frac{\mu_1 - \mu_2}{\sqrt{(\sigma_1^2 + \sigma_2^2)/2}}$ |
| M26 | Intervalo de Confiança 95% | `framework_investigativo_completo.py:calcular_ic95():L800-815` | `data`, `confidence=0.95` | $(CI_{low}, CI_{high})$ |
| M27 | Logging de Experimentos | `framework_investigativo_completo.py:salvar_resultados():L850-870` | `results_dict`, `output_path='resultados/'` | `resultados_experimento.json` |
| M28 | Checkpoint de Modelo | `framework_investigativo_completo.py:salvar_checkpoint():L880-900` | `theta`, `epoch`, `loss`, `path='checkpoints/'` | `modelo_epoca_{epoch}.pkl` |
| M29 | Seed Aleatória Global | `framework_investigativo_completo.py:definir_seed():L50-65` | `seed=42` | Seeds fixas para `numpy`, `random`, `torch` |
| M30 | Geração de Figuras | `framework_investigativo_completo.py:gerar_visualizacoes():L950-1020` | `results_df`, `output_path='figuras/'` | `.png` files |

## Mapa de Dependências

### Bibliotecas Core

| Biblioteca | Versão | Uso no Método | Arquivo |
|------------|--------|---------------|---------|
| PennyLane | 0.38.0 | Construção de circuitos quânticos, operadores | `framework_investigativo_completo.py` |
| NumPy | 1.24.3 | Operações matriciais, álgebra linear | Todos |
| scikit-learn | 1.3.0 | Datasets, métricas, normalização | `framework_investigativo_completo.py:L150-250` |
| SciPy | 1.11.1 | Testes estatísticos (ANOVA, t-test) | `framework_investigativo_completo.py:L700-800` |
| Matplotlib | 3.7.2 | Visualizações | `framework_investigativo_completo.py:L950-1020` |
| Pandas | 2.0.3 | Manipulação de resultados tabulares | `framework_investigativo_completo.py:L850-900` |

### Hardware e Simuladores

| Componente | Configuração | Local |
|------------|--------------|-------|
| Simulador PennyLane | `default.qubit` | `framework_investigativo_completo.py:L240` |
| Número de Qubits | 4 | Configurável via CLI |
| Shots | `None` (statevector) | `framework_investigativo_completo.py:L245` |

## Rastreabilidade de Resultados

### Arquivo de Resultados → Seção do Artigo

| Arquivo de Resultado | Métrica/Figura | Seção do Artigo | Linha do Código |
|----------------------|----------------|-----------------|-----------------|
| `resultados_experimento.json` | Tabela de acurácias | Results, Tabela 1 | `framework_investigativo_completo.py:L850` |
| `figura_comparacao_ansatze.png` | Figura 2: Comparação de ansätze | Results, Fig. 2 | `framework_investigativo_completo.py:L980` |
| `figura_ruido_benefico.png` | Figura 3: Regime benéfico | Results, Fig. 3 | `framework_investigativo_completo.py:L1000` |
| `anova_results.csv` | Teste ANOVA | Results, texto | `framework_investigativo_completo.py:L720` |
| `tabela_s1_configuracoes.csv` | Tabela S1: Todas as configs | Supplementary | `framework_investigativo_completo.py:L860` |

## Configurações de Execução

### Comandos de Execução Documentados

```bash
# Execução completa (todas as configurações)
python framework_investigativo_completo.py \
    --datasets moons circles blobs iris \
    --ansatze BasicEntangler StronglyEntangling RandomEntangling \
    --noise-types Depolarizing AmplitudeDamping PhaseDamping \
    --noise-params 0.0 0.01 0.05 0.1 \
    --schedules Constant Linear Cosine \
    --seeds 42 123 456 789 1024 \
    --output resultados_completos/

# Execução rápida (subset para teste)
python framework_investigativo_completo.py \
    --datasets moons \
    --ansatze BasicEntangler \
    --noise-types Depolarizing \
    --noise-params 0.0 0.05 \
    --schedules Constant \
    --seeds 42 \
    --output resultados_teste/
```

### Seeds Documentadas

| Seed | Propósito | Localização |
|------|-----------|-------------|
| 42 | Dataset split, inicialização | Global |
| 123 | Replicação 2 | Loop de experimentos |
| 456 | Replicação 3 | Loop de experimentos |
| 789 | Replicação 4 | Loop de experimentos |
| 1024 | Replicação 5 | Loop de experimentos |

## Verificação de Consistência

### Checklist

- [ ] Todos os componentes da metodologia têm entrada nesta tabela
- [ ] Todos os caminhos de arquivo estão corretos e verificados
- [ ] Todos os parâmetros correspondem aos valores no código
- [ ] Todas as equações correspondem às implementações
- [ ] Todos os artefatos listados existem nos diretórios de output

### Script de Verificação

```python
# verificar_rastreabilidade.py
import os
import re

def verificar_mapeamento(tabela_csv, codigo_path):
    """Verifica se todos os mappings estão corretos."""
    with open(tabela_csv) as f:
        for linha in f:
            arquivo, funcao, linha_num = extrair_localizacao(linha)
            if not verificar_existe(codigo_path, arquivo, funcao, linha_num):
                print(f"❌ ERRO: {arquivo}:{funcao}:{linha_num} não encontrado")
    print("✅ Verificação completa")

# Execute:
# python verificar_rastreabilidade.py tabela_codigo_metodo.csv framework_investigativo_completo.py
```

---

**Última atualização**: [Data]  
**Revisor**: [Nome]  
**Status**: [Draft | Verificado | Aprovado]
