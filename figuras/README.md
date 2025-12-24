# Figuras QUALIS A1 - Framework de Ruído Quântico Benéfico

Este diretório contém as visualizações científicas geradas pelo framework, todas em conformidade com padrões QUALIS A1.

## Especificações Técnicas

- **Resolução**: 300 DPI (4800×3000 pixels exportados como 1600×1000)
- **Formatos**: PNG, PDF (vetorial), SVG (vetorial), HTML (interativo)
- **Fonte**: Times New Roman (padrão científico)
- **Bordas**: 2px espelhadas em todos os eixos
- **Intervalos de Confiança**: 95% CI nas análises estatísticas

## Figuras Disponíveis

### Figura 2: Beneficial Noise Analysis
**Arquivo**: `figura2_beneficial_noise.png`
- Análise do impacto do ruído quântico na acurácia
- Demonstra região ótima de ruído benéfico

### Figura 2b: Beneficial Noise with 95% CI
**Arquivo**: `figura2b_beneficial_noise_ic95.png`
- Acurácia média com intervalos de confiança de 95%
- Rigor estatístico máximo

### Figura 3: Noise Types Comparison
**Arquivo**: `figura3_noise_types.png`
- Comparação entre 5 tipos de ruído quântico

### Figura 3b: Noise Types with 95% CI
**Arquivo**: `figura3b_noise_types_ic95.png`
- Comparação estatística entre tipos de ruído

### Figura 4: Initialization Strategies
**Arquivo**: `figura4_initialization.png`
- Impacto de diferentes estratégias de inicialização

### Figura 5: Architecture Trade-offs
**Arquivo**: `figura5_architecture_tradeoffs.png`
- Análise de trade-offs entre arquiteturas VQC

### Figura 7: Overfitting Analysis
**Arquivo**: `figura7_overfitting.png`
- Análise do gap treino-teste

## Notas

As figuras são geradas automaticamente durante a execução do framework.
Para regenerá-las, execute:

```bash
python framework_investigativo_completo.py --bayes --trials 5 --dataset-bayes moons
```

Os arquivos serão criados em `resultados_YYYY-MM-DD_HH-MM-SS/`.
