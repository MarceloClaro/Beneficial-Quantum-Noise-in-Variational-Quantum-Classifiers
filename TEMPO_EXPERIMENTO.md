# ⏱️ Registro de Tempo e Previsão do Experimento.

## 📊 Informações Gerais

- **Versão do Framework**: v7.2
- **Modo de Execução**: Completo (Grid Search + Otimização Bayesiana)
- **Comando**: `python framework_investigativo_completo.py --run-both`
- **Pasta de Resultados**: `resultados_2025-10-28_22-35-28/`

---

## 🚀 Início da Execução

- **Data/Hora de Início**: 28 de outubro de 2025, 22:35:28
- **Datasets Carregados**: 5 (moons, circles, iris, breast_cancer, wine)
  - moons: 280 treino, 120 teste
  - circles: 280 treino, 120 teste
  - iris: 70 treino, 30 teste
  - breast_cancer: 398 treino, 171 teste
  - wine: 91 treino, 39 teste

---

## 📈 Status Atual

**Última Atualização**: 29 de outubro de 2025, 05:18:32

### Progresso
- **Experimentos Completos**: 142 de 8280
- **Percentual Concluído**: 1.71%
- **Experimentos Restantes**: 8138

### Experimento Atual
- **ID**: [143/8280]
- **Dataset**: moons
- **Seed**: 44
- **Qubits**: 4
- **Camadas**: 2
- **Arquitetura**: basic_entangler
- **Inicialização**: matematico
- **Tipo de Ruído**: crosstalk
- **Nível de Ruído**: 0.0000

---

## ⏰ Análise de Tempo

### Tempo Decorrido
- **Início**: 28/out/2025 22:35:28
- **Última Atualização**: 29/out/2025 05:18:32
- **Tempo Decorrido**: ~6 horas e 43 minutos

### Tempos por Experimento (Amostra dos últimos experimentos)
- Exp 120: 148.8s (~2.5 min)
- Exp 121: 187.1s (~3.1 min)
- Exp 126: 186.1s (~3.1 min)
- Exp 131: 199.1s (~3.3 min)
- Exp 136: 189.5s (~3.2 min)
- Exp 140: 74.8s (~1.2 min) - convergência rápida
- **Exp 141: 372.5s (~6.2 min)** - experimento com crosstalk (mais pesado)
- **Exp 142: 343.0s (~5.7 min)** - experimento com crosstalk

### Estatísticas de Tempo
- **Tempo Médio por Experimento**: ~160 segundos (~2.7 minutos)
- **Tempo Mínimo Observado**: 74.0s
- **Tempo Máximo Observado**: 372.5s
- **Nota**: Experimentos com ruído crosstalk são significativamente mais lentos (300-375s vs 75-200s)

---

## 🔮 Previsão de Conclusão

### Cenários de Estimativa

#### Cenário Conservador (Média: 160s por experimento)
- **Tempo Total Estimado**: 8280 × 160s = 1,324,800 segundos
- **Conversão**: 368 horas = **15.3 dias**
- **Previsão de Conclusão**: ~13 de novembro de 2025

#### Cenário Otimista (Média: 120s por experimento)
- **Tempo Total Estimado**: 8280 × 120s = 993,600 segundos
- **Conversão**: 276 horas = **11.5 dias**
- **Previsão de Conclusão**: ~10 de novembro de 2025

#### Cenário Pessimista (Média: 200s por experimento)
- **Tempo Total Estimado**: 8280 × 200s = 1,656,000 segundos
- **Conversão**: 460 horas = **19.2 dias**
- **Previsão de Conclusão**: ~17 de novembro de 2025

### Marcos Esperados

Com base no tempo médio de 160s por experimento:

| Marco | Experimentos | % | Data/Hora Prevista | Dias Restantes |
|-------|-------------|---|-------------------|----------------|
| 25% | 2,070 | 25% | 02/nov/2025 ~07:00 | 3.8 dias |
| 50% | 4,140 | 50% | 06/nov/2025 ~15:30 | 7.6 dias |
| 75% | 6,210 | 75% | 11/nov/2025 ~00:00 | 11.4 dias |
| 100% | 8,280 | 100% | 13/nov/2025 ~08:30 | 15.2 dias |

---

## 📁 Artefatos Gerados (até agora)

### CSVs Individuais
- **Quantidade**: 142 arquivos
- **Localização**: `resultados_2025-10-28_22-35-28/experimentos_individuais/`
- **Formato**: exp_00001.csv até exp_00142.csv

### Visualizações: 
- **Circuitos PNG**: 142 arquivos (um por experimento)
- **Barren Plateau 3D**: 142 arquivos (um por experimento)
- **Localização**: `resultados_2025-10-28_22-35-28/circuitos/` e `barren_plateaus/`

---

## 🎯 Próximas Etapas Automáticas (ao finalizar)

1. **Consolidação de CSVs** → `resultados_completos_artigo.csv`
2. **Comparação de Baselines** → `comparacao_baselines.csv`
3. **Geração de Metadados** → `metadata_orchestrator.json`
4. **Análises Estatísticas Completas**
   - ANOVA, Kruskal-Wallis
   - Post-hoc: Bonferroni, Scheffé
   - Effect sizes: Cohen's d, Glass's Δ, Hedges' g
5. **Visualizações Interativas** (9 figuras HTML)
   - figura2: Beneficial Noise
   - figura2b: Beneficial Noise IC95
   - figura3: Noise Types
   - figura3b: Noise Types IC95
   - figura4: Initialization
   - figura5: Architecture Tradeoffs
   - figura6: Effect Sizes
   - figura7: Overfitting
   - Análises Profundas (PCA, Clustering, Correlação, Sensibilidade)
6. **Resumo Final** com melhores configurações por dataset

---

## 💡 Observações Técnicas

### Desempenho
- Experimentos com ruído **crosstalk** são 2-5× mais lentos que outros tipos
- Experimentos que convergem rapidamente: 74-80s
- Experimentos com convergência lenta: 180-200s
- Experimentos com crosstalk: 300-375s

### Recursos do Sistema
- **Python**: 3.9.23
- **Device**: PennyLane default.mixed (simulação)
- **CPU**: Processo ativo e saudável
- **Memória**: Estável (sem vazamentos observados)

### Monitoramento
- Monitor automático em `tools/monitor_progress.py`
- Verificação a cada 5 minutos
- Notificações em marcos de 25%, 50%, 75% e 100%

---

## 📝 Notas

- A execução está **estável e progredindo normalmente**
- Variações no tempo são esperadas devido a:
  - Diferentes tipos de ruído (crosstalk é mais pesado)
  - Velocidade de convergência do otimizador
  - Complexidade da superfície de loss
- O framework salva resultados incrementalmente (seguro contra interrupções)
- Todos os logs estão sendo registrados em tempo real

---

**Última Atualização Automática**: 29 de outubro de 2025, 05:18:32

*Este documento será atualizado periodicamente conforme o experimento progride.*
