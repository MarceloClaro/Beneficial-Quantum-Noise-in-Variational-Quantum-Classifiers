# Relatório de Conivência Código-Texto

**Código:** framework_investigativo_completo.py
**Artigo:** artigo_cientifico

## 🎯 Congruência Geral

**76.8%** (⚠️ REGULAR - Algumas inconsistências)

---

## 💻 Análise do Código

- **Linhas de Código:** 5661
- **Classes:** 24
- **Funções:** 95
- **Datasets:** 16
- **Modelos de Ruído:** 6
- **Ansätze:** 2
- **Métricas:** 3
- **Seeds:** [42]

---

## 📊 Verificação Componente por Componente

| Componente | Código | Artigo | Congruência | Status |
|------------|--------|--------|-------------|--------|
| Classes Implementadas | 24 classes | 11 mencionadas | 45.8% | ❌ |
| Datasets | 16: acuracia_teste, arquitetura, breast_cancer, circles, dataset, gap_abs, iris, markers, moons, nivel_ruido, nivel_ruido_cat, tipo_ruido, wine, x_test, x_train, y_test | 6: breast, circles, digits, iris, moons, wine | 25.0% | ❌ |
| Modelos de Ruído | 6: amplitude, bitflip, damping, depolarizing, phase, phaseflip | 6: amplitude, bitflip, damping, depolarizing, phase, phaseflip | 100.0% | ✅ |
| Ansätze Quânticos | 2 tipos | 9 mencionados | 100.0% | ✅ |
| Métricas de Avaliação | 3: accuracy, confusion_matrix, f1 | 6: accuracy, auc, f1, precision, recall, roc | 66.7% | ❌ |
| Bibliotecas Principais | 2: numpy, sklearn | 2: numpy, sklearn | 100.0% | ✅ |
| Seeds de Reprodutibilidade | [42] | [42] | 100.0% | ✅ |

---

## 📝 Detalhes das Verificações

### Classes Implementadas ❌

- **Código:** 24 classes
- **Artigo:** 11 mencionadas
- **Congruência:** 45.8%
- **Detalhes:** Classes: ConstantesFundamentais, ScheduleRuido, DetectorBarrenPlateau, MonitorEmaranhamento, OtimizadorAvancado

### Datasets ❌

- **Código:** 16: acuracia_teste, arquitetura, breast_cancer, circles, dataset, gap_abs, iris, markers, moons, nivel_ruido, nivel_ruido_cat, tipo_ruido, wine, x_test, x_train, y_test
- **Artigo:** 6: breast, circles, digits, iris, moons, wine
- **Congruência:** 25.0%
- **Detalhes:** Match: circles, iris, moons, wine

### Modelos de Ruído ✅

- **Código:** 6: amplitude, bitflip, damping, depolarizing, phase, phaseflip
- **Artigo:** 6: amplitude, bitflip, damping, depolarizing, phase, phaseflip
- **Congruência:** 100.0%
- **Detalhes:** Match: amplitude, bitflip, damping, depolarizing, phase, phaseflip

### Ansätze Quânticos ✅

- **Código:** 2 tipos
- **Artigo:** 9 mencionados
- **Congruência:** 100.0%
- **Detalhes:** Implementados: 2, Mencionados: 9

### Métricas de Avaliação ❌

- **Código:** 3: accuracy, confusion_matrix, f1
- **Artigo:** 6: accuracy, auc, f1, precision, recall, roc
- **Congruência:** 66.7%
- **Detalhes:** Match: accuracy, f1

### Bibliotecas Principais ✅

- **Código:** 2: numpy, sklearn
- **Artigo:** 2: numpy, sklearn
- **Congruência:** 100.0%
- **Detalhes:** Match: numpy, sklearn

### Seeds de Reprodutibilidade ✅

- **Código:** [42]
- **Artigo:** [42]
- **Congruência:** 100.0%
- **Detalhes:** Seeds críticos para reprodutibilidade

## 💡 Recomendações

### ⚠️ Inconsistências Detectadas

- **Classes Implementadas:** Revisar e alinhar código e texto
- **Datasets:** Revisar e alinhar código e texto
- **Métricas de Avaliação:** Revisar e alinhar código e texto

**Ação Recomendada:** Atualizar seções do artigo para refletir com precisão os componentes implementados no código.
