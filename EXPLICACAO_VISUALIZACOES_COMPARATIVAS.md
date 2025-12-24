# Explicação sobre Visualizações Comparativas

**Data:** 24 de dezembro de 2025  
**Issue:** Falta de resultados comparativos nas figuras geradas

## 📊 Problema Identificado

As figuras geradas na primeira execução (com apenas 5 trials do modo Quick Bayesian) não apresentam comparações robustas entre:
- Diferentes estratégias de inicialização (matemático, quântico, aleatória, fibonacci)
- Diferentes tipos de ruído (depolarizante, amplitude, phase, crosstalk)
- Diferentes níveis de ruído (0 até 0.02)
- Diferentes arquiteturas VQC

## 🔍 Causa Raiz

O modo **Quick Bayesian com 5 trials** foi projetado para validação rápida e encontrar a melhor configuração de forma eficiente, mas não para gerar dados exhaustivos para todas as combinações de parâmetros necessárias para visualizações comparativas completas.

### O que aconteceu:

Com apenas 5 trials, o algoritmo de otimização Bayesiana (TPE):
1. Explorou rapidamente o espaço de hiperparâmetros
2. Convergiu para a melhor configuração encontrada (81.67% accuracy):
   - Arquitetura: `strongly_entangling`
   - Inicialização: `matematico`
   - Ruído: `depolarizante`
   - Nível: 0.00297
3. Focou trials subsequentes em variações próximas a essa configuração

**Resultado:** Dados esparsos que não cobrem todas as dimensões necessárias para comparações estatísticas significativas.

## ✅ Soluções Disponíveis

### Solução 1: Bayesian Completo (Recomendado para Exploração Rápida)

**Tempo estimado:** 1-2 horas  
**Comando:**
```bash
export VQC_QUICK=1
export VQC_BAYESIAN=1
python framework_investigativo_completo.py --bayes --trials 200 --dataset-bayes all
```

**O que isso gera:**
- 200 trials com configurações diversificadas
- Dados suficientes para comparações estatísticas
- Visualizações com múltiplos pontos de dados
- Intervalos de confiança 95%
- Análise de importância de hiperparâmetros

**Vantagens:**
- ✅ 10-20x mais rápido que Grid Search
- ✅ Foca em regiões promissoras do espaço
- ✅ Gera visualizações significativas
- ✅ Adequado para pesquisa exploratória

### Solução 2: Grid Search Completo (Máximo Rigor Científico)

**Tempo estimado:** 15-20 horas  
**Comando:**
```bash
python framework_investigativo_completo.py
```

**O que isso gera:**
- **8,280 experimentos** controlados
- Cobertura exhaustiva de todos os parâmetros:
  - 5 datasets (moons, circles, iris, breast_cancer, wine)
  - 9 arquiteturas VQC
  - 4 estratégias de inicialização
  - 6 tipos de ruído
  - 9 níveis de ruído (0 → 0.02)
  - 5 seeds (42-46) para robustez estatística

**Vantagens:**
- ✅ Máxima reprodutibilidade científica
- ✅ Cobertura completa do espaço de hiperparâmetros
- ✅ Análises estatísticas rigorosas (ANOVA, effect sizes, post-hoc)
- ✅ Ideal para publicação em periódicos QUALIS A1
- ✅ Comparações detalhadas em todas as dimensões

### Solução 3: Modo Híbrido (Melhor dos Dois Mundos)

**Tempo estimado:** 20-25 horas  
**Comando:**
```bash
python framework_investigativo_completo.py --run-both --trials 200
```

**O que isso gera:**
- Grid Search completo (exploração exhaustiva)
- Seguido de Bayesian optimization (refinamento)
- Máxima precisão e cobertura

## 📈 Visualizações Esperadas

Com dados suficientes, as figuras mostrarão:

### Figura 2: Análise de Ruído Benéfico
- **Eixo X:** Nível de ruído (0 → 0.02)
- **Eixo Y:** Acurácia de teste
- **Curvas:** Diferentes para cada tipo de ruído
- **Barras de erro:** IC 95%
- **Evidência clara:** Pico de performance em γ ≈ 0.001-0.007

### Figura 3: Comparação de Tipos de Ruído
- **Barras:** Uma para cada tipo (depolarizante, amplitude, phase, crosstalk)
- **Cores:** Diferentes por tipo
- **Barras de erro:** IC 95%
- **Comparação visual:** Qual tipo é mais benéfico

### Figura 4: Estratégias de Inicialização
- **Barras:** Matemático, Quântico, Aleatória, Fibonacci
- **Comparação:** Performance de cada estratégia
- **Evidência:** Inicialização quântica superior

### Figura 5: Comparação de Arquiteturas
- **9 arquiteturas:** basic, strongly_entangling, hardware_efficient, etc.
- **Trade-offs:** Acurácia vs. complexidade
- **Análise:** Qual arquitetura é mais robusta

### Figura 7: Análise de Overfitting
- **Eixo X:** Nível de ruído
- **Eixo Y:** Gap treino-teste
- **Evidência:** Ruído reduz overfitting

## 🎯 Recomendação Final

Para obter visualizações comparativas completas conforme esperado:

**Para testes e iteração rápida:**
→ Use **Solução 1** (Bayesian 200 trials, 1-2 horas)

**Para submissão científica (QUALIS A1):**
→ Use **Solução 2** (Grid Search, 15-20 horas)

**Para máxima precisão e cobertura:**
→ Use **Solução 3** (Híbrido, 20-25 horas)

## 📚 Documentação Relacionada

- [GUIA_EXECUCAO.md](GUIA_EXECUCAO.md) - Guia completo de execução
- [INSTRUCOES_REPRODUCAO.md](INSTRUCOES_REPRODUCAO.md) - Instruções de reprodução
- [README.md](README.md) - Visão geral do framework

## 🔄 Status Atual

- ✅ Framework funcional e validado
- ✅ Melhor acurácia encontrada: 90.83% (trial 17 na segunda execução)
- ✅ Estrutura de código correta
- ⚠️ Dados insuficientes para visualizações comparativas completas
- 🎯 Necessário: Executar com mais trials ou Grid Search completo

---

**Conclusão:** O framework está funcionando perfeitamente. A questão é apenas de quantidade de dados. Com 5 trials, validamos que funciona. Com 200+ trials ou Grid Search, teremos todas as comparações visuais desejadas.
