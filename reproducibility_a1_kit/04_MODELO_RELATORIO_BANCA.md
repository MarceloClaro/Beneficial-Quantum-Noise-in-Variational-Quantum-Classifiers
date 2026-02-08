# 04 — Modelo de Relatório Técnico para Banca

## 1. Identificação do experimento

- Projeto:
- Data/hora da execução:
- Dataset alvo:
- Trials:
- Seed:
- Ambiente (SO, Python):

## 2. Objetivo

Descrever a hipótese principal e o endpoint primário (AUC).

## 3. Metodologia resumida

- preparação dos dados;
- otimização com validação cruzada;
- validação inferencial (teste unilateral, bootstrap e permutação);
- critérios objetivos de aceite.

## 4. Resultados

Preencher a partir de `02_validation_report.json`:

- AUC média:
- IC95%:
- p-valor unilateral:
- p-valor permutação:
- tamanho de efeito:
- estabilidade (CV):

## 5. Gates de qualidade

- gate_pvalue:
- gate_permutation:
- gate_ci_above_random:
- gate_stability_cv_lt_10pct:

## 6. Conclusão para banca

- Situação final: Aprovado / Ajustes / Reprovado.
- Riscos ou limitações observadas.
- Próximos passos para submissão Qualis A1.
