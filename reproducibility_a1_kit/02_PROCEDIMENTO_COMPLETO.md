# 02 — Procedimento Completo de Reprodutibilidade

## 1. Preparação do ambiente

1. Registrar versão de Python e sistema operacional.
2. Instalar pacote em modo editável.
3. Confirmar que o comando `vqc-drug-a1` está disponível.

## 2. Execução padronizada

Executar com parâmetros explícitos de validação:

```bash
vqc-drug-a1 \
  --target EGFR \
  --trials 500 \
  --seed 42 \
  --permutations 5000 \
  --alpha 0.05
```

## 3. Artefatos mandatórios

A execução deve gerar:
- pré-registro (`01_protocolo_pre_registrado_*.json`)
- snapshot de ambiente
- comparação de trials (`cv_comparison.csv`)
- validação estatística (`02_validation_report.json` + `.md`)
- checksums finais

## 4. Critérios de aceite técnico

No `02_validation_report.json`, verificar:
- `quality_gates.gate_pvalue == true`
- `quality_gates.gate_permutation == true`
- `quality_gates.gate_ci_above_random == true`
- `quality_gates.gate_stability_cv_lt_10pct == true`

## 5. Auditoria e rastreabilidade

1. Conferir existência de `checksums_final.sha256`.
2. Arquivar a pasta de resultados sem alterações.
3. Anexar checklist e ata de validação institucional.

## 6. Publicação (Qualis A1)

Recomendação para seção metodológica:
- explicitar hipótese primária vs baseline aleatório;
- relatar IC bootstrap, teste de permutação e tamanho de efeito;
- descrever seed, versão de ambiente e trilha de checksums.
