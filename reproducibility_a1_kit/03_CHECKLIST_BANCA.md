# 03 — Checklist de Validação (Banca CNPq)

Marque cada item com ✅ / ⚠️ / ❌.

## A. Ambiente e execução

- [ ] Versão de Python registrada.
- [ ] Pacote `vqc_drug_v10a1` instalado.
- [ ] Execução com `--seed` e `--permutations` explícitos.

## B. Artefatos obrigatórios

- [ ] `01_protocolo_pre_registrado_*.json` presente.
- [ ] `environment_snapshot.txt` presente.
- [ ] `cv_comparison.csv` presente.
- [ ] `02_validation_report.json` presente.
- [ ] `02_validation_report.md` presente.
- [ ] `checksums_final.sha256` presente.

## C. Validação estatística

- [ ] `gate_pvalue == true`.
- [ ] `gate_permutation == true`.
- [ ] `gate_ci_above_random == true`.
- [ ] `gate_stability_cv_lt_10pct == true`.

## D. Auditoria

- [ ] Pasta de resultados arquivada sem modificação.
- [ ] Relatório técnico preenchido.
- [ ] Ata de validação institucional anexada.
