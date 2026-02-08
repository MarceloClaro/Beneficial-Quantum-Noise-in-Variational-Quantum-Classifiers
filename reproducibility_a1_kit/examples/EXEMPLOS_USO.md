# Exemplos de uso

## 1) Execução rápida para banca (EGFR)

```bash
./reproducibility_a1_kit/scripts/run_reproducibility.sh EGFR 100 42 2000
```

## 2) Validação da pasta de resultados

```bash
python reproducibility_a1_kit/scripts/validate_artifacts.py \
  --results-dir vqc_drug_v10a1/results_EGFR_YYYY-MM-DD_hh-mm-ss
```

## 3) Execução com dataset HIV

```bash
./reproducibility_a1_kit/scripts/run_reproducibility.sh HIV 200 123 5000
```
