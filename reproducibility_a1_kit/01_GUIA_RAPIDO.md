# 01 — Guia Rápido (banca)

## Pré-requisitos

- Python 3.9+
- Dependências instaladas em `vqc_drug_v10a1`

## Execução em 3 comandos

```bash
cd vqc_drug_v10a1
pip install -e .
../reproducibility_a1_kit/scripts/run_reproducibility.sh EGFR 100 42
```

Parâmetros:
- `target`: `EGFR|HIV|Malaria|COVID`
- `trials`: número de trials Optuna
- `seed`: seed determinística

## Validação automática

```bash
python ../reproducibility_a1_kit/scripts/validate_artifacts.py --results-dir <DIRETORIO_RESULTADOS>
```

## Evidências mínimas para banca

- `01_protocolo_pre_registrado_*.json`
- `environment_snapshot.txt`
- `cv_comparison.csv`
- `02_validation_report.json`
- `02_validation_report.md`
- `checksums_final.sha256`

Se os arquivos existirem e os gates forem aprovados, o experimento está metodologicamente validado para apresentação.
