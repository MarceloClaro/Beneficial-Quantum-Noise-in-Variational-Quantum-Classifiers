#!/usr/bin/env bash
set -euo pipefail

TARGET="${1:-EGFR}"
TRIALS="${2:-100}"
SEED="${3:-42}"
PERMUTATIONS="${4:-2000}"

if [[ ! -d "vqc_drug_v10a1" ]]; then
  echo "[ERRO] Execute este script a partir da raiz do repositório."
  exit 1
fi

cd vqc_drug_v10a1

echo "[INFO] Executando pipeline: target=${TARGET}, trials=${TRIALS}, seed=${SEED}, permutations=${PERMUTATIONS}"
vqc-drug-a1 \
  --target "${TARGET}" \
  --trials "${TRIALS}" \
  --seed "${SEED}" \
  --permutations "${PERMUTATIONS}" \
  --alpha 0.05

echo "[INFO] Execução concluída. Use validate_artifacts.py no diretório de resultados gerado."
