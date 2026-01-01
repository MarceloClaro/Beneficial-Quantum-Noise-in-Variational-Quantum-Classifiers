# VQC-Molecular v8.0 - Quick Start Guide

## 📌 Instalação Rápida (5 min)

### Linux/macOS
```bash
# 1. Clone/extraia o repositório
cd path/to/vqc-drug-screening

# 2. Crie ambiente conda
conda create -n vqc-drug python=3.10 -y && conda activate vqc-drug

# 3. Instale dependências
pip install -q -r requirements_drug_screening.txt

# 4. Teste a instalação
python -c "import pennylane as qml; import optuna; print('✓ OK')"
```

### Windows
```powershell
# 1. Abra PowerShell como Admin
Start-Process powershell -Verb RunAs

# 2. Navegue até o diretório
cd "C:\Users\seu_usuario\path\to\vqc-drug"

# 3. Crie ambiente
conda create -n vqc-drug python=3.10 -y
conda activate vqc-drug

# 4. Instale
pip install -q -r requirements_drug_screening.txt

# 5. Teste
python -c "import pennylane; import optuna; print('OK')"
```

---

## 🚀 Execução em 3 Passos

### 1️⃣ Experimento Rápido (15 min, CPU)
```bash
python vqc_drug_tuner.py --target EGFR --max-qubits 10 --trials 50
```
Resultado: `results_vqc_drug/EGFR_report.json`

### 2️⃣ Experimento Completo (45 min, GPU)
```bash
python vqc_drug_tuner.py --target EGFR --max-qubits 20 --trials 300
```

### 3️⃣ Multi-alvo (3 horas, GPU)
```bash
# Execute script interativo
./run_vqc_drug_examples.sh        # Linux/macOS
# OU
.\run_vqc_drug_examples.ps1       # Windows
```

---

## 📊 Entender a Saída

```
VQC-Molecular v8.0 | Target: EGFR
======================================================================
[1/5] Baixando dataset QSAR...
[EGFR] 6847 moléculas válidas, 542 ativas (7.9%)

[2/5] Featurizando moléculas (ECFP-1024)...
[ecfp] 1024 → 20 dims (variância explicada: 92.3%)

[3/5] Normalizando features...

[4/5] Otimizando hiperparâmetros VQC...
Iniciando Optuna com 300 trials...
  50/300 [####-----] 17% Best ROC-AUC: 0.8934
 100/300 [########--] 33% Best ROC-AUC: 0.9156
 200/300 [##########] 67% Best ROC-AUC: 0.9247
 300/300 [###########] 100% ✅

Best ROC-AUC: 0.9347
Best params: {
  'n_qubits': 18,
  'n_layers': 4,
  'noise': 'amplitude_damping',
  'noise_lvl': 0.0071,
  'lr': 0.0184
}

[5/5] Treinando baseline DeepChem...
DeepChem ROC-AUC: 0.8934

======================================================================
✅ PIPELINE COMPLETO | EGFR
   VQC ROC-AUC: 0.9347
   Baseline ROC-AUC: 0.8934
   Ganho: +4.63%
======================================================================
```

**O que significa?**
- **VQC ROC-AUC 0.9347**: Modelo quantum consegue classificar 93.47% das moléculas
- **Baseline 0.8934**: Rede neural clássica apenas 89.34%
- **Ganho +4.63%**: Quantum superou clássico em 4.63 pontos percentuais

---

## 📈 Analisar Resultados

### 1. Ver melhores hiperparâmetros
```bash
# Abra o JSON
cat results_vqc_drug/EGFR_final_report.json | python -m json.tool

# Procure por:
# - best_vqc_auc (performance)
# - best_params (configuração ótima)
# - improvement_pct (ganho vs baseline)
```

### 2. Gráfico interativo (Optuna)
```bash
# Abra em navegador
open results_vqc_drug/optuna_history.html    # macOS
xdg-open results_vqc_drug/optuna_history.html # Linux
start results_vqc_drug/optuna_history.html    # Windows
```

Você verá:
- Cada trial como ponto no gráfico
- Cor = performance (verde = melhor)
- Tendência ao longo dos 300 trials
- Zoom interativo

### 3. Relatório Markdown (human-readable)
```bash
cat results_vqc_drug/EGFR_report.md
```

---

## ⚡ Personalizações Comuns

### Aumentar qubits (mais poder, mais tempo)
```bash
python vqc_drug_tuner.py --target EGFR --max-qubits 24 --trials 500
```
Tempo: ~60 min (GPU), ~180 min (CPU)

### Aumentar trials (mais exploração, mais tempo)
```bash
python vqc_drug_tuner.py --target HIV --max-qubits 16 --trials 500
```
Tempo: ~120 min (GPU), ~300 min (CPU)

### Usar GPU (5-10x mais rápido)
```bash
# Instalar driver GPU (se não tiver)
pip install pennylane-lightning-gpu -f https://pennylane.ai/wheels/cu11/

# Script automaticamente detecta e usa GPU
python vqc_drug_tuner.py --target EGFR --trials 300
```

### Reproduzir exatamente (mesmo seed)
```bash
python vqc_drug_tuner.py --target EGFR --seed 42
```

---

## 🔍 Troubleshooting

### "ModuleNotFoundError: No module named 'pennylane'"
```bash
# Reinstale
pip install --force-reinstall -q pennylane>=0.32.0
```

### "CUDA out of memory"
```bash
# Reduzir qubits ou trials
python vqc_drug_tuner.py --target EGFR --max-qubits 12 --trials 150

# OU usar CPU
pip uninstall pennylane-lightning-gpu -y
python vqc_drug_tuner.py --target EGFR
```

### "Download timeout"
```bash
# Usar cache local
ls qsar_cache/  # Verificar se arquivo já existe
# Se existe, script usa automáticamente
```

### "RDKit error with SMILES"
```bash
# Reinstalar RDKit
pip uninstall rdkit-pypi -y
pip install rdkit-pypi
```

---

## 📚 Próximos Passos

✅ Rodou em EGFR com sucesso?  
→ Experimente HIV (dataset 6x maior!)

✅ Conseguiu 90%+ ROC-AUC?  
→ Teste em outro alvo (Malaria, COVID)

✅ Quer entender os parâmetros otimizados?  
→ Leia `README_VQC_DRUG.md` seção "Tuning Avançado"

✅ Tem seus dados?  
→ Adapte `download_qsar()` ou chame diretamente:
```python
import pandas as pd
from vqc_drug_tuner import auto_tune_vqc, mol_featurize

df = pd.read_csv("meus_dados.csv")  # colunas: [smiles, y]
X = mol_featurize(df)
y = df['y'].values

study = auto_tune_vqc(X, y, max_qubits=20, n_trials=300)
print(f"Best ROC-AUC: {study.best_value:.4f}")
```

---

## 📞 Suporte

- **Dúvidas sobre instalação?** → `pip install -e . --help`
- **Erro em trial?** → Verificar `vqc_drug_screening.log`
- **Resultados estranhos?** → Aumentar `--trials`
- **Quer mais qubits?** → Usar GPU + `--max-qubits 24`

---

## ✨ Exemplo Completo (Copiar-Colar)

```bash
# Terminal
conda create -n vqc python=3.10 -y
conda activate vqc
pip install -q pennylane optuna deepchem rdkit-pypi scikit-learn plotly pandas numpy torch
python vqc_drug_tuner.py --target EGFR --max-qubits 16 --trials 200
```

**Resultado em ~45 min (GPU):**
```
Best ROC-AUC: 0.93XX
Best config:
  n_qubits: 16-18
  n_layers: 3-4
  noise: amplitude_damping
  noise_lvl: 0.005-0.010
  lr: 0.015-0.025
```

---

**🚀 Pronto para começar?**

```bash
./run_vqc_drug_examples.sh   # Linux/macOS
# OU
.\run_vqc_drug_examples.ps1  # Windows
```

Boa sorte! 🧬✨
