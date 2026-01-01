# VQC-Molecular v10.0-A1 - Quick Installation Guide

## 🚀 1-Minute Setup

### Option 1: pip (Recommended)
```bash
cd vqc_drug_v10a1
pip install -e .

# Test installation
vqc-drug-a1 --help
```

### Option 2: Conda (Isolated Environment)
```bash
cd vqc_drug_v10a1
conda env create -f docker/environment.yml
conda activate vqc-a1
pip install -e .

# Test
vqc-drug-a1 --help
```

### Option 3: Docker (100% Reproducible)
```bash
cd vqc_drug_v10a1
docker build -t vqc-a1:latest docker/
docker run vqc-a1:latest vqc-drug-a1 --help
```

---

## ⚡ Quick Test (30 seconds)

```bash
# Run 10-trial smoke test
vqc-drug-a1 --target EGFR --trials 10 --out-dir test_results

# Check outputs
ls test_results/
```

Expected output:
```
test_results/
├── 01_protocolo_pre_registrado_*.json
├── cv_comparison.csv
├── fig1_auc_heatmap.png
├── checksums_final.sha256
└── ...
```

---

## 🔬 Run Full Pipeline

```bash
# EGFR with 500 trials (~30 min CPU, ~10 min GPU)
vqc-drug-a1 --target EGFR --max-qubits 18 --trials 500

# HIV with extended configuration
vqc-drug-a1 --target HIV --max-qubits 16 --trials 1000 --power 0.85
```

---

## 🐛 Troubleshooting

### Error: "No module named 'pennylane'"
**Solution**: Install dependencies
```bash
pip install pennylane torch optuna rdkit scikit-learn pandas matplotlib seaborn click
```

### Error: "CUDA not available" (GPU)
**Solution**: This is normal! Framework auto-falls back to CPU (`lightning.qubit`)

### Error: "Dataset download failed"
**Solution**: Check internet connection or use cached data:
```bash
mkdir -p .cache
# Place pre-downloaded CSV files in .cache/
```

---

## 📊 Expected Runtime

| Configuration | CPU Time | GPU Time | RAM |
|--------------|----------|----------|-----|
| 10 trials (test) | 30s | 10s | 2 GB |
| 100 trials | 5 min | 2 min | 4 GB |
| 500 trials | 30 min | 10 min | 8 GB |
| 2000 trials | 2 hours | 40 min | 16 GB |

---

## ✅ Verification Checklist

- [ ] `pip install -e .` runs without errors
- [ ] `vqc-drug-a1 --help` shows CLI options
- [ ] `pytest tests/ -v` passes all tests
- [ ] 10-trial smoke test completes successfully
- [ ] Output files generated in results directory

---

## 📖 Next Steps

1. **Read full documentation**: [README.md](README.md)
2. **Run extended pipeline**: 500+ trials for publication
3. **Generate figures**: Check `fig1_*.png` outputs (600 dpi)
4. **Submit to Zenodo**: Get DOI for reproducibility
5. **Prepare manuscript**: Use `latex_boilerplate.tex` as template

---

**Version**: 10.0-A1  
**Author**: Marcelo Claro Laranjeira  
**Contact**: marceloclaro@gmail.com  
**Status**: ✅ Production Ready
