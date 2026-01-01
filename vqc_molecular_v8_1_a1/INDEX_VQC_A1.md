# 🗂️ ÍNDICE E NAVEGAÇÃO - VQC-Molecular v8.1-A1

**Framework Quântico para Descoberta de Fármacos com Conformidade Qualis A1**

---

## 📍 Comece Aqui

### Se tem 5 minutos
👉 Leia [QUICKSTART_VQC_A1.md](QUICKSTART_VQC_A1.md)
- Setup (conda)
- Executar estudo piloto
- Revisar resultados

### Se tem 30 minutos
👉 Leia [README_VQC_DRUG_A1.md](README_VQC_DRUG_A1.md)
- Visão geral completa
- Módulos Qualis A1 explicados
- Exemplos de uso
- FAQ

### Se quer verificar implementação
👉 Leia [IMPLEMENTATION_VERIFICATION_VQC_A1.md](IMPLEMENTATION_VERIFICATION_VQC_A1.md)
- Checklist de funcionalidades
- Estrutura de saída
- Testes recomendados

---

## 📚 ESTRUTURA DE ARQUIVOS

### 🐍 Código Python (Core)

```
vqc_drug_qualis_a1.py          ← MAIN: Pipeline completo Qualis A1
├── preregister.py              ← Pré-registro com SHA-256
├── audit.py                    ← Checksums integridade
├── power_analysis.py           ← Análise de poder
├── statistics.py               ← Testes múltiplos + effect sizes
├── figures.py                  ← Figuras 600 dpi
└── supp_tables.py              ← Tabelas suplementares Excel
```

**Como usar:**
```bash
python vqc_drug_qualis_a1.py --target EGFR --trials 300
```

---

### 📋 Documentação

```
README_VQC_DRUG_A1.md               ← Referência COMPLETA
QUICKSTART_VQC_A1.md                ← Guia RÁPIDO (1 hora)
IMPLEMENTATION_VERIFICATION_VQC_A1.md ← Verificação
INDEX_VQC_A1.md                     ← Este arquivo (navegação)
```

**Qual ler primeiro?**
- Novo usuário? → `QUICKSTART_VQC_A1.md`
- Quer detalhes? → `README_VQC_DRUG_A1.md`
- Quer verificar? → `IMPLEMENTATION_VERIFICATION_VQC_A1.md`

---

### ⚙️ Configuração

```
environment.yml                 ← Conda environment (recomendado)
Dockerfile                      ← Docker (reprodutibilidade)
requirements_drug_screening.txt ← Pip requirements (alternativo)
```

**Como usar:**
```bash
# Opção 1: Conda
conda env create -f environment.yml && conda activate vqc-a1

# Opção 2: Docker
docker build -t vqc-a1:8.1 . && docker run -it vqc-a1:8.1

# Opção 3: Pip
pip install -r requirements_drug_screening.txt
```

---

### 📊 Scripts Executáveis (v8.0 legado, ainda disponível)

```
run_vqc_drug_examples.sh        ← Menu interativo (Bash/Linux/macOS)
run_vqc_drug_examples.ps1       ← Menu interativo (PowerShell/Windows)
```

**Uso:**
```bash
bash run_vqc_drug_examples.sh      # Linux/macOS
./run_vqc_drug_examples.ps1        # Windows PowerShell
```

---

## 🎯 FLUXO DE TRABALHO TÍPICO

### 1️⃣ SETUP (5 minutos)
```bash
# Instalar ambiente
conda env create -f environment.yml
conda activate vqc-a1

# Verificar
python -c "import pennylane; print('✅')"
```

### 2️⃣ PRÉ-EXPERIMENTO (2 minutos)
```python
# Revisar power analysis
from power_analysis import required_sample_size
n = required_sample_size(effect_size=0.35, alpha=0.05, power=0.8)
# → 151 amostras por grupo
```

### 3️⃣ EXECUTAR ESTUDO (45 minutos GPU / 120 min CPU)
```bash
# Estudo completo (publicação)
python vqc_drug_qualis_a1.py --target EGFR --trials 300

# Ou teste rápido
python vqc_drug_qualis_a1.py --target EGFR --trials 50
```

### 4️⃣ REVISAR RESULTADOS (5 minutos)
```bash
# Ver resultado principal
cat results_*/final_report_EGFR.json | python -m json.tool

# Visualizar figuras
open results_*/04_figuras_publicacao/fig*.png

# Checar integridade
python audit.py
```

### 5️⃣ PREPARAR PUBLICAÇÃO (30 minutos)
```bash
# Gerar DOI em Zenodo
tar czf vqc_molecular_v8.1.tgz results_*/
# Upload em zenodo.org

# Escrever paper usando
# - JSON como dados (results_*/final_report_*.json)
# - Figuras como figures (04_figuras_publicacao/)
# - Tabelas como supplementary (05_tabelas_suplementares/)
```

---

## 🔬 MÓDULOS QUALIS A1 - GUIA TÉCNICO

### `preregister.py` - Pré-Registro
**O quê**: Cria protocolo imutável com SHA-256  
**Por quê**: Detecta mudanças de hipóteses post-hoc  
**Como usar**:
```python
from preregister import pre_register
proto_file, hash = pre_register(target="EGFR", n_trials=300)
```

---

### `audit.py` - Auditoria
**O quê**: Checksum SHA-256 de todos arquivos  
**Por quê**: Verifica integridade de dados (detecção de tampering)  
**Como usar**:
```python
from audit import hash_all_files, verify_checksums
hash_all_files(".", "checksums.sha256")
is_valid = verify_checksums("checksums.sha256")[0]
```

---

### `power_analysis.py` - Análise de Poder
**O quê**: Calcula tamanho amostral necessário  
**Por quê**: Garante poder estatístico suficiente (β < 0.2)  
**Como usar**:
```python
from power_analysis import required_sample_size
n = required_sample_size(effect_size=0.35, alpha=0.05, power=0.8)
```

---

### `statistics.py` - Testes Múltiplos + Effect Sizes
**O quê**: Bonferroni-Holm, FDR, Cohen d com IC 95% bootstrap  
**Por quê**: Controla FPR, quantifica tamanho de efeito com incerteza  
**Como usar**:
```python
from statistics import ttest_with_correction, cohen_d_with_bootstrap_ci
results = ttest_with_correction(vqc_aucs, baseline_aucs)
d_info = cohen_d_with_bootstrap_ci(vqc_aucs, baseline_aucs)
```

---

### `figures.py` - Figuras 600 dpi
**O quê**: Gera 4 figuras automáticas em alta resolução  
**Por quê**: Pronto para Nature/Quantum (sem edição)  
**Como usar**:
```python
from figures import fig_power_curve, fig_roc_comparison
fig_power_curve(output_file="fig1_power.png")  # 600 dpi PNG
```

---

### `supp_tables.py` - Tabelas Suplementares
**O quê**: Gera 4 tabelas Excel formatadas  
**Por quê**: Pronto para publicação (trials, statísticas, effects, params)  
**Como usar**:
```python
from supp_tables import generate_all_supplementary_tables
tables = generate_all_supplementary_tables(study, results_df, effect_sizes)
```

---

## 📁 ESTRUTURA DE SAÍDA (após execução)

```
results_2025-01-15_14-30-45/
│
├── 01_protocolo_pre_registrado/
│   └── 01_protocolo_pre_registrado_a1b2c3d4.json
│       → Protocolo imutável com hash SHA-256
│
├── 02_dados_brutos/
│   ├── raw_EGFR.csv
│   └── checksums.sha256
│       → Dados originais + checksums para auditoria
│
├── 03_analises_estatisticas/
│   ├── cv_comparison.csv
│   ├── bonferroni_holm.csv
│   └── effect_sizes_bootstrap.csv
│       → Testes múltiplos, effect sizes, IC 95%
│
├── 04_figuras_publicacao/
│   ├── fig1_power_curve.png (600 dpi)
│   ├── fig2_roc_comparison.png (600 dpi)
│   ├── fig3_forest_effect_sizes.png (600 dpi)
│   └── fig4_optuna_trials.png (600 dpi)
│       → Pronto para Nature (nenhuma edição necessária)
│
├── 05_tabelas_suplementares/
│   ├── supp_table1_trials_complete.xlsx
│   ├── supp_table2_statistical_tests.xlsx
│   ├── supp_table3_effect_sizes.xlsx
│   └── supp_table4_best_hyperparameters.xlsx
│       → Tabelas Excel formatadas para publicação
│
├── 06_reprodutibilidade/
│   ├── environment.yml
│   ├── dockerfile
│   └── codigo_hash.sha256
│       → Ferramentas para reprodução exata
│
├── 07_log_execucao/
│   ├── vqc_execution_20250115_143045.log
│   ├── audit_report.json
│   └── resumo_qualis_a1.json
│       → Logs DEBUG completos + relatório auditoria
│
├── final_report_EGFR.json
│   → RESULTADO PRINCIPAL (dados estruturados)
│
├── optuna_history.html
│   → Visualização interativa (Plotly)
│
└── checksums_final.sha256
    → Checksums finais para integridade
```

---

## 🎓 EXEMPLOS DE CÓDIGO

### Exemplo 1: Estudo Piloto Rápido (10 minutos)
```bash
python vqc_drug_qualis_a1.py --target EGFR --trials 50
```

### Exemplo 2: Estudo Completo (45 minutos GPU)
```bash
python vqc_drug_qualis_a1.py --target EGFR --trials 300 --max-qubits 20
```

### Exemplo 3: Todos os Targets
```bash
for target in EGFR HIV Malaria COVID; do
    python vqc_drug_qualis_a1.py --target $target --trials 300
done
```

### Exemplo 4: Com Saída Customizada
```bash
python vqc_drug_qualis_a1.py \
    --target EGFR \
    --trials 500 \
    --max-qubits 20 \
    --seed 42 \
    --out-dir my_results
```

### Exemplo 5: Reproduzir via Docker
```bash
docker build -t vqc-a1:final .
docker run -v /results:/app/results vqc-a1:final \
    python vqc_drug_qualis_a1.py --target EGFR
```

---

## 📞 TROUBLESHOOTING RÁPIDO

| Problema | Solução |
|----------|---------|
| `ModuleNotFoundError: pennylane` | `pip install pennylane` |
| GPU out of memory | Reduza `--max-qubits 12` |
| Estudo muito lento | Use `--trials 50` para teste |
| RDKit not found | `conda install -c conda-forge rdkit` |
| DeepChem erro | Instale: `pip install deepchem` |

Veja `README_VQC_DRUG_A1.md` seção FAQ para mais detalhes.

---

## ✅ CHECKLIST PRÉ-SUBMISSÃO

- [ ] Executar `python vqc_drug_qualis_a1.py --target EGFR --trials 300`
- [ ] Verificar `final_report_EGFR.json` (Cohen d > 0.3?)
- [ ] Revisar `07_log_execucao/audit_report.json` (VALID?)
- [ ] Conferir 4 figuras: `04_figuras_publicacao/fig*.png`
- [ ] Revisar 4 tabelas: `05_tabelas_suplementares/`
- [ ] Testar Docker: `docker build . && docker run ...`
- [ ] Gerar DOI em Zenodo: `tar czf vqc.tgz results_*/ && zenodo upload`
- [ ] Escrever cover letter mencionando:
  - [ ] Pre-registration with audit trail
  - [ ] Power analysis (power ≥ 0.8)
  - [ ] Multiple comparison corrections (Bonferroni-Holm)
  - [ ] Effect sizes with 95% bootstrap CI (10,000 iterations)
  - [ ] Reproducible via Docker
  - [ ] DOI at Zenodo

---

## 🚀 PRÓXIMOS PASSOS

1. **Hoje**: Instalar + teste piloto (1 hora)
2. **Amanhã**: Executar estudo completo + revisar (2 horas)
3. **Esta semana**: Escrever paper + cover letter (2 dias)
4. **Próxima semana**: Submeter para Nature/Quantum (1 hora)

---

## 📖 RECURSOS ADICIONAIS

- **PennyLane Docs**: https://pennylane.ai
- **Optuna Docs**: https://optuna.org
- **DeepChem Docs**: https://deepchem.io
- **Nature Publishing**: https://www.nature.com/nature
- **Zenodo**: https://zenodo.org

---

## 📄 Informações

**Framework**: VQC-Molecular v8.1-A1  
**Versão**: 8.1-A1  
**Data**: 2025-01-15  
**Status**: ✅ Pronto para Produção  
**Licença**: MIT (código) + CC-BY (dados)  

**Autores**: Quantum Drug Discovery Team

---

**Última atualização**: 2025-01-15  
**Próxima versão**: 8.2-A1 (MPNN hybrid + multi-task learning)
