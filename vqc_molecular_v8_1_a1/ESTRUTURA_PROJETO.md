# VQC-Molecular v8.1-A1 - Estrutura do Projeto

## 📁 Organização de Diretórios

```
vqc_molecular_v8_1_a1/
├── run_vqc_a1.py                              # 🚀 Ponto de entrada principal
├── vqc_drug_qualis_a1.py                      # 📦 Pipeline integrado
├── preregister.py                              # 📝 Pré-registro SHA-256
├── audit.py                                    # 🔐 Auditoria de checksums
├── power_analysis.py                           # 📊 Análise de poder
├── statistics.py                               # 📈 Testes estatísticos
├── figures.py                                  # 🎨 Figuras 600 dpi
├── supp_tables.py                              # 📋 Tabelas suplementares
├── environment.yml                             # 🐍 Conda environment
├── Dockerfile                                  # 🐳 Docker reproducibilidade
├── README_VQC_DRUG_A1.md                      # 📖 Guia de referência
├── QUICKSTART_VQC_A1.md                       # ⚡ Início rápido 1h
├── INDEX_VQC_A1.md                            # 🗂️  Navegação
├── IMPLEMENTATION_VERIFICATION_VQC_A1.md      # ✅ Checklist
└── ESTRUTURA_PROJETO.md                       # 📍 Este arquivo

results_EGFR/                                   # 📂 Saídas após execução
├── 01_protocolo_pre_registrado/
├── 02_analise_poder/
├── 03_testes_estatisticos/
├── 04_figuras_publicacao/                     # Figuras 600 dpi PNG
├── 05_tabelas_suplementares/                  # Tabelas Excel
├── 06_dados_json/
└── 07_log_execucao/
```

## 🚀 Como Executar

### Opção 1: Execução Rápida (Recomendado)
```bash
cd vqc_molecular_v8_1_a1/
python run_vqc_a1.py --target EGFR --trials 300
```

### Opção 2: Execução Direta
```bash
cd vqc_molecular_v8_1_a1/
python vqc_drug_qualis_a1.py --target EGFR --trials 300
```

### Opção 3: Via Conda (Ambiente Isolado)
```bash
cd vqc_molecular_v8_1_a1/
conda env create -f environment.yml
conda activate vqc-a1
python run_vqc_a1.py --target EGFR --trials 300
```

### Opção 4: Via Docker (Reprodutibilidade Total)
```bash
cd vqc_molecular_v8_1_a1/
docker build -t vqc-a1:latest .
docker run -v $(pwd)/results:/app/results vqc-a1:latest python run_vqc_a1.py --target EGFR
```

## 📚 Documentação

| Arquivo | Público-Alvo | Tempo |
|---------|--------------|--------|
| **QUICKSTART_VQC_A1.md** | Usuários Impatientes | 1 hora |
| **README_VQC_DRUG_A1.md** | Referência Completa | 2-3 horas |
| **INDEX_VQC_A1.md** | Navegação + Exemplos | 30 minutos |
| **IMPLEMENTATION_VERIFICATION_VQC_A1.md** | Checklist de Testes | 45 minutos |

## 🔧 Módulos Principais

### `vqc_drug_qualis_a1.py` (600+ linhas)
Pipeline principal com 6 estágios:
1. Pré-registro SHA-256
2. Download e preparação de dados QSAR
3. Análise de poder estatístico
4. Otimização VQC com Optuna
5. Treinamento do modelo
6. Relatórios Qualis A1

### `preregister.py` (184 linhas)
- `pre_register()` - Criar protocolo imutável
- `validate_preregistration()` - Validar integridade

### `audit.py` (280 linhas)
- `hash_file()` - SHA-256 para arquivo único
- `hash_all_files()` - SHA-256 recursivo para diretório
- `verify_checksums()` - Verificar integridade
- `audit_report()` - Gerar relatório JSON

### `power_analysis.py` (270 linhas)
- `required_sample_size()` - Calcular n necessário
- `power_curve()` - Curva de poder
- `plot_power_curve()` - Visualização 600 dpi
- `sample_size_table()` - Tabela ASCII

### `statistics.py` (307 linhas)
- `cohen_d_with_bootstrap_ci()` - Cohen d + IC 95%
- `hedges_g()` - Hedges g (não-enviesado)
- `glass_delta()` - Glass Δ
- `bonferroni_holm_correction()` - Correção múltipla
- `fdr_benjamini_hochberg()` - Controle FDR
- `ttest_with_correction()` - T-teste completo

### `figures.py` (340 linhas)
Gera 4 tipos de figuras publicáveis:
1. Curva de Poder (`fig_power_curve()`)
2. ROC Comparison (`fig_roc_comparison()`)
3. Forest Plot (`fig_forest_plot()`)
4. Optuna History (`fig_optuna_history()`)

Todas em 600 dpi PNG automático.

### `supp_tables.py` (318 linhas)
Gera 4 tabelas suplementares Excel:
1. Todas as tentativas Optuna (300+ linhas)
2. Testes estatísticos VQC vs Baseline
3. Effect sizes com IC 95%
4. Melhores hiperparâmetros

## 📊 Saídas Esperadas

Após executar `run_vqc_a1.py --target EGFR --trials 300`:

```
results_EGFR/
├── 01_protocolo_pre_registrado/
│   └── protocolo_pre_registrado_<HASH>.json
├── 02_analise_poder/
│   ├── poder_analise.png (600 dpi)
│   ├── poder_relatorio.json
│   └── amostra_size_table.txt
├── 03_testes_estatisticos/
│   ├── comparacao_vqc_baseline.json
│   └── effect_sizes_completo.json
├── 04_figuras_publicacao/
│   ├── fig1_poder.png (600 dpi)
│   ├── fig2_roc_comparison.png (600 dpi)
│   ├── fig3_forest_plot.png (600 dpi)
│   └── fig4_optuna_history.png (600 dpi)
├── 05_tabelas_suplementares/
│   ├── supp_table1_trials_complete.xlsx
│   ├── supp_table2_statistical_tests.xlsx
│   ├── supp_table3_effect_sizes.xlsx
│   └── supp_table4_hyperparameters_best.xlsx
├── 06_dados_json/
│   ├── dataset_egfr.json
│   ├── best_hyperparameters.json
│   └── final_metrics.json
└── 07_log_execucao/
    ├── execution_log_qualis_a1.log
    └── audit_report.json (status: VALID)
```

## 🎯 Fluxo de Trabalho Típico

```
1. Setup (5 min)
   └─ conda env create -f environment.yml

2. Execução (45 min GPU / 120 min CPU)
   └─ python run_vqc_a1.py --target EGFR --trials 300

3. Revisão (5 min)
   ├─ Verificar JSON outputs
   ├─ Revisar figuras 600 dpi
   ├─ Inspecionar tabelas Excel
   └─ Validar audit_report.json (VALID?)

4. Submissão (2 horas)
   ├─ Escrever paper usando JSON
   ├─ Incluir figuras de 04_figuras_publicacao/
   ├─ Anexar tabelas de 05_tabelas_suplementares/
   └─ Citar pré-registro de 01_protocolo_pre_registrado/
```

## 🔍 Verificação Rápida

Para validar a instalação:
```bash
cd vqc_molecular_v8_1_a1/
python -c "
import vqc_drug_qualis_a1
import preregister
import audit
import power_analysis
import statistics
import figures
import supp_tables
print('✅ Todos os módulos carregados com sucesso!')
"
```

## ❓ Troubleshooting

### Erro: "ModuleNotFoundError: No module named 'pennylane'"
**Solução**: 
```bash
conda env create -f environment.yml
conda activate vqc-a1
```

### Erro: "ImportError: cannot import name 'VQCMolecularAudit'"
**Solução**: 
Certifique-se de estar no diretório `vqc_molecular_v8_1_a1/` antes de executar.

### Execução Muito Lenta (CPU)
**Solução**: 
```bash
# Usar GPU (se disponível)
export PENNYLANE_PLUGIN_PATH=qiskit_aer
# Ou reduzir número de tentativas
python run_vqc_a1.py --target EGFR --trials 100
```

## 📖 Próximos Passos

1. **Ler QUICKSTART_VQC_A1.md** (1 hora)
2. **Executar `run_vqc_a1.py` piloto** (45 min)
3. **Revisar outputs** (5 min)
4. **Escrever paper** (2 horas)
5. **Submeter para Qualis A1** ✨

---

**Version**: v8.1-A1  
**Last Updated**: 30 de Dezembro de 2025  
**Status**: ✅ Production Ready
