# MegaPrompt v2.0 Tools

Automation tools for scientific article generation with Qualis A1 compliance.

## 📦 Available Tools

### 1. generate_s1.py
Generate Table S1 with complete experimental configuration matrix.

**Usage:**

```bash
python generate_s1.py --config config.json --output artigo_cientifico/fase5_suplementar/tabela_s1_configuracoes.csv

```text

#### Features:
- Calculates all factor combinations
- Generates CSV with configuration details
- Validates total count against expected value
- Adds computed fields (framework, n_qubits, etc.)


**Example Output:**

```

Total configurations: 2,688
Columns: 15
✓ Table S1 generated: artigo_cientifico/fase5_suplementar/tabela_s1_configuracoes.csv

```text

### 2. check_consistency.py
Verify consistency between code implementation and manuscript text.

**Usage:**

```bash
python check_consistency.py --manuscript artigo_cientifico/ --code . --output consistency_report.md

```text

#### Checks Performed:
- Numeric claims without evidence
- Missing citations
- Referenced code files not found
- Integrity markers ([INFORMAÇÃO AUSENTE], etc.)


**Example Output:**

```

Consistency: 96.5%
✅ Status: EXCELENTE - Meta atingida (≥95%)

```text

### 3. build_paper.sh
Consolidate all sections into final manuscript.

**Usage:**

```bash
bash build_paper.sh

```text

#### Features:
- Concatenates all sections in order
- Adds supplementary material references
- Generates executive summary
- Validates output
- Supports MODE_A (English/LaTeX) and MODE_B (Portuguese/ABNT)


#### Output:
- `manuscrito_internacional_final.md` (MODE_A)
- `artigo_abnt_final.md` (MODE_B)
- `sumario_executivo.md` (always generated)


### 4. audit_checklist.py
Comprehensive 0-100 point audit checklist.

**Usage:**

```bash
python audit_checklist.py --config config.json --manuscript artigo_cientifico/

```text

**Categories:**
1. **Reproducibility (30 pts)**
   - Environment documented
   - Seeds fixed and reported
   - Pipeline executable


2. **Traceability (30 pts)**
   - Traceability table complete
   - Code→Method mapping complete


3. **Statistical Rigor (20 pts)**
   - Appropriate tests
   - Multiple comparison correction
   - Confidence intervals
   - Effect sizes


4. **Transparency (20 pts)**
   - Code publicly available
   - Data availability
   - Limitations discussed


#### Scoring:
- ≥ 90: EXCELENTE - Ready for submission
- 80-89: BOM - Minor adjustments needed
- 70-79: ACEITÁVEL - Significant improvements needed
- < 70: INSUFICIENTE - Complete revision required


## 🚀 Quick Start Workflow

```bash

# 1. Configure your project
nano config.json

# 2. Generate Table S1
python tools/megaprompt_v2/generate_s1.py

# 3. Build consolidated manuscript
bash tools/megaprompt_v2/build_paper.sh

# 4. Check consistency
python tools/megaprompt_v2/check_consistency.py

# 5. Run audit checklist
python tools/megaprompt_v2/audit_checklist.py

```text

## 📋 Requirements

All tools require Python 3.7+ and standard library only (no additional dependencies).

## 🔧 Advanced Usage

### Verify Configuration Count Only

```bash
python generate_s1.py --verify-only

```text

### Custom Output Paths

```bash
python check_consistency.py \

    --manuscript artigo_cientifico/ \
    --code . \
    --output reports/consistency_$(date +%Y%m%d).md

```text

### Automated Pipeline

Create a `Makefile`:

```makefile
.PHONY: all s1 build check audit

all: s1 build check audit

s1:
	python tools/megaprompt_v2/generate_s1.py

build:
	bash tools/megaprompt_v2/build_paper.sh

check:
	python tools/megaprompt_v2/check_consistency.py

audit:
	python tools/megaprompt_v2/audit_checklist.py

clean:
	rm -f artigo_cientifico/fase6_consolidacao/*.md
	rm -f artigo_cientifico/fase5_suplementar/*.csv

```text

Then run:

```bash
make all

```text

## 📊 Output Files

After running all tools, you'll have:

```

artigo_cientifico/
├── fase5_suplementar/
│   └── tabela_s1_configuracoes.csv       (generate_s1.py)
├── fase6_consolidacao/
│   ├── manuscrito_internacional_final.md  (build_paper.sh)
│   ├── sumario_executivo.md               (build_paper.sh)
│   ├── relatorio_consistencia.md          (check_consistency.py)
│   └── checklist_auditoria_100pts.md      (audit_checklist.py)

```text

## 🐛 Troubleshooting

### "config.json not found"
Make sure you're running from the repository root directory, or specify the full path:

```bash
python tools/megaprompt_v2/generate_s1.py --config /path/to/config.json

```text

### "Manuscript directory not found"
Ensure the `artigo_cientifico/` directory exists and has the expected structure.

### Build script fails
Check that all required section files exist:

```bash
ls -la artigo_cientifico/fase4_secoes/

```text

### Low consistency score
Review the consistency report for specific issues:

```bash
cat artigo_cientifico/fase6_consolidacao/relatorio_consistencia.md

```

## 📚 Related Documentation

- Main framework: [../README.md](../../README.md)
- MegaPrompt v2.0 overview: [../MEGAPROMPT_V2_README.md](../../MEGAPROMPT_V2_README.md)
- Configuration guide: [../config.json](../../config.json)
- FAQ: [../FAQ_TROUBLESHOOTING.md](../../FAQ_TROUBLESHOOTING.md)
- Glossary: [../GLOSSARIO.md](../../GLOSSARIO.md)


## 🤝 Contributing

To add new tools:

1. Follow the same structure as existing tools
2. Add command-line argument parsing
3. Include comprehensive docstrings
4. Update this README
5. Add to the workflow section


## 📄 License

Same as parent project (see main LICENSE file).

---


**MegaPrompt v2.0 Framework**  
*Tools for Qualis A1 Scientific Article Generation*

