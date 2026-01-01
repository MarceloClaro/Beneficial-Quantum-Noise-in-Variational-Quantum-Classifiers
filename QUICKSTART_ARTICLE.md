# 🚀 Quick Start - Article Generation Framework

## One-Command Validation

```bash

# Check if your article is ready for submission
python tools/validate_qualis_a1.py --article artigo_cientifico/ --report VALIDATION.md && \
python tools/verify_code_text_congruence.py --code framework_investigativo_completo.py --article artigo_cientifico/ --output CONGRUENCE.md

```text

## Interpretation

### Validation Score
- **90-100:** 🥇 Submit to Nature/Science/Quantum immediately
- **80-89:** 🥈 Submit to Qualis A1 journals  
- **70-79:** 🥉 Minor revisions needed
- **<70:** ❌ Substantial revision required


### Congruence Score
- **≥95%:** ✅ Perfect reproducibility
- **85-94%:** ✅ Good alignment
- **70-84%:** ⚠️ Some inconsistencies
- **<70%:** ❌ Major alignment issues


## Common Issues & Quick Fixes

### Low Validation Score

**Issue:** Word counts outside expected ranges


**Fix:**

```bash

# Check which sections are problematic
cat VALIDATION.md | grep "⚠️\|❌"

# Edit the corresponding files in artigo_cientifico/fase4_secoes/
# Expand short sections or condense long ones

```text

### Low Congruence Score

**Issue:** Code components not mentioned in article


**Fix:**

```bash

# Check specific inconsistencies
cat CONGRUENCE.md | grep "❌"

# Update Methodology section to include all implemented components
nano artigo_cientifico/fase4_secoes/metodologia_completa.md

```text

### Missing References

**Issue:** Less than 35 references


**Fix:**

```bash

# Add references to bibliografia
nano artigo_cientifico/fase2_bibliografia/referencias_compiladas.md

# Update citations in article sections
nano artigo_cientifico/fase4_secoes/revisao_literatura_completa.md

```text

## Full Documentation

- **Complete Guide:** [WORKFLOW_ARTIGO.md](WORKFLOW_ARTIGO.md)
- **Specification:** [MEGA_PROMPT_QUALIS_A1.md](MEGA_PROMPT_QUALIS_A1.md)
- **Article Folder:** [artigo_cientifico/](artigo_cientifico/)


## Target Journals

1. **npj Quantum Information** ⭐⭐⭐ (Recommended)
2. **Nature Communications** ⭐⭐⭐
3. **Quantum** ⭐⭐
4. **Physical Review A** ⭐⭐


## Need Help?

```bash

# View all available commands
python tools/validate_qualis_a1.py --help
python tools/verify_code_text_congruence.py --help
python gerador_artigo_completo.py --help

```

---


**Framework Version:** 1.0  
**Status:** Production Ready  
**License:** MIT

