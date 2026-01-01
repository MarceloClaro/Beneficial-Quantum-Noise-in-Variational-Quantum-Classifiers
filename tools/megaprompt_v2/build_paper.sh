#!/bin/bash
# Build Paper Script - Consolidate all sections into final manuscript
# MegaPrompt v2.0 Framework
# Date: 2025-12-26

set -e  # Exit on error

echo "======================================================================"
echo "MEGAPROMPT V2.0 - BUILD PAPER SCRIPT"
echo "======================================================================"
echo ""

# Configuration
ARTIGO_DIR="artigo_cientifico"
OUTPUT_DIR="${ARTIGO_DIR}/fase6_consolidacao"
CONFIG_FILE="config.json"

# Check if config exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "❌ Error: config.json not found"
    exit 1
fi

# Read output mode from config
OUTPUT_MODE=$(python3 -c "import json; print(json.load(open('$CONFIG_FILE'))['output_mode'])" 2>/dev/null || echo "MODE_A")

echo "Configuration:"
echo "  Output Mode: $OUTPUT_MODE"
echo "  Article Directory: $ARTIGO_DIR"
echo "  Output Directory: $OUTPUT_DIR"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Step 1: Consolidating sections..."

# Define section files
SECTIONS=(
    "fase4_secoes/resumo_abstract.md"
    "fase4_secoes/introducao_completa.md"
    "fase4_secoes/revisao_literatura_completa.md"
    "fase4_secoes/metodologia_completa.md"
    "fase4_secoes/resultados_completo.md"
    "fase4_secoes/discussao_completa.md"
    "fase4_secoes/conclusao_completa.md"
    "fase4_secoes/agradecimentos_referencias.md"
)

# Check all sections exist
MISSING=0
for section in "${SECTIONS[@]}"; do
    if [ ! -f "${ARTIGO_DIR}/${section}" ]; then
        echo "  ⚠️  Missing: $section"
        MISSING=$((MISSING + 1))
    else
        echo "  ✓ Found: $section"
    fi
done

if [ $MISSING -gt 0 ]; then
    echo ""
    echo "❌ Error: $MISSING section(s) missing"
    exit 1
fi

echo ""
echo "Step 2: Building consolidated manuscript..."

# Create consolidated manuscript
if [ "$OUTPUT_MODE" = "MODE_A" ]; then
    # MODE_A: English/LaTeX
    OUTPUT_FILE="${OUTPUT_DIR}/manuscrito_internacional_final.md"
    echo "% Beneficial Quantum Noise in Variational Quantum Classifiers" > "$OUTPUT_FILE"
    echo "% Research Team" >> "$OUTPUT_FILE"
    echo "% $(date +%Y-%m-%d)" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
else
    # MODE_B: Portuguese/ABNT
    OUTPUT_FILE="${OUTPUT_DIR}/artigo_abnt_final.md"
    echo "# RUÍDO QUÂNTICO BENÉFICO EM CLASSIFICADORES QUÂNTICOS VARIACIONAIS" > "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
fi

# Concatenate all sections
for section in "${SECTIONS[@]}"; do
    section_name=$(basename "$section" .md)
    echo "  📄 Adding: $section_name"
    
    echo "" >> "$OUTPUT_FILE"
    echo "<!-- Section: $section_name -->" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
    cat "${ARTIGO_DIR}/${section}" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
done

echo "  ✓ Consolidated manuscript: $OUTPUT_FILE"

echo ""
echo "Step 3: Adding supplementary material references..."

# Add reference to supplementary material
echo "<!-- Supplementary Material -->" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"
echo "## Supplementary Material" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"
echo "See separate supplementary files:" >> "$OUTPUT_FILE"
echo "- Table S1-S5: \`fase5_suplementar/tabelas_suplementares.md\`" >> "$OUTPUT_FILE"
echo "- Figures S1-S8: \`fase5_suplementar/figuras_suplementares.md\`" >> "$OUTPUT_FILE"
echo "- Additional Notes: \`fase5_suplementar/notas_metodologicas_adicionais.md\`" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"

echo "  ✓ Supplementary references added"

echo ""
echo "Step 4: Generating statistics..."

# Count words
WORD_COUNT=$(wc -w < "$OUTPUT_FILE")
LINE_COUNT=$(wc -l < "$OUTPUT_FILE")

echo "  Words: $WORD_COUNT"
echo "  Lines: $LINE_COUNT"

echo ""
echo "Step 5: Creating executive summary..."

# Create executive summary
SUMMARY_FILE="${OUTPUT_DIR}/sumario_executivo.md"

cat > "$SUMMARY_FILE" << EOF
# Sumário Executivo - Manuscrito Final

**Data de Geração**: $(date +"%Y-%m-%d %H:%M:%S")

## 📊 Estatísticas do Manuscrito

- **Arquivo Principal**: \`$(basename "$OUTPUT_FILE")\`
- **Total de Palavras**: $WORD_COUNT
- **Total de Linhas**: $LINE_COUNT
- **Seções Incluídas**: ${#SECTIONS[@]}
- **Modo de Saída**: $OUTPUT_MODE

## 📁 Estrutura do Artigo

$(for i in "${!SECTIONS[@]}"; do
    num=$((i + 1))
    section_name=$(basename "${SECTIONS[$i]}" .md)
    echo "$num. ${section_name//_/ }"
done)

## 📎 Material Suplementar

- Tabelas S1-S5 (5 tabelas)
- Figuras S1-S8 (8 figuras)
- Notas Metodológicas Adicionais
- Tabela de Rastreabilidade Completa

## ✅ Próximos Passos

1. **Revisão Final**: Ler o manuscrito completo
2. **Verificação de Consistência**: Executar \`check_consistency.py\`
3. **Validação LaTeX**: Compilar template LaTeX (se MODE_A)
4. **Preparação para Submissão**: Criar ZIP com todos os arquivos

## 🎯 Checklist Pré-Submissão

- [ ] Manuscrito revisado e aprovado
- [ ] Consistência código-texto ≥ 95%
- [ ] Todas as figuras e tabelas incluídas
- [ ] Referências completas e formatadas
- [ ] Material suplementar organizado
- [ ] README e instruções de reprodução incluídos
- [ ] Código disponibilizado publicamente
- [ ] Checklist de auditoria preenchido (≥90/100 pts)

---

*Gerado automaticamente por MegaPrompt v2.0 Build Script*
EOF

echo "  ✓ Executive summary: $SUMMARY_FILE"

echo ""
echo "Step 6: Validating output..."

# Basic validation
if [ ! -f "$OUTPUT_FILE" ]; then
    echo "  ❌ Error: Output file not created"
    exit 1
fi

if [ $WORD_COUNT -lt 5000 ]; then
    echo "  ⚠️  Warning: Word count seems low ($WORD_COUNT < 5000)"
fi

echo "  ✓ Validation passed"

echo ""
echo "======================================================================"
echo "✅ BUILD COMPLETE!"
echo "======================================================================"
echo ""
echo "Output Files:"
echo "  📄 Manuscript: $OUTPUT_FILE"
echo "  📋 Summary: $SUMMARY_FILE"
echo ""
echo "Next Steps:"
echo "  1. Review the consolidated manuscript"
echo "  2. Run consistency checker:"
echo "     python3 tools/megaprompt_v2/check_consistency.py"
echo "  3. Generate Table S1:"
echo "     python3 tools/megaprompt_v2/generate_s1.py"
echo "  4. Complete final audit checklist"
echo ""
echo "For LaTeX compilation (MODE_A):"
echo "  cd $ARTIGO_DIR/latex_template"
echo "  pdflatex npj_qi_submission.tex"
echo ""
