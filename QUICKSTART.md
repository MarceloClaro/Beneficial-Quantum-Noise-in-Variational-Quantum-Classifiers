# MegaPrompt v2.0 - Quick Start Guide

**5 minutos para começar!**

## 🚀 Setup Rápido

```bash
# 1. Clone e configure
git clone https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers.git
cd Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers
pip install -r requirements.txt

# 2. Configure seu projeto
nano config.json  # Edite conforme necessário

# 3. Verifique estrutura existente
ls artigo_cientifico/fase*/
```

## ⚡ Comandos Essenciais

### Gerar Tabela S1
```bash
python tools/megaprompt_v2/generate_s1.py
```

### Consolidar Manuscrito
```bash
bash tools/megaprompt_v2/build_paper.sh
```

### Verificar Qualidade
```bash
# Consistência código-texto
python tools/megaprompt_v2/check_consistency.py

# Auditoria Qualis A1
python tools/megaprompt_v2/audit_checklist.py
```

### Workflow Completo
```bash
# Tudo de uma vez
python tools/megaprompt_v2/generate_s1.py && \
bash tools/megaprompt_v2/build_paper.sh && \
python tools/megaprompt_v2/check_consistency.py && \
python tools/megaprompt_v2/audit_checklist.py
```

## 📊 Verificar Resultados

```bash
# Ver manuscrito final
cat artigo_cientifico/fase6_consolidacao/manuscrito_internacional_final.md | head -100

# Ver sumário executivo
cat artigo_cientifico/fase6_consolidacao/sumario_executivo.md

# Ver checklist de auditoria
cat artigo_cientifico/fase6_consolidacao/checklist_auditoria_100pts.md

# Ver relatório de consistência
cat artigo_cientifico/fase6_consolidacao/relatorio_consistencia.md
```

## 🎯 Metas de Qualidade

- ✅ Consistência código-texto: **≥ 95%**
- ✅ Auditoria Qualis A1: **≥ 90/100 pontos**
- ✅ Manuscrito final: **≥ 20,000 palavras**
- ✅ Referências compiladas: **≥ 30 referências**

## 📚 Documentação Completa

- 📖 [README Principal](MEGAPROMPT_V2_README.md) - Visão geral completa
- 💡 [Exemplos Práticos](EXEMPLOS_PRATICOS.md) - Casos reais
- 🔧 [Workflow Exemplo](WORKFLOW_EXEMPLO.md) - Passo a passo detalhado
- 🛠️ [Documentação de Ferramentas](tools/megaprompt_v2/README.md)
- ❓ [FAQ](FAQ_TROUBLESHOOTING.md)
- 📚 [Glossário](GLOSSARIO.md)

## 🆘 Ajuda Rápida

### Problema: "config.json not found"
```bash
cp config.json.example config.json  # Se existir
# Ou crie manualmente
```

### Problema: Consistência baixa
```bash
# Veja o relatório detalhado
cat artigo_cientifico/fase6_consolidacao/relatorio_consistencia.md
# Corrija os issues apontados
```

### Problema: Build falha
```bash
# Verifique se todas as seções existem
ls -la artigo_cientifico/fase4_secoes/
# Se faltando, veja FAQ_TROUBLESHOOTING.md
```

## 🎓 Próximos Passos

1. **Leia**: [MEGAPROMPT_V2_README.md](MEGAPROMPT_V2_README.md)
2. **Veja**: [EXEMPLOS_PRATICOS.md](EXEMPLOS_PRATICOS.md)
3. **Execute**: [WORKFLOW_EXEMPLO.md](WORKFLOW_EXEMPLO.md)
4. **Submeta**: Seu artigo Qualis A1! 🚀

---

**Tempo estimado**: 6-10 dias úteis para artigo completo  
**Qualidade**: Pronto para Nature, Quantum, Physical Review

*MegaPrompt v2.0 - Framework para Artigos Científicos Qualis A1*
