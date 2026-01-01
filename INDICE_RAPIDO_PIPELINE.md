# 🔗 Índice Rápido: Pipeline Automatizada QAOA

**Comece aqui!** Escolha seu caminho:

---

## 🏃 Tenho 1 minuto

```powershell
python test_pipeline_verificacao.py
```

**Resultado esperado**: ✅ PIPELINE FUNCIONANDO CORRETAMENTE!

Veja: [QUICKSTART_PIPELINE_QAOA.md](QUICKSTART_PIPELINE_QAOA.md) (2 min)

---

## 🚶 Tenho 5 minutos

1. Ler: [VISAO_GERAL_PIPELINE.md](VISAO_GERAL_PIPELINE.md)
2. Rodar: `python test_pipeline_verificacao.py`
3. Entender: 5 etapas de automatização

---

## 📖 Tenho 15 minutos

1. Ler: [PLAYBOOK_REPRODUCAO_QAOA.md](PLAYBOOK_REPRODUCAO_QAOA.md) (passo a passo)
2. Executar: Pipeline completa (4-7 min)
3. Verificar: Gráficos em `auditoria_qaoa/auditoria_qaoa_energia.html`

---

## 🔧 Preciso entender tudo

**Leia nesta ordem:**

1. **[QUICKSTART_PIPELINE_QAOA.md](QUICKSTART_PIPELINE_QAOA.md)** (2 min)
   - O que cada script faz
   - Onde estão os resultados

2. **[VISAO_GERAL_PIPELINE.md](VISAO_GERAL_PIPELINE.md)** (5 min)
   - Diagrama ASCII
   - Checklist de uso
   - Tempo esperado

3. **[PLAYBOOK_REPRODUCAO_QAOA.md](PLAYBOOK_REPRODUCAO_QAOA.md)** (10 min)
   - Reprodução passo a passo
   - Saídas esperadas
   - Troubleshooting

4. **[PIPELINE_AUTOMATIZADA_QAOA.md](PIPELINE_AUTOMATIZADA_QAOA.md)** (15 min)
   - Detalhes técnicos
   - Exemplos de execução
   - Modificações possíveis

**Total**: ~45 minutos para entender tudo completamente

---

## 📊 Preciso dos dados

### Dados Brutos
```
resultados_qaoa_otimizado/resultados_20251228_145335.csv
```

### Consolidado
```
auditoria_qaoa/auditoria_qaoa_master.csv
```

### Gráficos Interativos
```
auditoria_qaoa/auditoria_qaoa_energia.html    ← Abra no navegador!
auditoria_qaoa/auditoria_qaoa_tempo.html      ← Abra no navegador!
```

### Rastreabilidade
```
auditoria_qaoa/manifest_codigo.json           ← Hashes SHA-256
```

---

## 🎓 Para Artigo QUALIS A1

### Seção Methods
→ Cite: [PIPELINE_AUTOMATIZADA_QAOA.md](PIPELINE_AUTOMATIZADA_QAOA.md)

### Material Suplementar
→ Inclua:
- `auditoria_qaoa/auditoria_qaoa_master.csv` (Tabela S1)
- `auditoria_qaoa/README_AUDITORIA_QAOA.md` (Methodology)
- `auditoria_qaoa/manifest_codigo.json` (Code Reproducibility)

### Figuras
→ Use: `auditoria_qaoa/auditoria_qaoa_energia.html` (interativo)

Veja: [RESUMO_EXECUTIVO_PIPELINE_QAOA.md](RESUMO_EXECUTIVO_PIPELINE_QAOA.md)

---

## 🐛 Algo deu errado?

### Pipeline não executa
1. Ativar venv: `& ".\.venv\Scripts\Activate.ps1"`
2. Verificar deps: `pip list | grep qiskit`
3. Consultar: [PLAYBOOK_REPRODUCAO_QAOA.md](PLAYBOOK_REPRODUCAO_QAOA.md) → Troubleshooting

### Hashes divergem
- Código foi modificado
- Execute `python calculador_hashes_qaoa.py` para novo manifest

### Gráficos não aparecem
- Normal em headless
- Verifique arquivos `.html` no navegador

---

## 📋 Checklist Rápido

- [ ] Tenho Python 3.13+ com Qiskit? → [PLAYBOOK_REPRODUCAO_QAOA.md](PLAYBOOK_REPRODUCAO_QAOA.md#pré-requisitos)
- [ ] Executei `test_pipeline_verificacao.py`? → Resultado deve ser 5/5 ✅
- [ ] Abri `auditoria_qaoa/auditoria_qaoa_energia.html`? → Veja os gráficos!
- [ ] Li metodologia? → [auditoria_qaoa/README_AUDITORIA_QAOA.md](auditoria_qaoa/README_AUDITORIA_QAOA.md)
- [ ] Pronto para publicar? → [RESUMO_EXECUTIVO_PIPELINE_QAOA.md](RESUMO_EXECUTIVO_PIPELINE_QAOA.md)

---

## 🎯 Mapa de Documentação

```
COMEÇAR AQUI
    ↓
[QUICKSTART_PIPELINE_QAOA.md] (2 min, TL;DR)
    ↓
[VISAO_GERAL_PIPELINE.md] (5 min, entender arquitetura)
    ↓
┌─────────────────────┴──────────────────────┐
│                                            │
Preciso EXECUTAR                    Preciso ESTUDAR
    │                                    │
    ↓                                    ↓
[PLAYBOOK_REPRODUCAO_QAOA.md]  [PIPELINE_AUTOMATIZADA_QAOA.md]
(passo a passo)                (detalhes técnicos)
    ↓                                    ↓
Run: python experimento...    Modificar/customizar
    ↓                                    ↓
Ver gráficos                    Entender código
    ↓                                    ↓
Validar: QUALIS A1            Design decisions
    ↓                                    ↓
PUBLICAR                       CONTRIBUIR
```

---

## 🔗 Links Principais

| Preciso de... | Ler... |
| --- | --- |
| Começar em 30 seg | [QUICKSTART_PIPELINE_QAOA.md](QUICKSTART_PIPELINE_QAOA.md) |
| Entender a arquitetura | [VISAO_GERAL_PIPELINE.md](VISAO_GERAL_PIPELINE.md) |
| Reproduzir manualmente | [PLAYBOOK_REPRODUCAO_QAOA.md](PLAYBOOK_REPRODUCAO_QAOA.md) |
| Detalhes técnicos | [PIPELINE_AUTOMATIZADA_QAOA.md](PIPELINE_AUTOMATIZADA_QAOA.md) |
| Status geral e checklist | [RESUMO_EXECUTIVO_PIPELINE_QAOA.md](RESUMO_EXECUTIVO_PIPELINE_QAOA.md) |
| Inventário completo | [INVENTARIO_PIPELINE_QAOA.md](INVENTARIO_PIPELINE_QAOA.md) |
| Metodologia de auditoria | [auditoria_qaoa/README_AUDITORIA_QAOA.md](auditoria_qaoa/README_AUDITORIA_QAOA.md) |
| Resumo da sessão | [SUMARIO_SESSAO_PIPELINE.md](SUMARIO_SESSAO_PIPELINE.md) |

---

## 🚀 Próximos Passos

### Imediato (agora mesmo)
```powershell
python test_pipeline_verificacao.py
```

### Curto Prazo (próximas horas)
1. Ler [VISAO_GERAL_PIPELINE.md](VISAO_GERAL_PIPELINE.md)
2. Executar `python experimento_qaoa_otimizado.py`
3. Abrir gráficos no navegador

### Médio Prazo (para publicação)
1. Integrar com artigo QUALIS A1
2. Incluir material suplementar
3. Descrever em Methods

---

## 📞 Suporte Rápido

### P: Onde estão os resultados?
**R**: `resultados_qaoa_otimizado/` e `auditoria_qaoa/`

### P: Como verificar se tudo funciona?
**R**: `python test_pipeline_verificacao.py` → esperado 5/5 ✅

### P: Como ver os gráficos?
**R**: Abra `auditoria_qaoa/auditoria_qaoa_energia.html` no navegador

### P: Como reproduzir?
**R**: Siga [PLAYBOOK_REPRODUCAO_QAOA.md](PLAYBOOK_REPRODUCAO_QAOA.md) (7 passos)

### P: Está pronto para publicar?
**R**: Sim! Veja [RESUMO_EXECUTIVO_PIPELINE_QAOA.md](RESUMO_EXECUTIVO_PIPELINE_QAOA.md)

---

## ✨ Resumão (30 segundos)

✅ Pipeline com 5 etapas automatizadas  
✅ Teste rápido: `python test_pipeline_verificacao.py` (5/5 PASS)  
✅ Dados em: `auditoria_qaoa/auditoria_qaoa_master.csv`  
✅ Gráficos: `auditoria_qaoa/auditoria_qaoa_energia.html`  
✅ Rastreabilidade: `auditoria_qaoa/manifest_codigo.json` (SHA-256)  
✅ Conforme QUALIS A1  
✅ Pronto para publicar  

**Comece com**: [QUICKSTART_PIPELINE_QAOA.md](QUICKSTART_PIPELINE_QAOA.md)

---

**Desenvolvido em**: 28/12/2025  
**Status**: ✅ Operacional  
**Próximo**: Publicação QUALIS A1
