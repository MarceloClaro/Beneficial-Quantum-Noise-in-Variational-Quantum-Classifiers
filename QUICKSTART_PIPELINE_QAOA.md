# 🚀 Guia Rápido (30 segundos)

**Quer apenas usar a pipeline? Aqui está:**

## TL;DR (Resumido demais?)

```powershell
cd seu\projeto\qaoa
& ".\.venv\Scripts\Activate.ps1"
python test_pipeline_verificacao.py
```

**Resultado**: ✅ Se passar 5/5, tudo está funcionando!

---

## Um pouco mais de detalhe

### Opção A: Apenas Verificar

```powershell
python test_pipeline_verificacao.py
```

Leva ~30 segundos. Mostra status de todos os artefatos.

### Opção B: Rodar Experimento Novo

```powershell
# Experimento
python experimento_qaoa_otimizado.py       # ~2-3 min

# Consolidar tudo
python auditoria_qaoa_resultados.py        # ~1-2 min

# Validar
python validar_auditoria_qaoa.py           # ~10 seg

# Verificar rastreabilidade
python calculador_hashes_qaoa.py           # ~5 seg
```

**Total**: ~4-7 minutos. Gera todos os artefatos.

---

## O que cada script faz?

| Script | O que faz | Tempo |
| --- | --- | --- |
| `experimento_qaoa_otimizado.py` | Roda QAOA com 4 ruídos | 2-3 min |
| `enriquecer_resultados_qaoa.py` | Adiciona metadados + métricas | 30 seg |
| `auditoria_qaoa_resultados.py` | Consolida + gera gráficos | 1-2 min |
| `validar_auditoria_qaoa.py` | Valida QUALIS A1 | 10 seg |
| `calculador_hashes_qaoa.py` | SHA-256 dos scripts | 5 seg |
| `test_pipeline_verificacao.py` | Testa se tudo existe | 30 seg |

---

## Onde estão os resultados?

```
resultados_qaoa_otimizado/
└── resultados_YYYYMMDD_HHMMSS.csv    ← Seu experimento

auditoria_qaoa/
├── auditoria_qaoa_master.csv         ← Tabela consolidada
├── auditoria_qaoa_energia.html       ← Abra no navegador! 📊
├── auditoria_qaoa_tempo.html         ← Abra no navegador! ⏱️
├── manifest_codigo.json              ← Hashes SHA-256
└── README_AUDITORIA_QAOA.md         ← Metodologia
```

---

## Preciso de mais informações?

- **Entender o que a pipeline faz**: [VISAO_GERAL_PIPELINE.md](VISAO_GERAL_PIPELINE.md) (5 min de leitura)
- **Como reproduzir**: [PLAYBOOK_REPRODUCAO_QAOA.md](PLAYBOOK_REPRODUCAO_QAOA.md) (passo a passo)
- **Detalhes técnicos**: [PIPELINE_AUTOMATIZADA_QAOA.md](PIPELINE_AUTOMATIZADA_QAOA.md) (completo)
- **Status geral**: [RESUMO_EXECUTIVO_PIPELINE_QAOA.md](RESUMO_EXECUTIVO_PIPELINE_QAOA.md) (checklist)

---

## Está tudo funcionando?

```powershell
python test_pipeline_verificacao.py
```

Se retornar:
```
✅ PIPELINE FUNCIONANDO CORRETAMENTE!
```

→ Então sim! 🎉

Se retornar:
```
⚠️  PIPELINE INCOMPLETA
```

→ Então execute `python experimento_qaoa_otimizado.py` para gerar os dados que faltam.

---

## Posso modificar?

Claro! Edite `experimento_qaoa_otimizado.py`:

```python
CONFIG = ConfigQAOA(
    n_qubits=8,          # ← mude de 6 para 8 para mais qubits
    p_layers=1,
    shots=1024,          # ← ou aumente shots
    max_iter=100,
    seed=42
)
```

Depois re-execute a pipeline. Tudo será consolidado automaticamente.

---

## Pronto!

Enjoy! 🚀✨
