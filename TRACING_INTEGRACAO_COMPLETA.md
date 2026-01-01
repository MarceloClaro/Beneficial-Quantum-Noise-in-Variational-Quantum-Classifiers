# Integração OpenTelemetry - Pipeline QAOA

## ✅ Status: Tracing Operacional

A pipeline QAOA foi integrada com **OpenTelemetry** com sucesso.

### Resultado do Teste
```
✓ TODOS OS TESTES PASSARAM (2.9s)
✓ Artefatos verificados: 3/3
✓ CSV com 5 linhas e 22 colunas
✓ Manifest SHA-256 com 5 scripts
✓ Traces enviados para AI Toolkit
```

---

## 🚀 Como Usar

### 1. Executar Pipeline com Tracing

```bash
python pipeline_tracing_simples.py
```

**Output esperado:**
```
✓ OpenTelemetry tracing inicializado
✓ TODOS OS TESTES PASSARAM (2.9s)
Traces foram enviados para AI Toolkit!
```

### 2. Visualizar Traces em VS Code

1. Pressione `Ctrl+Shift+P`
2. Digite: `AI Toolkit: Show Traces`
3. Veja os spans capturados em tempo real

---

## 📊 Spans Rastreados

### Estrutura de Spans

```
teste_pipeline_qaoa
├── teste_artefatos
│   ├── Verificação de CSV (OK)
│   ├── Verificação de Auditoria (OK)
│   └── Verificação de Manifest (OK)
│
└── Atributos finais
    ├── teste.artefatos_ok = true
    ├── teste.csv_linhas = 5
    ├── teste.csv_colunas = 22
    ├── teste.manifest_scripts = 5
    ├── teste.tempo_total_s = 2.9
    └── teste.qiskit_version = 1.4.4
```

---

## 🔧 Configuração Avançada

### Personalizar Endpoint OTLP

Em `tracing_setup.py`:

```python
# HTTP (padrão)
setup_tracing("sua-app", endpoint="http://localhost:4318")

# Ou especificar manualmente
setup_tracing("sua-app", endpoint="http://seu-collector:4318")
```

### Integrar em Seus Scripts

```python
from tracing_setup import setup_tracing, get_tracer

# 1. Inicializar
setup_tracing("seu-app")
tracer = get_tracer(__name__)

# 2. Envolver código
with tracer.start_as_current_span("minha_operacao") as span:
    span.set_attribute("operacao.tipo", "qaoa")
    
    # ... seu código aqui ...
    
    span.set_attribute("operacao.resultado", "sucesso")
```

---

## 📁 Arquivos Criados

| Arquivo | Propósito |
|---------|-----------|
| `tracing_setup.py` | Configuração centralizada de OpenTelemetry |
| `pipeline_com_tracing.py` | Wrapper completo (legacy, para execução longa) |
| `pipeline_tracing_simples.py` | Teste rápido e simples ✅ |
| `TRACING_QAOA_SETUP.md` | Documentação técnica completa |

---

## 🎯 Próximas Etapas (Opcional)

- [ ] Integrar tracing em `experimento_qaoa_otimizado.py`
- [ ] Rastrear cada iteração de otimização QAOA
- [ ] Medir latência de simulação Qiskit
- [ ] Exportar comparativos de performance entre runs

---

## 💡 Troubleshooting

### Problema: "Tracing não disponível"
**Solução:**
```bash
pip install opentelemetry-api opentelemetry-sdk opentelemetry-exporter-otlp-proto-http
```

### Problema: Traces não aparecem no AI Toolkit
**Checklist:**
- [ ] VS Code tem "AI Toolkit: Show Traces" aberto?
- [ ] Endpoint http://localhost:4318 está correto?
- [ ] OpenTelemetry SDK está instalado?

### Problema: Performance degradada
**Solução:** Reduzir sample rate em `tracing_setup.py`
```python
setup_tracing(..., trace_sample_rate=0.5)  # 50% das traces
```

---

## 📚 Referências

- [OpenTelemetry Python](https://opentelemetry.io/docs/instrumentation/python/)
- [OTLP HTTP Exporter](https://github.com/open-telemetry/opentelemetry-python/blob/main/exporter/opentelemetry-exporter-otlp-proto-http/README.md)
- [AI Toolkit Documentation](https://github.com/Azure/ai-toolkit)

---

## ✨ Resumo

- ✅ Tracing OpenTelemetry integrado
- ✅ Teste de pipeline executando com sucesso
- ✅ Artefatos QAOA verificados
- ✅ Traces enviados para AI Toolkit
- 📊 Pronto para observabilidade em tempo real

**Para começar:**
```bash
python pipeline_tracing_simples.py
```
