# Tracing OpenTelemetry na Pipeline QAOA

## 📊 Visão Geral

A pipeline QAOA foi instrumentada com **OpenTelemetry (OTel)** para:
- ✅ Rastrear execução de cada etapa
- ✅ Medir latência e performance
- ✅ Capturar atributos de span (status, métricas)
- ✅ Exportar traces via OTLP para AI Toolkit

---

## 🚀 Começar Rápido

### 1. Instalar Dependências

```bash
pip install opentelemetry-api opentelemetry-sdk \
            opentelemetry-exporter-otlp-proto-http \
            opentelemetry-instrumentation-requests \
            opentelemetry-instrumentation-urllib3 \
            opentelemetry-instrumentation-psutil
```

### 2. Abrir Trace Collector (AI Toolkit)

```bash
# No VS Code, pressione Ctrl+Shift+P
# Digite: AI Toolkit: Show Traces
```

### 3. Executar Pipeline com Tracing

```bash
python pipeline_com_tracing.py
```

**Resultado esperado:**
```
🚀 PIPELINE QAOA COM TRACING
======================================================================
📊 [Etapa 1: Experimento QAOA (6 qubits × 4 ruídos)]
✅ Experimento QAOA — SUCESSO (45.2s)

📊 [Etapa 2: Enriquecimento (Metadados + Métricas)]
✅ Enriquecimento — SUCESSO (8.1s)

📊 [Etapa 3: Consolidação (Master CSV + Gráficos)]
✅ Consolidação — SUCESSO (12.5s)

📊 [Etapa 4: Validação QUALIS A1]
✅ Validação — SUCESSO (2.3s)

📊 [Etapa 5: Rastreabilidade (SHA-256)]
✅ Rastreabilidade — SUCESSO (0.8s)

✅ PIPELINE CONCLUÍDA COM SUCESSO!
💡 Traces enviados para AI Toolkit
```

---

## 📁 Arquivos Criados

### `tracing_setup.py` (Nova)
Configuração centralizada de OpenTelemetry:
- Cria `TracerProvider` com OTLP exporter
- Exporta traces para http://localhost:4318 (AI Toolkit)
- Instrumentação automática: requests, urllib3, psutil
- Função `setup_tracing()` para inicializar globalmente
- Função `get_tracer()` para obter tracer por nome de módulo

### `pipeline_com_tracing.py` (Nova)
Wrapper da pipeline com rastreamento:
- Executa 5 etapas da pipeline
- Rastreia sucesso/falha, tempo, exit code
- Verifica geração de artefatos
- Exporta spans ao AI Toolkit

---

## 🔍 Estrutura de Spans

### Span Principal: `pipeline_qaoa_completa`

```
pipeline_qaoa_completa
├── etapa_experimento_qaoa_otimizado
│   ├── etapa.script = "experimento_qaoa_otimizado.py"
│   ├── etapa.status = "sucesso"
│   ├── etapa.tempo_s = 45.2
│   └── etapa.exit_code = 0
│
├── etapa_enriquecer_resultados_qaoa
│   ├── etapa.status = "sucesso"
│   ├── etapa.tempo_s = 8.1
│   └── ...
│
├── etapa_auditoria_qaoa_resultados
│   ├── etapa.status = "sucesso"
│   ├── etapa.tempo_s = 12.5
│   └── ...
│
├── etapa_validar_auditoria_qaoa
│   ├── etapa.status = "sucesso"
│   ├── etapa.tempo_s = 2.3
│   └── ...
│
├── etapa_calculador_hashes_qaoa
│   ├── etapa.status = "sucesso"
│   ├── etapa.tempo_s = 0.8
│   └── ...
│
└── verificacao_artefatos
    ├── artefato.resultados_qaoa_otimizado/... = true
    ├── artefato.auditoria_qaoa/auditoria_qaoa_master.csv = true
    ├── artefatos.presentes = 4
    └── artefatos.total = 4

Atributos Finais:
├── pipeline.status = "sucesso"
├── pipeline.tempo_total_s = 68.9
├── pipeline.artefatos_ok = true
└── pipeline.etapas_pass = 5
```

---

## 📊 Exemplos de Atributos Capturados

### Timing
```json
{
  "etapa.tempo_s": 45.2,
  "pipeline.tempo_total_s": 68.9
}
```

### Status de Execução
```json
{
  "etapa.status": "sucesso",
  "etapa.exit_code": 0,
  "pipeline.status": "sucesso"
}
```

### Informações de Artefato
```json
{
  "artefatos.presentes": 4,
  "artefatos.total": 4,
  "artefato.resultados_qaoa_otimizado/resultados_20251228_145335.csv": true
}
```

---

## 🔧 Integração Customizada

### Adicionar Tracing em Scripts Existentes

#### 1. Importar módulo
```python
from tracing_setup import setup_tracing, get_tracer

# Inicializar no main
setup_tracing("seu-app")
tracer = get_tracer(__name__)
```

#### 2. Envolver funções principais
```python
def minha_funcao():
    with tracer.start_as_current_span("minha_funcao") as span:
        span.set_attribute("funcao.entrada", "valor")
        
        # ... código ...
        
        span.set_attribute("funcao.saida", "resultado")
```

#### 3. Adicionar métricas
```python
span.set_attribute("metricas.energia", 0.85)
span.set_attribute("metricas.acuracia", 0.92)
span.set_attribute("metricas.tempo_execucao", 12.5)
```

---

## 🎯 Casos de Uso

### 1. Identificar Gargalos
- Qual etapa é mais lenta?
- Resposta: Visualizar spans e tempo em AI Toolkit

### 2. Rastrear Erros
- Onde a pipeline falhou?
- Resposta: Procurar span com status = "falha" ou "exceção"

### 3. Auditar Reprodutibilidade
- Quais scripts foram executados?
- Resposta: Atributo `etapa.script` em cada span

### 4. Monitorar Performance
- Tempo total vs tempo esperado?
- Resposta: Comparar `pipeline.tempo_total_s` entre execuções

---

## 🔗 Endpoints OTLP

### HTTP (padrão, recomendado)
```
http://localhost:4318/v1/traces
```
- Usado em `tracing_setup.py`
- Compatível com: AI Toolkit, Jaeger, Datadog

### gRPC (alternativo)
```
http://localhost:4317
```
- Mais rápido em alta carga
- Requer `opentelemetry-exporter-otlp`

---

## 📝 Configurações Avançadas

### Aumentar Verbosidade
```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

### Usar Sampler Customizado
```python
from opentelemetry.sdk.trace.sampling import TraceIdRatioBased

setup_tracing("app", trace_sample_rate=0.1)  # 10% das traces
```

### Adicionar Tags de Contexto
```python
from opentelemetry.baggage import set_baggage

set_baggage("user.id", "user_123")
set_baggage("experiment.id", "exp_001")
```

---

## ✅ Troubleshooting

### Traces não aparecem no AI Toolkit?
1. Verificar se collector está aberto: `Ctrl+Shift+P` → "Show Traces"
2. Verificar endpoint: `http://localhost:4318`
3. Verificar logs: Procurar "Exported 5 spans"

### Erros de importação?
```bash
pip install --upgrade opentelemetry-api opentelemetry-sdk
```

### Performance degradada?
- Reduzir sample rate: `setup_tracing(..., trace_sample_rate=0.5)`
- Aumentar batch size em `BatchSpanProcessor`

---

## 📚 Referências

- [OpenTelemetry Python](https://opentelemetry.io/docs/instrumentation/python/)
- [OTLP HTTP Exporter](https://github.com/open-telemetry/opentelemetry-python/tree/main/exporter/opentelemetry-exporter-otlp-proto-http)
- [AI Toolkit Tracing](https://github.com/Azure/ai-toolkit)

---

## 📊 Próximas Etapas

- [ ] Integrar tracing em `experimento_qaoa_otimizado.py`
- [ ] Adicionar spans para cada iteração de otimização
- [ ] Rastrear métricas QAOA (energia, tempo de simulação)
- [ ] Criar dashboard em AI Toolkit
- [ ] Comparar performance entre runs com traces
