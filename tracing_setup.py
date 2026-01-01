"""
tracing_setup.py
Configuração de OpenTelemetry Tracing para Pipeline QAOA
Rastreia execução de experimentos e consolidação de dados
"""

from opentelemetry import trace, metrics
from opentelemetry.exporter.otlp.proto.http.trace_exporter import OTLPSpanExporter
from opentelemetry.exporter.otlp.proto.http.metric_exporter import OTLPMetricExporter
from opentelemetry.sdk.trace import TracerProvider
from opentelemetry.sdk.trace.export import BatchSpanProcessor
from opentelemetry.sdk.metrics import MeterProvider
from opentelemetry.sdk.metrics.export import PeriodicExportingMetricReader
from opentelemetry.sdk.resources import Resource
import logging

# Instrumentação opcional - comentada para evitar dependências extras
# from opentelemetry.instrumentation.requests import RequestsInstrumentor
# from opentelemetry.instrumentation.urllib3 import URLLib3Instrumentor

# Configurar logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def setup_tracing(service_name: str = "qaoa-pipeline", endpoint: str = "http://localhost:4318"):
    """
    Configurar OpenTelemetry Tracing para a pipeline QAOA

    Args:
        service_name: Nome do serviço para tracing
        endpoint: Endpoint OTLP (padrão: AI Toolkit local)

    Returns:
        tracer: Tracer configurado para uso
    """
    try:
        # Criar recurso com atributos da aplicação
        resource = Resource.create({
            "service.name": service_name,
            "service.version": "1.0.0",
        })

        # Configurar exportador OTLP (envia para AI Toolkit)
        otlp_exporter = OTLPSpanExporter(
            endpoint=f"{endpoint}/v1/traces",
        )

        # Criar TracerProvider com exportador
        tracer_provider = TracerProvider(resource=resource)
        tracer_provider.add_span_processor(BatchSpanProcessor(otlp_exporter))

        # Definir como tracer global
        trace.set_tracer_provider(tracer_provider)

        # Instrumentação automática comentada (requer pacotes extras)
        # RequestsInstrumentor().instrument()
        # URLLib3Instrumentor().instrument()

        # Configurar métricas (opcional)
        metric_reader = PeriodicExportingMetricReader(
            OTLPMetricExporter(
                endpoint=f"{endpoint}/v1/metrics",
            )
        )
        meter_provider = MeterProvider(resource=resource, metric_readers=[metric_reader])
        metrics.set_meter_provider(meter_provider)

        logger.info(f"✅ Tracing configurado: {service_name} → {endpoint}")
        return trace.get_tracer(__name__)

    except Exception as e:
        logger.warning(f"⚠️  Tracing não disponível: {e}")
        logger.warning("💡 Certifique-se que AI Toolkit trace collector está rodando")
        # Retornar no-op tracer
        return trace.get_tracer(__name__)


def get_tracer(module_name: str | None = None):
    """
    Obter tracer para um módulo específico

    Args:
        module_name: Nome do módulo (ex: "experimento_qaoa_otimizado")

    Returns:
        tracer: Tracer para o módulo
    """
    return trace.get_tracer(module_name or __name__)


# Exemplo de uso com context manager
if __name__ == "__main__":
    # Configurar tracing na inicialização
    tracer = setup_tracing("qaoa-pipeline-test")

    # Exemplo: Rastrear uma operação
    with tracer.start_as_current_span("teste_operacao") as span:
        span.set_attribute("operacao.tipo", "teste")
        span.set_attribute("operacao.status", "iniciado")
        logger.info("Operação de teste")
        # ... código aqui
        span.set_attribute("operacao.status", "concluído")

    logger.info("✅ Exemplo de tracing concluído")
    logger.info("💡 Abra http://localhost:4318 para visualizar traces")
