#!/usr/bin/env python3
"""
Monitor de Progresso do Framework v7.2
Verifica o progresso da execução e notifica em marcos importantes.
"""
import time
from datetime import datetime
from pathlib import Path


def extrair_progresso_terminal(pasta_resultados: str) -> dict:
    """Extrai progresso analisando os CSVs individuais gerados."""
    pasta_exp = Path(pasta_resultados) / "experimentos_individuais"
    
    if not pasta_exp.exists():
        return {
            'experimentos_completos': 0,
            'total_esperado': 8280,
            'percentual': 0.0,
            'status': 'Aguardando início'
        }
    
    # Contar CSVs gerados
    csvs = list(pasta_exp.glob("exp_*.csv"))
    n_completos = len(csvs)
    total = 8280
    percentual = (n_completos / total) * 100
    
    return {
        'experimentos_completos': n_completos,
        'total_esperado': total,
        'percentual': percentual,
        'status': 'Em execução' if n_completos < total else 'Concluído'
    }


def formatar_tempo_estimado(experimentos_restantes: int, tempo_medio_por_exp: float) -> str:
    """Formata tempo estimado para conclusão."""
    segundos_restantes = experimentos_restantes * tempo_medio_por_exp
    horas = segundos_restantes / 3600
    dias = horas / 24
    
    if dias >= 1:
        return f"{dias:.1f} dias"
    elif horas >= 1:
        return f"{horas:.1f} horas"
    else:
        return f"{segundos_restantes/60:.0f} minutos"


def monitorar(pasta_resultados: str, intervalo_segundos: int = 300):
    """
    Monitora progresso e exibe notificações em marcos.
    
    Args:
        pasta_resultados: Caminho para a pasta de resultados
        intervalo_segundos: Intervalo entre verificações (padrão: 5 minutos)
    """
    marcos_notificados = set()
    marcos = [25, 50, 75, 100]
    tempo_medio_exp = 150  # segundos por experimento (estimativa conservadora)
    
    print("="*80)
    print(" MONITOR DE PROGRESSO - FRAMEWORK v7.2")
    print("="*80)
    print(f"Pasta de resultados: {pasta_resultados}")
    print(f"Intervalo de verificação: {intervalo_segundos}s ({intervalo_segundos/60:.1f} min)")
    print(f"Marcos de notificação: {', '.join(f'{m}%' for m in marcos)}")
    print("="*80)
    print("\nMonitoramento iniciado. Pressione Ctrl+C para encerrar.\n")
    
    inicio_monitoramento = time.time()
    info = {'percentual': 0.0, 'experimentos_completos': 0, 'total_esperado': 8280, 'status': 'Inicializando'}
    
    try:
        while True:
            # Buscar pasta de resultados mais recente se não especificada
            if not pasta_resultados or not Path(pasta_resultados).exists():
                pasta_base = Path.cwd()
                pastas = sorted(pasta_base.glob("resultados_*"), reverse=True)
                if pastas:
                    pasta_resultados = str(pastas[0])
            
            # Extrair progresso
            info = extrair_progresso_terminal(pasta_resultados)
            
            # Calcular tempo estimado
            restantes = info['total_esperado'] - info['experimentos_completos']
            tempo_estimado = formatar_tempo_estimado(restantes, tempo_medio_exp)
            
            # Timestamp
            agora = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            
            # Exibir status
            print(f"[{agora}] Status: {info['status']}")
            print(f"  Progresso: {info['experimentos_completos']}/{info['total_esperado']} "
                  f"({info['percentual']:.2f}%)")
            print(f"  Restantes: {restantes} experimentos")
            print(f"  Tempo estimado: {tempo_estimado}")
            
            # Verificar marcos
            percentual_atual = int(info['percentual'])
            for marco in marcos:
                if percentual_atual >= marco and marco not in marcos_notificados:
                    marcos_notificados.add(marco)
                    print("\n" + "="*80)
                    print(f" 🎯 MARCO ATINGIDO: {marco}% COMPLETO!")
                    print("="*80)
                    print(f"  Experimentos: {info['experimentos_completos']}/{info['total_esperado']}")
                    print(f"  Tempo decorrido: {(time.time() - inicio_monitoramento)/3600:.1f}h")
                    print(f"  Tempo estimado restante: {tempo_estimado}")
                    print("="*80 + "\n")
            
            # Finalizado?
            if info['percentual'] >= 100:
                print("\n" + "="*80)
                print(" ✅ EXECUÇÃO COMPLETA FINALIZADA COM SUCESSO!")
                print("="*80)
                print(f"  Total de experimentos: {info['experimentos_completos']}")
                print(f"  Tempo total: {(time.time() - inicio_monitoramento)/3600:.1f}h")
                print(f"  Pasta de resultados: {pasta_resultados}")
                print("="*80)
                print("\n📊 Próximas etapas automáticas:")
                print("  1. Consolidação de CSVs → resultados_completos_artigo.csv")
                print("  2. Comparação de baselines → comparacao_baselines.csv")
                print("  3. Metadados → metadata_orchestrator.json")
                print("  4. Análises estatísticas")
                print("  5. Visualizações interativas (9 figuras HTML)")
                print("  6. Resumo final com melhores configurações")
                print("\nVerifique os logs do terminal principal para detalhes.")
                break
            
            # Aguardar próximo ciclo
            print(f"  Próxima verificação em {intervalo_segundos/60:.1f} minutos...\n")
            time.sleep(intervalo_segundos)
    
    except KeyboardInterrupt:
        print("\n\n⚠️  Monitoramento interrompido pelo usuário.")
        print(f"Último status: {info['percentual']:.2f}% completo")
        print("A execução do framework continua em background no terminal principal.\n")


if __name__ == "__main__":
    import sys
    
    # Argumentos: pasta_resultados (opcional) e intervalo (opcional)
    pasta = sys.argv[1] if len(sys.argv) > 1 else ""
    intervalo = int(sys.argv[2]) if len(sys.argv) > 2 else 300  # 5 minutos padrão
    
    monitorar(pasta, intervalo)
