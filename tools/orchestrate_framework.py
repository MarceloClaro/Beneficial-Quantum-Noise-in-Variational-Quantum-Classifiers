import os
import sys
import json
import time
import shutil
import argparse
import subprocess
from typing import Optional, Dict, Any

# Permite import relativo do consolidator
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..'))
CONSOLIDATOR_PATH = os.path.join(SCRIPT_DIR, 'consolidate_results.py')
FRAMEWORK_PATH = os.path.join(REPO_ROOT, 'framework_investigativo_completo.py')

# Import dinâmico do consolidator como módulo
try:
    sys.path.insert(0, SCRIPT_DIR)
    from consolidate_results import consolidate_results
except Exception:
    consolidate_results = None


def _list_result_dirs() -> list:
    return [os.path.join(REPO_ROOT, d) for d in os.listdir(REPO_ROOT)
            if d.startswith('resultados_') and os.path.isdir(os.path.join(REPO_ROOT, d))]


def _detect_new_results_dir(before: set) -> Optional[str]:
    after = set(_list_result_dirs())
    new_dirs = list(after - before)
    if not new_dirs:
        return None
    # Se houver mais de um, escolhe o mais recente por mtime
    new_dirs.sort(key=lambda p: os.path.getmtime(p), reverse=True)
    return new_dirs[0]


def _write_metadata(results_dir: str, extra: Optional[Dict[str, Any]] = None) -> str:
    """Gera um metadata_orchestrator.json com sumário de arquivos e estatísticas básicas."""
    results_dir = os.path.abspath(results_dir)
    meta = {
        'tipo': 'metadata_orchestrator',
        'versao': '1.0',
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'arquivos': sorted([f for f in os.listdir(results_dir) if os.path.isfile(os.path.join(results_dir, f))])
    }
    if extra:
        meta.update(extra)
    out_path = os.path.join(results_dir, 'metadata_orchestrator.json')
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    return out_path


def _run_framework(mode: str, quick: bool, trials: Optional[int], epocas_bayes: Optional[int], dataset_bayes: Optional[str], timeout: Optional[int]) -> Dict[str, Any]:
    """Executa o framework em subprocesso e retorna o diretório de resultados detectado e logs."""
    if not os.path.isfile(FRAMEWORK_PATH):
        raise FileNotFoundError(f'framework_investigativo_completo.py não encontrado em {FRAMEWORK_PATH}')

    env = os.environ.copy()
    if quick:
        env['VQC_QUICK'] = '1'

    args = [sys.executable, FRAMEWORK_PATH]
    if mode == 'bayes':
        args.append('--bayes')
    elif mode == 'both':
        args.append('--bayes-after-grid')

    # parâmetros bayes opcionais
    if trials is not None:
        args += ['--trials', str(trials)]
    if epocas_bayes is not None:
        args += ['--epocas-bayes', str(epocas_bayes)]
    if dataset_bayes is not None:
        args += ['--dataset-bayes', str(dataset_bayes)]

    # Detectar diretórios antes
    before = set(_list_result_dirs())

    # Executa e captura saída para um arquivo temporário
    tmp_log = os.path.join(SCRIPT_DIR, f'orchestrator_temp_{int(time.time())}.log')
    with open(tmp_log, 'w', encoding='utf-8') as logf:
        proc = subprocess.Popen(args, cwd=REPO_ROOT, env=env, stdout=logf, stderr=subprocess.STDOUT)
        try:
            proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            return {'status': 'timeout', 'log_path': tmp_log, 'results_dir': None}

    code = proc.returncode
    new_results_dir = _detect_new_results_dir(before)

    return {
        'status': 'ok' if code == 0 else f'exit_{code}',
        'code': code,
        'log_path': tmp_log,
        'results_dir': new_results_dir
    }


def run_post(results_dir: str) -> Dict[str, Any]:
    """Executa apenas a pós-processamento: consolidação e metadados."""
    if consolidate_results is None:
        # fallback: executar consolidator via subprocess
        args = [sys.executable, CONSOLIDATOR_PATH, results_dir]
        proc = subprocess.run(args, cwd=REPO_ROOT)
        if proc.returncode != 0:
            raise RuntimeError('Falha ao consolidar resultados via subprocess.')
        summary = {'results_dir': os.path.abspath(results_dir)}
    else:
        summary = consolidate_results(results_dir, verbose=True)

    # Metadados
    meta_path = _write_metadata(summary['results_dir'], extra={
        'consolidated_csv': summary.get('consolidated_path'),
        'comparacao_baselines': summary.get('comparacao_path'),
        'rows_consolidated': summary.get('rows_consolidated'),
        'columns': summary.get('columns')
    })
    summary['metadata_orchestrator'] = meta_path
    return summary


def main():
    parser = argparse.ArgumentParser(description='Orquestrador do Framework Investigativo (grid/bayes/both/post)')
    parser.add_argument('--mode', choices=['grid', 'bayes', 'both', 'post'], default='post', help='Modo de execução')
    parser.add_argument('--quick', action='store_true', help='Modo rápido (VQC_QUICK=1)')
    parser.add_argument('--results-dir', default=None, help='Diretório de resultados (necessário no modo post)')
    parser.add_argument('--trials', type=int, default=None, help='Trials para o modo bayesano')
    parser.add_argument('--epocas-bayes', type=int, default=None, help='Épocas por trial no modo bayesano')
    parser.add_argument('--dataset-bayes', type=str, default=None, help='Dataset para o modo bayesano (moons, circles, iris, etc.)')
    parser.add_argument('--timeout', type=int, default=None, help='Timeout (segundos) para execução do framework')
    args = parser.parse_args()

    if args.mode == 'post':
        if not args.results_dir:
            parser.error('--results-dir é obrigatório no modo post')
        summary = run_post(args.results_dir)
        print(json.dumps(summary, indent=2, ensure_ascii=False))
        return

    # Execução do framework
    outcome = _run_framework(
        mode=args.mode,
        quick=args.quick,
        trials=args.trials,
        epocas_bayes=args.epocas_bayes,
        dataset_bayes=args.dataset_bayes,
        timeout=args.timeout
    )

    if outcome['status'] == 'timeout':
        print('⚠️ Execução atingiu timeout. Logs em:', outcome['log_path'])
        return

    if outcome['results_dir'] is None:
        print('⚠️ Não foi possível detectar o diretório de resultados. Logs em:', outcome['log_path'])
        return

    # Copiar log para dentro do diretório de resultados
    try:
        final_log = os.path.join(outcome['results_dir'], 'orchestrator.log')
        shutil.copyfile(outcome['log_path'], final_log)
    except Exception:
        final_log = outcome['log_path']

    print('Resultado da execução do framework:', outcome['status'], '| dir:', outcome['results_dir'])

    # Pós-processamento
    post_summary = run_post(outcome['results_dir'])
    post_summary['orchestrator_log'] = final_log
    print(json.dumps(post_summary, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
