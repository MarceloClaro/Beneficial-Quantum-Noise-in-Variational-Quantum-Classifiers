import os
import sys
import json
from typing import Dict, Any, Optional
import pandas as pd

def consolidate_results(results_dir: str, verbose: bool = True) -> Dict[str, Any]:
    """Consolida CSVs individuais em um único arquivo e gera comparações de baseline.

    Parâmetros:
    - results_dir: caminho absoluto/relativo para a pasta de resultados (contendo 'experimentos_individuais').
    - verbose: se True, imprime logs amigáveis.

    Retorna:
    - dict com caminhos gerados e estatísticas básicas.
    """
    results_dir = os.path.abspath(results_dir)
    individual_dir = os.path.join(results_dir, 'experimentos_individuais')

    if not os.path.isdir(individual_dir):
        msg = f"Pasta de experimentos individuais não encontrada: {individual_dir}"
        if verbose:
            print('ERRO:', msg)
        raise FileNotFoundError(msg)

    csv_files = sorted([os.path.join(individual_dir, f) for f in os.listdir(individual_dir) if f.endswith('.csv')])
    if len(csv_files) == 0:
        msg = f"Nenhum CSV individual encontrado em {individual_dir}"
        if verbose:
            print(msg)
        raise FileNotFoundError(msg)

    if verbose:
        print(f'Encontrados {len(csv_files)} CSVs individuais. Consolidando...')

    dfs = []
    for p in csv_files:
        try:
            df = pd.read_csv(p)
            dfs.append(df)
        except Exception as e:
            if verbose:
                print('Falha ao ler', p, str(e)[:200])

    if len(dfs) == 0:
        msg = 'Nenhum CSV lido com sucesso.'
        if verbose:
            print(msg)
        raise RuntimeError(msg)

    df_all = pd.concat(dfs, ignore_index=True)
    consolidated_path = os.path.join(results_dir, 'resultados_completos_artigo.csv')
    df_all.to_csv(consolidated_path, index=False)
    if verbose:
        print('CSV consolidado salvo em:', consolidated_path)
        print('Linhas no consolidado:', len(df_all))
        print('Colunas:', list(df_all.columns))
        try:
            print('Amostra (5 linhas):')
            print(df_all.head().to_string())
        except Exception:
            pass

    # Gerar comparacao_baselines.csv se colunas existirem
    comp_path: Optional[str] = None
    required = {'dataset', 'arquitetura', 'tipo_ruido', 'nivel_ruido', 'acuracia_teste'}
    if required.issubset(df_all.columns):
        try:
            vqc_melhor = df_all[df_all['tipo_ruido'] != 'classico'].groupby('dataset')['acuracia_teste'].max()
            vqc_sem_ruido = df_all[(df_all['tipo_ruido'].isin(['sem_ruido', 'none'])) | (df_all['nivel_ruido'] == 0.0)].groupby('dataset')['acuracia_teste'].mean()
            svm = df_all[df_all['arquitetura'].str.lower().str.contains('svm', na=False)].groupby('dataset')['acuracia_teste'].mean()
            rf = df_all[df_all['arquitetura'].str.lower().str.contains('randomforest|random_forest|rf', na=False)].groupby('dataset')['acuracia_teste'].mean()

            comparacao = pd.DataFrame({
                'dataset': vqc_melhor.index,
                'vqc_melhor': vqc_melhor.values,
                'vqc_sem_ruido_media': vqc_sem_ruido.reindex(vqc_melhor.index).values,
                'svm': svm.reindex(vqc_melhor.index).values,
                'rf': rf.reindex(vqc_melhor.index).values
            })
            comparacao['delta_vqc_svm'] = comparacao['vqc_melhor'] - comparacao['svm']
            comparacao['delta_vqc_rf'] = comparacao['vqc_melhor'] - comparacao['rf']
            comp_path = os.path.join(results_dir, 'comparacao_baselines.csv')
            comparacao.to_csv(comp_path, index=False)
            if verbose:
                print('comparacao_baselines.csv salvo em:', comp_path)
        except Exception as e:
            if verbose:
                print('Falha ao gerar comparacao_baselines.csv:', str(e))
    else:
        if verbose:
            print('Colunas necessárias para gerar comparacao_baselines.csv não encontradas.')

    # Retornar sumário útil (para orquestrador)
    return {
        'results_dir': results_dir,
        'num_csvs_individuais': len(csv_files),
        'consolidated_path': consolidated_path,
        'comparacao_path': comp_path,
        'rows_consolidated': len(df_all),
        'columns': list(df_all.columns)
    }


def _main():
    import argparse
    parser = argparse.ArgumentParser(description='Consolida CSVs de experimentos individuais em um único arquivo.')
    parser.add_argument('results_dir', nargs='?', default=None, help='Diretório de resultados (contendo "experimentos_individuais/")')
    parser.add_argument('--json', dest='json_out', default=None, help='Se informado, salva um JSON com o sumário neste caminho')
    args = parser.parse_args()

    ROOT = os.path.dirname(os.path.abspath(__file__))
    # Se não informado, tenta subir um nível e detectar um resultados_* mais recente (fallback)
    results_dir = args.results_dir
    if results_dir is None:
        parent = os.path.abspath(os.path.join(ROOT, '..'))
        candidatos = [os.path.join(parent, d) for d in os.listdir(parent) if d.startswith('resultados_') and os.path.isdir(os.path.join(parent, d))]
        results_dir = max(candidatos, key=os.path.getmtime) if candidatos else None
        if results_dir is None:
            print('ERRO: Nenhum diretório resultados_* encontrado automaticamente. Informe results_dir.')
            sys.exit(1)

    try:
        summary = consolidate_results(results_dir, verbose=True)
        if args.json_out:
            with open(args.json_out, 'w', encoding='utf-8') as f:
                json.dump(summary, f, indent=2, ensure_ascii=False)
    except Exception as e:
        print('Falha na consolidação:', str(e))
        sys.exit(1)


if __name__ == '__main__':
    _main()

