"""
Exemplo de Uso Programático do Framework v7.2
==============================================

Este script demonstra como usar o framework de forma programática,
incluindo execução personalizada, consolidação de resultados e
análise de dados.
"""

import os
import json
import pandas as pd
from framework_investigativo_completo import (
    consolidar_e_gerar_metadados,
    carregar_datasets,
)


def exemplo_1_consolidar_resultados_existentes():
    """Exemplo 1: Consolidar um diretório de resultados já existente."""
    print("=" * 80)
    print("EXEMPLO 1: Consolidar Resultados Existentes")
    print("=" * 80)

    # Substitua pelo caminho do seu diretório de resultados
    results_dir = "resultados_2025-10-28_17-09-33"

    if not os.path.isdir(results_dir):
        print(f"⚠️ Diretório não encontrado: {results_dir}")
        print("   Execute o framework primeiro para gerar resultados.")
        return

    # Consolidar tudo automaticamente
    resultado = consolidar_e_gerar_metadados(results_dir, verbose=True)

    # Acessar informações
    print(f"\n✅ Status: {resultado['consolidacao']['status']}")
    print(f"📦 CSVs processados: {resultado['consolidacao']['num_csvs_individuais']}")
    print(f"📊 Linhas consolidadas: {resultado['consolidacao']['rows_consolidated']}")
    print(f"💾 Arquivo: {resultado['consolidacao']['consolidated_path']}")

    return resultado


def exemplo_2_analisar_resultados_consolidados():
    """Exemplo 2: Analisar resultados consolidados com Pandas."""
    print()
    print("=" * 80)
    print("EXEMPLO 2: Analisar Resultados com Pandas")
    print("=" * 80)

    # Substitua pelo caminho do seu CSV consolidado
    csv_path = "resultados_2025-10-28_17-09-33/resultados_completos_artigo.csv"

    if not os.path.isfile(csv_path):
        print(f"⚠️ Arquivo não encontrado: {csv_path}")
        return

    # Carregar dados
    df = pd.read_csv(csv_path)
    print(f"\n📊 Carregados {len(df)} experimentos")
    print(f"📋 Colunas: {list(df.columns)}")

    # Análise 1: Melhor configuração global
    print("\n🏆 MELHOR CONFIGURAÇÃO GLOBAL:")
    melhor = df.loc[df['acuracia_teste'].idxmax()]
    print(f"   Dataset: {melhor['dataset']}")
    print(f"   Arquitetura: {melhor['arquitetura']}")
    print(f"   Ruído: {melhor['tipo_ruido']}")
    print(f"   Nível: {melhor['nivel_ruido']:.3f}")
    print(f"   Acurácia: {melhor['acuracia_teste']:.4f}")

    # Análise 2: Comparar ruído vs. sem ruído
    print("\n📈 IMPACTO DO RUÍDO:")
    sem_ruido = df[df['tipo_ruido'] == 'sem_ruido']['acuracia_teste'].mean()
    print(f"   Média sem ruído: {sem_ruido:.4f}")

    for tipo in ['depolarizante', 'amplitude_damping', 'phase_damping']:
        com_ruido = df[df['tipo_ruido'] == tipo]['acuracia_teste'].mean()
        if not pd.isna(com_ruido):
            delta = com_ruido - sem_ruido
            status = "✓ BENÉFICO" if delta > 0 else "✗ Prejudicial"
            print(f"   {tipo:20s}: {com_ruido:.4f} (Δ={delta:+.4f}) {status}")

    # Análise 3: Por dataset
    print("\n📊 DESEMPENHO POR DATASET:")
    for dataset in df['dataset'].unique():
        df_ds = df[df['dataset'] == dataset]
        media = df_ds['acuracia_teste'].mean()
        melhor_ds = df_ds['acuracia_teste'].max()
        print(f"   {dataset:15s}: média={media:.4f}, melhor={melhor_ds:.4f}")

    return df


def exemplo_3_comparar_com_baselines():
    """Exemplo 3: Comparar VQC com SVM e Random Forest."""
    print()
    print("=" * 80)
    print("EXEMPLO 3: Comparação com Baselines Clássicos")
    print("=" * 80)

    # Substitua pelo caminho do seu arquivo de comparação
    comp_path = "resultados_2025-10-28_17-09-33/comparacao_baselines.csv"

    if not os.path.isfile(comp_path):
        print(f"⚠️ Arquivo não encontrado: {comp_path}")
        print("   Execute consolidar_e_gerar_metadados() primeiro.")
        return

    # Carregar comparação
    comp = pd.read_csv(comp_path)
    print(f"\n📊 Comparação para {len(comp)} datasets")

    # Datasets onde VQC vence
    vence_svm = comp[comp['delta_vqc_svm'] > 0]
    vence_rf = comp[comp['delta_vqc_rf'] > 0]

    print(f"\n🏆 VQC vence SVM em {len(vence_svm)} de {len(comp)} datasets:")
    for _, row in vence_svm.iterrows():
        print(f"   {row['dataset']:15s}: VQC={row['vqc_melhor']:.4f}, SVM={row['svm']:.4f}, Δ={row['delta_vqc_svm']:+.4f}")

    print(f"\n🏆 VQC vence RF em {len(vence_rf)} de {len(comp)} datasets:")
    for _, row in vence_rf.iterrows():
        print(f"   {row['dataset']:15s}: VQC={row['vqc_melhor']:.4f}, RF={row['rf']:.4f}, Δ={row['delta_vqc_rf']:+.4f}")

    # Vantagem média
    media_delta_svm = comp['delta_vqc_svm'].mean()
    media_delta_rf = comp['delta_vqc_rf'].mean()
    print("\n📊 Vantagem média:")
    print(f"   vs. SVM: {media_delta_svm:+.4f}")
    print(f"   vs. RF:  {media_delta_rf:+.4f}")

    return comp


def exemplo_4_carregar_datasets():
    """Exemplo 4: Carregar e explorar os datasets usados no framework."""
    print()
    print("=" * 80)
    print("EXEMPLO 4: Explorar Datasets")
    print("=" * 80)

    # Carregar datasets
    datasets = carregar_datasets(seed=42)

    print(f"\n📦 Carregados {len(datasets)} datasets:")
    for nome, data in datasets.items():
        n_train = len(data['y_train'])
        n_test = len(data['y_test'])
        n_features = data['X_train'].shape[1]
        n_classes = len(set(data['y_train']))
        print(f"   {nome:15s}: {n_train} treino, {n_test} teste, {n_features} features, {n_classes} classes")

    return datasets


def exemplo_5_metadados():
    """Exemplo 5: Ler e interpretar metadados."""
    print("\n" + "=" * 80)
    print("EXEMPLO 5: Interpretar Metadados")
    print("=" * 80)

    meta_path = "resultados_2025-10-28_17-09-33/metadata_orchestrator.json"

    if not os.path.isfile(meta_path):
        print(f"⚠️ Arquivo não encontrado: {meta_path}")
        return

    # Carregar metadados
    with open(meta_path, 'r', encoding='utf-8') as f:
        meta = json.load(f)

    print("\n📋 Metadados da execução:")
    print(f"   Versão do framework: {meta.get('versao_framework', 'N/A')}")
    print(f"   Timestamp: {meta.get('timestamp', 'N/A')}")
    print(f"   Status consolidação: {meta.get('consolidacao', {}).get('status', 'N/A')}")
    print(f"   CSVs processados: {meta.get('consolidacao', {}).get('num_csvs_individuais', 0)}")
    print(f"   Linhas consolidadas: {meta.get('consolidacao', {}).get('rows_consolidated', 0)}")

    print(f"\n📁 Arquivos gerados ({len(meta.get('arquivos_raiz', []))}):")
    for arquivo in meta.get('arquivos_raiz', [])[:10]:  # Primeiros 10
        print(f"   - {arquivo}")

    print(f"\n📂 Subpastas ({len(meta.get('subpastas', []))}):")
    for pasta in meta.get('subpastas', []):
        print(f"   - {pasta}/")

    return meta


def main():
    """Função principal que executa todos os exemplos."""
    print("\n" + "=" * 80)
    print("EXEMPLOS DE USO PROGRAMÁTICO - Framework v7.2")
    print("=" * 80)

    # Execute os exemplos que desejar
    # Comente/descomente conforme necessário

    # exemplo_1_consolidar_resultados_existentes()
    # exemplo_2_analisar_resultados_consolidados()
    # exemplo_3_comparar_com_baselines()
    # exemplo_4_carregar_datasets()
    # exemplo_5_metadados()

    print("\n" + "=" * 80)
    print("✅ Exemplos disponíveis:")
    print("   1. Consolidar resultados existentes")
    print("   2. Analisar resultados com Pandas")
    print("   3. Comparar VQC com baselines clássicos")
    print("   4. Explorar datasets")
    print("   5. Interpretar metadados")
    print("\nDescomente as funções em main() para executar.")
    print("=" * 80)


if __name__ == "__main__":
    main()
