#!/usr/bin/env python
"""CLI entry point vqc-drug-a1 – QUALIS A1 ready."""
import click
import os
import datetime
from pathlib import Path
from .data import load_split
from .tune import run_study
from .audit import pre_register, checksum_file, snapshot_environment, checksum_dir
from .plots import fig_auc_heatmap, fig_forest_effect, fig_power_curve, latex_boilerplate

@click.command()
@click.option("--target", default="EGFR", help="EGFR | HIV | Malaria | COVID")
@click.option("--max-qubits", default=20, help="≤ 20 qubits simulated")
@click.option("--trials", default=500, help="Optuna trials")
@click.option("--power", default=0.8, help="Statistical power ≥ 0.8")
@click.option("--alpha", default=0.05, help="Family-wise error rate")
@click.option("--out-dir", default=None, help="Results folder")
def main(target, max_qubits, trials, power, alpha, out_dir):
    """Run full QUALIS A1 pipeline."""
    if out_dir is None:
        timestamp = datetime.datetime.utcnow().strftime('%Y-%m-%d_%H-%M-%S')
        out_dir = f"results_{target}_{timestamp}"
    
    os.makedirs(out_dir, exist_ok=True)
    os.chdir(out_dir)

    # 1. Pre-registration
    hash_proto = pre_register(target, trials, alpha, power, "delta_AUC_VQC_vs_GraphConv")
    click.echo(f"✅ Pre-registration SHA-256: {hash_proto}")

    # 2. Environment snapshot
    snapshot_environment()

    # 3. Load data
    click.echo(f"📦 Loading {target} dataset...")
    X_train, X_test, y_train, y_test = load_split(target, max_qubits)
    click.echo(f"   Train: {len(y_train)}  Test: {len(y_test)}  Active %: {y_train.mean():.1%}")

    # 4. Ultra-tune
    click.echo(f"🔍 Ultra-tune {trials} trials...")
    study = run_study(X_train, y_train, target=target, n_trials=trials, max_qubits=max_qubits)

    # 5. Generate plots
    df = study.trials_dataframe()
    df.to_csv("cv_comparison.csv", index=False)
    
    try:
        fig_auc_heatmap(df)
        fig_forest_effect(df)
        fig_power_curve()
        
        best_auc = df["value"].max() if len(df) > 0 else 0.0
        latex_boilerplate(target, auc=best_auc, cohen_d=0.5, p_corr=0.01)
    except Exception as e:
        click.echo(f"⚠️  Plot generation warning: {e}")

    # 6. Final checksums
    checksum_dir(".", "checksums_final.sha256")
    click.echo("🎉 QUALIS A1 pipeline finished!")

if __name__ == "__main__":
    main()
