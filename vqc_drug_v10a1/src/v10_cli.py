"""
CLI estendida v10.1-A1 (k8-cli)
Orquestra módulos Next-Gen: Fisher, Lindblad, Meta, QASR, WY, DD.
"""
from __future__ import annotations
import os
import datetime
import click
from .audit import pre_register, checksum_dir, snapshot_environment
from .data import load_split
from .tune import run_study
from .v10_stats import qualis_report_from_trials


@click.command()
@click.option("--target", default="HIV", help="EGFR | HIV | Malaria | COVID")
@click.option("--trials", default=200, help="Número de trials Optuna")
@click.option("--max-qubits", default=20, help="Limite de qubits na featurização")
@click.option("--out-dir", default=None, help="Diretório de saída")
@click.option("--alpha", default=0.05, help="FWER")
@click.option("--power", default=0.8, help="Power mínimo")
def main(target: str, trials: int, max_qubits: int, out_dir: str | None, alpha: float, power: float):
    if out_dir is None:
        out_dir = f"results_{target}_{datetime.datetime.utcnow().strftime('%Y%m%d_%H%M%S')}"
    os.makedirs(out_dir, exist_ok=True)
    os.chdir(out_dir)

    click.echo(f"🚀 VQC-Molecular v10.1-A1 – target={target}, trials={trials}")

    # 1. Pre-register
    proto_hash = pre_register(target, trials, alpha, power, "delta_AUC_v10.1")
    click.echo(f"✅ Pre-registration SHA-256: {proto_hash}")

    # 2. Snapshot env
    snapshot_environment()

    # 3. Load data
    click.echo("📦 Carregando dataset...")
    X_train, X_test, y_train, y_test = load_split(target, max_qubits)
    click.echo(f"   Train={len(y_train)} | Test={len(y_test)} | Ativos={y_train.mean():.2%}")

    # 4. Study
    click.echo("🔍 Otimização (QASR + Fisher + Meta + Lindblad)...")
    study = run_study(X_train, y_train, target=target, n_trials=trials, max_qubits=max_qubits)

    # 5. Qualis report
    df = study.trials_dataframe()
    df.to_csv("trials_v10_1.csv", index=False)
    qreport = qualis_report_from_trials(df)
    qreport.to_csv("qualis_a1_stats.csv", index=False)

    best_auc = df["value"].max() if not df.empty else 0.0
    click.echo(f"🎯 Best AUC={best_auc:.4f}")

    # 6. Checksums
    checksum_dir(".", "checksums_v10_1.sha256")
    click.echo("🎉 Pipeline v10.1 finalizado")


if __name__ == "__main__":
    main()

