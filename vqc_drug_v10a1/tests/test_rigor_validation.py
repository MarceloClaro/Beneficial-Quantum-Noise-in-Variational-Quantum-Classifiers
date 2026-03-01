import pandas as pd

from src.rigor import build_validation_report, write_validation_artifacts


def test_build_validation_report_structure(tmp_path):
    df = pd.DataFrame({"value": [0.71, 0.75, 0.73, 0.74, 0.72]})
    report = build_validation_report(df, alpha=0.05, n_permutations=500, seed=7)

    assert report["summary"]["n_trials"] == 5
    assert report["summary"]["auc_mean"] > 0.5
    assert "quality_gates" in report

    write_validation_artifacts(report, tmp_path)
    assert (tmp_path / "02_validation_report.json").exists()
    assert (tmp_path / "02_validation_report.md").exists()
