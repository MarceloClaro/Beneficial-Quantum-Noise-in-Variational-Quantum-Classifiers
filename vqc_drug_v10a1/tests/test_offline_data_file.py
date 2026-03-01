import pandas as pd

from src.data import load_split


def test_load_split_from_data_file(tmp_path):
    df = pd.DataFrame({
        "smiles": ["CCO", "CCN", "CCC", "CCCl", "c1ccccc1", "CC(=O)O", "CCS", "COC", "CCBr", "CCF"],
        "activity": [1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
    })
    path = tmp_path / "custom.csv"
    df.to_csv(path, index=False)

    X_train, X_test, y_train, y_test = load_split("EGFR", n_qubits=4, seed=42, data_file=str(path))
    assert X_train.shape[1] == 4
    assert len(y_train) > 0 and len(y_test) > 0
