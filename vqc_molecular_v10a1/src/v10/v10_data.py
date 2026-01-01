"""
Utilidades de dados v10.1-A1 (Qualis A1)

Função principal: ``download_chembl_latest`` para extrair conjuntos ChEMBL 32
por família (ex.: kinases) com IC50 < limite, retornando CSV padronizado
(smiles, y, split) e checksum SHA-256.

Dependências opcionais: ``chembl_webresource_client`` e ``rdkit``.
"""
from __future__ import annotations

import hashlib
import logging
import tempfile
from pathlib import Path
from typing import Iterable, List, Tuple

import pandas as pd

logger = logging.getLogger(__name__)


def _load_chembl_clients():
    try:
        from chembl_webresource_client.new_client import new_client  # type: ignore
    except ImportError as exc:  # pragma: no cover - aviso amigável
        raise ImportError(
            "download_chembl_latest requer chembl_webresource_client. "
            "Instale via `pip install chembl_webresource_client`."
        ) from exc
    try:
        from rdkit import Chem  # type: ignore
    except ImportError as exc:  # pragma: no cover - aviso amigável
        raise ImportError(
            "download_chembl_latest requer RDKit. Sugestão: `conda install -c conda-forge rdkit`."
        ) from exc
    return new_client, Chem


def _compute_checksum(path: Path) -> str:
    sha = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(8192)
            if not chunk:
                break
            sha.update(chunk)
    return sha.hexdigest()


def download_chembl_latest(
    target_family: str = "kinase",
    max_mol: int = 15_000,
    ic50_threshold_nM: float = 100.0,
    train_frac: float = 0.8,
    random_state: int = 42,
    out_dir: str | Path | None = None,
) -> Path:
    """
    Baixa ChEMBL 32 filtrando por família e IC50, devolvendo CSV (smiles, y, split).

    Parameters
    ----------
    target_family : str
        Família-alvo (ex.: "kinase", "ppar").
    max_mol : int
        Número máximo de moléculas a coletar.
    ic50_threshold_nM : float
        Limite para label ativo (IC50 < limiar => y=1).
    train_frac : float
        Fração para split de treino (restante vira teste).
    random_state : int
        Semente para embaralhar.
    out_dir : str | Path | None
        Diretório de saída; se None, usa diretório temporário.

    Returns
    -------
    Path
        Caminho do CSV gerado.
    """
    new_client, Chem = _load_chembl_clients()

    targets = list(new_client.target.filter(family__icontains=target_family)[:50])
    target_ids: List[str] = [t.get("target_chembl_id") for t in targets if t.get("target_chembl_id")]
    if not target_ids:
        raise ValueError(f"Nenhum target encontrado para família '{target_family}'.")

    activities = new_client.activity.filter(
        target_chembl_id__in=target_ids,
        standard_type="IC50",
        standard_units="nM",
    )

    rows: List[Tuple[str, int]] = []
    for act in activities:
        if len(rows) >= max_mol:
            break
        mol_id = act.get("molecule_chembl_id")
        if not mol_id:
            continue
        mol_data = new_client.molecule.get(mol_id)
        smi = mol_data.get("molecule_structures", {}).get("canonical_smiles") if mol_data else None
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        ic50_val = act.get("standard_value") or act.get("ic50_value")
        try:
            ic50_num = float(ic50_val)
        except (TypeError, ValueError):
            continue
        label = int(ic50_num < ic50_threshold_nM)
        rows.append((smi, label))

    if not rows:
        raise RuntimeError("Nenhuma molécula coletada; verifique filtros ou conectividade com ChEMBL.")

    df = pd.DataFrame(rows, columns=["smiles", "y"]).dropna()
    df = df.sample(frac=1.0, random_state=random_state).reset_index(drop=True)

    n_train = max(1, int(len(df) * train_frac))
    df["split"] = "train"
    df.loc[df.index[n_train:], "split"] = "test"

    out_dir = Path(out_dir) if out_dir is not None else Path(tempfile.gettempdir())
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"chembl_{target_family}_latest.csv"
    df.to_csv(csv_path, index=False)

    sha256 = _compute_checksum(csv_path)
    logger.info("ChEMBL %s baixado: %d mols, SHA-256=%s", target_family, len(df), sha256)
    return csv_path


__all__ = ["download_chembl_latest"]
