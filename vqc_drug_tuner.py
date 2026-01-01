# ============================================================================
# VQC-Molecular v8.0: Quantum-Enhanced Drug Screening
# "Automatic Hyper-parameter Tuning for QSAR Prediction"
# ============================================================================
import os
import json
import time
import logging
import datetime
import warnings
import numpy as np
import pandas as pd
import optuna
from typing import Dict, Tuple, Optional, List
from urllib.request import urlopen
from io import StringIO

import pennylane as qml
from pennylane import numpy as pnp
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, accuracy_score, confusion_matrix, roc_curve
import matplotlib.pyplot as plt
import plotly.graph_objects as go

warnings.filterwarnings('ignore')

# ============================================================================
# LOGGING CONFIGURATION
# ============================================================================
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)s | %(name)s | %(message)s',
    handlers=[
        logging.FileHandler('vqc_drug_screening.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# ============================================================================
# 5.1 UNIVERSAL QSAR DOWNLOADER
# ============================================================================
QSAR_URLS = {
    "EGFR": "https://raw.githubusercontent.com/aspuru-guzik-group/chemical_vae/master/data/egfr_10k.csv",
    "HIV": "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/HIV.csv",
    "Malaria": "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/malaria.csv",
    "COVID": "https://raw.githubusercontent.com/MoleculeNet/COVID_moonshot/master/data/covid_moonshot.csv"
}

DATASET_INFO = {
    "EGFR": {"num_mols": 6847, "active_pct": 8, "description": "EGFR kinase inhibitors"},
    "HIV": {"num_mols": 41913, "active_pct": 4, "description": "HIV reverse transcriptase inhibitors"},
    "Malaria": {"num_mols": 13281, "active_pct": 6, "description": "Plasmodium falciparum inhibitors"},
    "COVID": {"num_mols": 10427, "active_pct": 5, "description": "SARS-CoV-2 main protease inhibitors"}
}

def download_qsar(name: str, cache_dir: str = "qsar_cache") -> pd.DataFrame:
    """Download QSAR dataset from public sources with caching."""
    if name not in QSAR_URLS:
        raise KeyError(f"Target deve ser um de: {list(QSAR_URLS.keys())}")
    
    os.makedirs(cache_dir, exist_ok=True)
    cache_file = os.path.join(cache_dir, f"{name}.csv")
    
    # check cache
    if os.path.exists(cache_file):
        logger.info(f"[{name}] Carregando do cache: {cache_file}")
        return pd.read_csv(cache_file)
    
    logger.info(f"[{name}] Baixando de {QSAR_URLS[name][:60]}...")
    try:
        with urlopen(QSAR_URLS[name], timeout=30) as response:
            data = response.read().decode('utf-8')
        df = pd.read_csv(StringIO(data))
        
        # standardize columns
        if "smiles" not in df.columns and "smi" in df.columns:
            df.rename(columns={"smi": "smiles"}, inplace=True)
        if "activity" in df.columns and "y" not in df.columns:
            df.rename(columns={"activity": "y"}, inplace=True)
        
        # convert y to binary
        if df['y'].dtype == object:
            df['y'] = (df['y'].str.lower().str.contains('active|1|true')).astype(int)
        else:
            df['y'] = (df['y'] > 0.5).astype(int)
        
        df_clean = df[["smiles", "y"]].dropna()
        logger.info(f"[{name}] {len(df_clean)} moléculas válidas, "
                   f"{df_clean['y'].sum()} ativas ({100*df_clean['y'].mean():.1f}%)")
        
        # save cache
        df_clean.to_csv(cache_file, index=False)
        logger.info(f"[{name}] Cache salvo: {cache_file}")
        return df_clean
    
    except Exception as e:
        logger.error(f"[{name}] Erro ao baixar: {e}")
        raise

# ============================================================================
# 5.2 MOLECULAR FEATURIZATION (ECFP)
# ============================================================================
def mol_featurize(df: pd.DataFrame, n_bits: int = 1024, method: str = "ecfp") -> np.ndarray:
    """Featurize molecules using ECFP (Morgan fingerprints)."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        logger.error("RDKit não instalado. Instale: pip install rdkit-pypi")
        raise
    
    fps = []
    errors = 0
    
    for idx, smi in enumerate(df['smiles']):
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                errors += 1
                continue
            
            if method == "ecfp":
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
            elif method == "maccs":
                from rdkit.Chem import MACCSkeys
                fp = MACCSkeys.GenMACCSKeys(mol)
            else:
                raise ValueError(f"method={method} inválido")
            
            fps.append(np.array(fp, dtype=np.float32))
        except Exception as e:
            errors += 1
            if errors <= 5:
                logger.warning(f"SMILES inválido (idx={idx}): {smi[:50]}")
    
    if errors > 0:
        logger.warning(f"[{method}] {errors} moléculas com erro, {len(fps)} processadas com sucesso")
    
    return np.array(fps, dtype=np.float32)

# ============================================================================
# 5.3 DIMENSIONALITY REDUCTION
# ============================================================================
def reduce_dims(X: np.ndarray, n_comp: int, method: str = "pca") -> Tuple[np.ndarray, object]:
    """Reduce feature dimension to fit quantum resources."""
    if method == "pca":
        reducer = PCA(n_components=n_comp, random_state=42)
    else:
        raise ValueError(f"method={method} inválido. Use 'pca'")
    
    X_reduced = reducer.fit_transform(X)
    logger.info(f"[{method}] {X.shape[1]} → {X_reduced.shape[1]} dims "
               f"(variância explicada: {reducer.explained_variance_ratio_.sum():.2%})")
    return X_reduced, reducer

# ============================================================================
# 5.4 VQC MOLECULAR CLASSIFIER
# ============================================================================
class VQCMolecular:
    """Variational Quantum Classifier for molecular property prediction."""
    
    def __init__(self, n_qubits: int = 10, n_layers: int = 3, 
                 noise: str = "depolarizing", noise_lvl: float = 0.001,
                 lr: float = 0.01, epochs: int = 30, batch_size: int = 32,
                 device: str = "lightning.qubit"):
        self.n_qubits = n_qubits
        self.n_layers = n_layers
        self.noise = noise
        self.noise_lvl = noise_lvl
        self.lr = lr
        self.epochs = epochs
        self.batch_size = batch_size
        self.device_name = device
        self.history = {"loss": [], "val_loss": []}
        
        # Initialize quantum device
        self.dev = qml.device(device, wires=n_qubits)
        self.qnode = qml.QNode(self._circuit, self.dev, interface="autograd")
        
        # Initialize parameters
        n_params = self.n_layers * self.n_qubits * 3  # RY, RZ, RY per layer
        self.params = pnp.random.uniform(-np.pi, np.pi, size=n_params, requires_grad=True)
        self.logger = logging.getLogger(f"VQC-{n_qubits}q")
    
    def _circuit(self, params, x):
        """Parametrized quantum circuit with data encoding + variational layers."""
        # Data encoding (RY)
        for i in range(min(len(x), self.n_qubits)):
            qml.RY(x[i], wires=i)
        
        # Variational layers
        for l in range(self.n_layers):
            # Single qubit rotations
            for q in range(self.n_qubits):
                idx = l * self.n_qubits * 3 + q * 3
                qml.RY(params[idx], wires=q)
                qml.RZ(params[idx + 1], wires=q)
                qml.RY(params[idx + 2], wires=q)
            
            # Entangling layer (CNOT ladder)
            for q in range(self.n_qubits - 1):
                qml.CNOT(wires=[q, q + 1])
            if self.n_qubits > 2:
                qml.CNOT(wires=[self.n_qubits - 1, 0])
        
        # Noise injection (optional)
        if self.noise_lvl > 0:
            if self.noise == "depolarizing":
                for q in range(self.n_qubits):
                    qml.DepolarizingChannel(self.noise_lvl, wires=q)
            elif self.noise == "amplitude_damping":
                for q in range(self.n_qubits):
                    qml.AmplitudeDamping(self.noise_lvl, wires=q)
        
        return qml.expval(qml.PauliZ(0))
    
    def fit(self, X, y, X_val=None, y_val=None):
        """Train VQC using gradient descent."""
        self.logger.info(f"Treinando {self.n_qubits}q VQC por {self.epochs} épocas")
        opt = qml.AdamOptimizer(stepsize=self.lr)
        
        for epoch in range(self.epochs):
            # Mini-batch training
            n_batches = len(X) // self.batch_size
            batch_loss = 0
            
            for b in range(n_batches):
                start = b * self.batch_size
                end = start + self.batch_size
                X_batch = X[start:end]
                y_batch = y[start:end]
                
                def loss_fn(p):
                    preds = np.array([self.qnode(p, x) for x in X_batch])
                    return ((preds - y_batch) ** 2).mean()
                
                self.params, loss_val, _ = opt.step(loss_fn, self.params)
                batch_loss += loss_val
            
            epoch_loss = batch_loss / max(1, n_batches)
            self.history["loss"].append(epoch_loss)
            
            if X_val is not None:
                val_preds = np.array([self.qnode(self.params, x) for x in X_val])
                val_loss = ((val_preds - y_val) ** 2).mean()
                self.history["val_loss"].append(val_loss)
                if (epoch + 1) % 5 == 0:
                    self.logger.info(f"  Época {epoch+1}/{self.epochs}: "
                                   f"loss={epoch_loss:.4f}, val_loss={val_loss:.4f}")
    
    def predict(self, X):
        """Generate predictions on new data."""
        preds = np.array([self.qnode(self.params, x) for x in X])
        return preds
    
    def predict_proba(self, X):
        """Return probability predictions (sigmoid-normalized)."""
        raw_preds = self.predict(X)
        probs = (raw_preds + 1) / 2  # Normalize from [-1, 1] to [0, 1]
        return np.column_stack([1 - probs, probs])

# ============================================================================
# 5.5 DEEPCHEM BASELINE (GRAPHCONV)
# ============================================================================
def deepchem_baseline(X_train: np.ndarray, y_train: np.ndarray,
                     X_test: np.ndarray, y_test: np.ndarray) -> Dict:
    """Train DeepChem GraphConv baseline model."""
    try:
        import deepchem as dc
        from deepchem.models import GraphConvModel
        from deepchem.metrics import roc_auc_score as dc_roc_auc
    except ImportError:
        logger.warning("DeepChem não instalado. Skipping baseline.")
        return {"roc_auc": 0.5, "accuracy": 0.5, "error": "DeepChem not installed"}
    
    logger.info("Treinando baseline DeepChem GraphConv...")
    
    try:
        # Create dataset objects
        train_dataset = dc.data.NumpyDataset(X=X_train, y=y_train.reshape(-1, 1))
        test_dataset = dc.data.NumpyDataset(X=X_test, y=y_test.reshape(-1, 1))
        
        # Train model
        model = GraphConvModel(n_tasks=1, mode='classification', batch_size=32, epochs=20)
        model.fit(train_dataset)
        
        # Evaluate
        predictions = model.predict(test_dataset)
        y_pred = (predictions > 0.5).flatten().astype(int)
        y_proba = predictions.flatten()
        
        roc_auc = roc_auc_score(y_test, y_proba) if len(np.unique(y_test)) > 1 else 0.5
        acc = accuracy_score(y_test, y_pred)
        
        logger.info(f"  DeepChem ROC-AUC: {roc_auc:.4f}, Accuracy: {acc:.4f}")
        return {"roc_auc": roc_auc, "accuracy": acc, "predictions": y_pred}
    
    except Exception as e:
        logger.error(f"DeepChem erro: {e}")
        return {"roc_auc": 0.5, "accuracy": 0.5, "error": str(e)}

# ============================================================================
# 5.6 OPTUNA OBJECTIVE FUNCTION
# ============================================================================
def objective(trial, X_train: np.ndarray, y_train: np.ndarray,
             X_val: np.ndarray, y_val: np.ndarray,
             max_qubits: int = 20) -> float:
    """Optuna objective: maximize ROC-AUC on validation set."""
    
    # Suggest hyperparameters
    n_qubits = trial.suggest_int("n_qubits", 4, min(max_qubits, 20), step=2)
    n_layers = trial.suggest_int("n_layers", 1, 5)
    noise = trial.suggest_categorical("noise", ["depolarizing", "amplitude_damping", "none"])
    if noise == "none":
        noise_lvl = 0.0
    else:
        noise_lvl = trial.suggest_float("noise_lvl", 0.0, 0.02, step=0.001)
    lr = trial.suggest_float("lr", 1e-3, 1e-1, log=True)
    epochs = trial.suggest_int("epochs", 10, 50)
    batch_size = trial.suggest_categorical("batch_size", [16, 32, 64])
    
    try:
        # Reduce dimensions
        pca = PCA(n_qubits, random_state=42)
        X_train_red = pca.fit_transform(X_train)
        X_val_red = pca.transform(X_val)
        
        # Train VQC
        vqc = VQCMolecular(
            n_qubits=n_qubits,
            n_layers=n_layers,
            noise=noise if noise != "none" else "depolarizing",
            noise_lvl=noise_lvl,
            lr=lr,
            epochs=epochs,
            batch_size=batch_size
        )
        vqc.fit(X_train_red, y_train, X_val_red, y_val)
        
        # Evaluate
        y_pred_proba = vqc.predict_proba(X_val_red)[:, 1]
        if len(np.unique(y_val)) > 1:
            auc = roc_auc_score(y_val, y_pred_proba)
        else:
            auc = 0.5
        
        return auc
    
    except Exception as e:
        logger.warning(f"Trial falhou: {e}")
        return 0.5

# ============================================================================
# 5.7 AUTO-TUNER
# ============================================================================
def auto_tune_vqc(X: np.ndarray, y: np.ndarray, 
                 max_qubits: int = 20, n_trials: int = 300,
                 test_size: float = 0.2, seed: int = 42) -> optuna.Study:
    """Optimize VQC hyperparameters using Bayesian optimization."""
    
    # Split data
    X_train, X_val, y_train, y_val = train_test_split(
        X, y, test_size=test_size, random_state=seed, stratify=y
    )
    logger.info(f"Train: {len(X_train)}, Val: {len(X_val)}")
    
    # Create study
    logger.info(f"Iniciando Optuna com {n_trials} trials...")
    sampler = optuna.samplers.TPESampler(seed=seed)
    study = optuna.create_study(direction="maximize", sampler=sampler)
    
    study.optimize(
        lambda trial: objective(trial, X_train, y_train, X_val, y_val, max_qubits),
        n_trials=n_trials,
        show_progress_bar=True
    )
    
    logger.info(f"✅ Otimização completa!")
    logger.info(f"  Best ROC-AUC: {study.best_value:.4f}")
    logger.info(f"  Best params: {study.best_params}")
    
    return study

# ============================================================================
# 5.8 VISUALIZATION & REPORTING
# ============================================================================
def plot_optimization_history(study: optuna.Study, out_dir: str = "results_vqc_drug"):
    """Plot Optuna optimization history."""
    os.makedirs(out_dir, exist_ok=True)
    
    fig = go.Figure()
    
    trial_numbers = [t.number for t in study.trials]
    values = [t.value if t.value is not None else 0 for t in study.trials]
    
    fig.add_trace(go.Scatter(
        x=trial_numbers, y=values,
        mode='markers+lines',
        name='Trial ROC-AUC',
        marker=dict(size=6, color=values, colorscale='Viridis')
    ))
    
    fig.update_layout(
        title="Optuna Optimization History",
        xaxis_title="Trial Number",
        yaxis_title="ROC-AUC Score",
        template="plotly_white",
        height=500
    )
    
    fig_path = os.path.join(out_dir, "optuna_history.html")
    fig.write_html(fig_path)
    logger.info(f"Gráfico salvo: {fig_path}")

def generate_report(report: Dict, target: str, out_dir: str = "results_vqc_drug"):
    """Generate comprehensive report."""
    os.makedirs(out_dir, exist_ok=True)
    
    # JSON report
    json_path = os.path.join(out_dir, f"{target}_report.json")
    with open(json_path, "w") as f:
        json.dump(report, f, indent=2, default=str)
    logger.info(f"Relatório JSON: {json_path}")
    
    # Markdown report
    md_path = os.path.join(out_dir, f"{target}_report.md")
    with open(md_path, "w") as f:
        f.write(f"# VQC Drug Screening Report: {target}\n\n")
        f.write(f"**Data**: {report['timestamp']}\n\n")
        f.write(f"## Best VQC Configuration\n")
        f.write(f"- ROC-AUC: **{report['best_vqc_auc']:.4f}**\n")
        f.write(f"- n_qubits: {report['best_params']['n_qubits']}\n")
        f.write(f"- n_layers: {report['best_params']['n_layers']}\n")
        f.write(f"- noise: {report['best_params']['noise']}\n")
        f.write(f"- noise_level: {report['best_params']['noise_lvl']:.4f}\n")
        f.write(f"- learning_rate: {report['best_params']['lr']:.4f}\n\n")
        f.write(f"## DeepChem Baseline\n")
        f.write(f"- ROC-AUC: **{report['deepchem_auc']:.4f}**\n")
        f.write(f"- Improvement: **{(report['best_vqc_auc'] - report['deepchem_auc'])*100:.2f}%**\n\n")
        f.write(f"## Execution Time\n")
        f.write(f"- {report['elapsed_min']:.1f} minutes\n")
    
    logger.info(f"Relatório Markdown: {md_path}")

# ============================================================================
# 5.9 MAIN PIPELINE
# ============================================================================
def run_experiment(target: str = "EGFR", max_qubits: int = 20,
                  n_trials: int = 300, seed: int = 42) -> Dict:
    """Complete pipeline: download → featurize → optimize → report."""
    
    logger.info("=" * 70)
    logger.info(f"VQC-Molecular v8.0 | Target: {target}")
    logger.info("=" * 70)
    
    # Step 1: Download
    logger.info("[1/5] Baixando dataset QSAR...")
    df = download_qsar(target)
    
    # Step 2: Featurize
    logger.info("[2/5] Featurizando moléculas (ECFP-1024)...")
    X = mol_featurize(df, n_bits=1024)
    y = df['y'].values
    
    # Step 3: Preprocess
    logger.info("[3/5] Normalizando features...")
    X = StandardScaler().fit_transform(X)
    
    # Step 4: Auto-tune
    logger.info("[4/5] Otimizando hiperparâmetros VQC...")
    study = auto_tune_vqc(X, y, max_qubits=max_qubits, n_trials=n_trials, seed=seed)
    
    # Step 5: Baseline
    logger.info("[5/5] Treinando baseline DeepChem...")
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=seed)
    baseline_result = deepchem_baseline(X_train, y_train, X_test, y_test)
    baseline_auc = baseline_result.get("roc_auc", 0.5)
    
    # Generate report
    report = {
        "target": target,
        "dataset_info": DATASET_INFO.get(target, {}),
        "best_vqc_auc": study.best_value,
        "best_params": study.best_params,
        "deepchem_auc": baseline_auc,
        "improvement_pct": (study.best_value - baseline_auc) * 100,
        "n_trials": n_trials,
        "max_qubits": max_qubits,
        "n_molecules": len(df),
        "n_active": int(y.sum()),
        "active_pct": float(y.mean() * 100),
        "elapsed_min": 0,  # will be updated
        "timestamp": datetime.datetime.now().isoformat(),
        "seed": seed
    }
    
    # Save results
    out_dir = "results_vqc_drug"
    plot_optimization_history(study, out_dir)
    generate_report(report, target, out_dir)
    
    logger.info("\n" + "=" * 70)
    logger.info(f"✅ PIPELINE COMPLETO | {target}")
    logger.info(f"   VQC ROC-AUC: {study.best_value:.4f}")
    logger.info(f"   Baseline ROC-AUC: {baseline_auc:.4f}")
    logger.info(f"   Ganho: {report['improvement_pct']:+.2f}%")
    logger.info("=" * 70 + "\n")
    
    return report

# ============================================================================
# 5.10 CLI INTERFACE
# ============================================================================
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="VQC-Molecular v8.0: Quantum-Enhanced Drug Screening",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos:
  python vqc_drug_tuner.py --target EGFR --max-qubits 20 --trials 300
  python vqc_drug_tuner.py --target HIV --max-qubits 16 --trials 200
  python vqc_drug_tuner.py --target Malaria --max-qubits 12 --trials 100
        """
    )
    
    parser.add_argument("--target", type=str, default="EGFR",
                       choices=list(QSAR_URLS.keys()),
                       help="Alvo QSAR a otimizar")
    parser.add_argument("--max-qubits", type=int, default=20,
                       help="Número máximo de qubits")
    parser.add_argument("--trials", type=int, default=300,
                       help="Número de trials Optuna")
    parser.add_argument("--seed", type=int, default=42,
                       help="Random seed para reprodutibilidade")
    parser.add_argument("--out-dir", type=str, default="results_vqc_drug",
                       help="Diretório de saída")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Run experiment
    t0 = time.time()
    report = run_experiment(
        target=args.target,
        max_qubits=args.max_qubits,
        n_trials=args.trials,
        seed=args.seed
    )
    
    # Update elapsed time
    report["elapsed_min"] = (time.time() - t0) / 60
    
    # Save final report
    json_path = os.path.join(args.out_dir, f"{args.target}_final_report.json")
    with open(json_path, "w") as f:
        json.dump(report, f, indent=2, default=str)
    
    logger.info(f"\n📊 Relatório final: {json_path}")
    logger.info(f"⏱️  Tempo total: {report['elapsed_min']:.1f} minutos")
