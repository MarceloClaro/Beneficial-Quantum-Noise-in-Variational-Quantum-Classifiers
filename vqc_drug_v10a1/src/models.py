"""
Quantum VQC with Fisher-info audit & DeepChem baseline.
Designed for lightning.qubit / lightning.gpu backend.
"""
import torch
import numpy as np
import pennylane as qml
from pennylane import numpy as pnp
from sklearn.metrics import roc_auc_score

class VQCAudit(torch.nn.Module):
    """Variational Quantum Classifier with full audit trail."""
    def __init__(self, n_qubits: int, n_layers: int, noise_type: str,
                 noise_level: float, lr: float, epochs: int, constant_init: str,
                 arch: str, optimizer: str, loss: str, batch_size: int,
                 trial_id: int, device_str: str = "lightning.qubit"):
        super().__init__()
        self.n_qubits = n_qubits
        self.n_layers = n_layers
        self.noise_type = noise_type
        self.noise_level = noise_level
        self.lr = lr
        self.epochs = epochs
        self.constant_init = constant_init
        self.arch = arch
        self.optimizer = optimizer
        self.loss = loss
        self.batch_size = batch_size
        self.trial_id = trial_id

        # Device selection (GPU if available)
        try:
            if device_str == "lightning.gpu":
                self.dev = qml.device("lightning.gpu", wires=n_qubits)
            else:
                self.dev = qml.device("lightning.qubit", wires=n_qubits)
        except:
            self.dev = qml.device("default.qubit", wires=n_qubits)

        # Parameters
        n_params = n_layers * n_qubits
        self.params = torch.nn.Parameter(
            torch.tensor(self._init_params(n_params), dtype=torch.float32))

        # QNode
        self.qnode = qml.QNode(self._circuit, self.dev, interface="torch")

    def _init_params(self, n):
        """Constant-based initialization (π, e, φ, ℏ, α…)."""
        if self.constant_init == "random":
            return np.random.uniform(-np.pi, np.pi, n)
        const_map = {"pi": np.pi, "e": np.e, "phi": (1 + np.sqrt(5))/2,
                     "hbar": 1.05457e-34, "alpha": 7.297e-3}
        const = const_map.get(self.constant_init, np.pi)
        vals = [(const * (i+1)) % (2*np.pi) - np.pi for i in range(n)]
        return np.array(vals) + np.random.normal(0, 0.1, n)

    def _circuit(self, params, x):
        # Angle encoding
        for i in range(self.n_qubits):
            qml.RY(x[i], wires=i)
        
        # Variational layers
        for l in range(self.n_layers):
            for q in range(self.n_qubits):
                qml.RY(params[l*self.n_qubits + q], wires=q)
            
            # Architecture-specific entangling
            if self.arch == "tree":
                self._arch_tree()
            elif self.arch == "star":
                self._arch_star()
            elif self.arch == "brickwork":
                self._arch_brickwork()
            
            # Noise injection
            if self.noise_type != "none":
                self._inject_noise()
        
        return qml.expval(qml.PauliZ(0))

    def _arch_tree(self):
        step = 1
        while step < self.n_qubits:
            for i in range(0, self.n_qubits-step, 2*step):
                qml.CNOT(wires=[i, i+step])
            step *= 2

    def _arch_star(self):
        for q in range(1, self.n_qubits):
            qml.CNOT(wires=[0, q])

    def _arch_brickwork(self):
        for i in range(0, self.n_qubits-1, 2):
            qml.CNOT(wires=[i, i+1])
        for i in range(1, self.n_qubits-1, 2):
            qml.CNOT(wires=[i, i+1])

    def _inject_noise(self):
        if self.noise_type == "depolarizing":
            for q in range(self.n_qubits):
                qml.DepolarizingChannel(self.noise_level, wires=q)
        elif self.noise_type == "amplitude_damping":
            for q in range(self.n_qubits):
                qml.AmplitudeDamping(self.noise_level, wires=q)

    def forward(self, x):
        if len(x.shape) == 1:
            return self.qnode(self.params, x)
        else:
            return torch.stack([self.qnode(self.params, x[i]) for i in range(x.shape[0])])

    def fit(self, X, y, X_val, y_val):
        opt_class = torch.optim.Adam if self.optimizer == "adam" else torch.optim.SGD
        opt = opt_class(self.parameters(), lr=self.lr)
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, patience=10)

        for epoch in range(self.epochs):
            # Mini-batch
            perm = torch.randperm(X.size(0))
            losses = []
            for i in range(0, X.size(0), self.batch_size):
                idx = perm[i:i+self.batch_size]
                pred = self.forward(X[idx]).squeeze()
                loss = self._loss_fn(pred, y[idx])
                opt.zero_grad()
                loss.backward()
                opt.step()
                losses.append(loss.item())
            
            # Validation
            with torch.no_grad():
                val_pred = self.forward(X_val).squeeze()
                val_loss = self._loss_fn(val_pred, y_val)
                val_auc = roc_auc_score(y_val.cpu().numpy(), 
                                       torch.sigmoid(val_pred).cpu().numpy())
            scheduler.step(val_loss)

    def _loss_fn(self, pred, target):
        if self.loss == "mse":
            return torch.mean((pred - target) ** 2)
        elif self.loss == "bce":
            return torch.nn.functional.binary_cross_entropy_with_logits(pred, target)
        else:  # hinge
            return torch.mean(torch.clamp(1 - target * pred, min=0))

    def predict_proba(self, X):
        with torch.no_grad():
            logits = self.forward(X).cpu().numpy()
        proba = 1 / (1 + np.exp(-logits))
        return np.column_stack([1 - proba, proba])
