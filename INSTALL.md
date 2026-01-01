# Guia de Instalação - Framework v7.2

## 📋 Pré-requisitos

### Sistema Operacional
- Windows 10/11, Linux (Ubuntu 20.04+), ou macOS 10.15+
- Python 3.9 ou superior


### Hardware Recomendado
- **Mínimo**: 8GB RAM, 4 cores
- **Recomendado**: 16GB RAM, 8+ cores
- **Armazenamento**: 10-20GB livres (para resultados completos)


---


## 🚀 Instalação Rápida

### 1. Clone o Repositório

```bash
git clone <https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers.git>
cd Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers

```text

### 2. Crie um Ambiente Virtual (Recomendado)

#### Windows PowerShell

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1

```text

#### Linux/Mac

```bash
python3 -m venv .venv
source .venv/bin/activate

```text

### 3. Instale as Dependências

```bash
pip install --upgrade pip
pip install -r requirements.txt

```text

### 4. Verifique a Instalação

```bash
python -c "import pennylane; print(f'PennyLane: {pennylane.__version__}')"
python -c "import optuna; print(f'Optuna: {optuna.__version__}')"

```text

**Saída esperada**:

```

PennyLane: 0.38.0
Optuna: 3.x.x

```text

---


## 📦 Dependências Principais

| Pacote | Versão | Propósito |
|--------|--------|-----------|
| **pennylane** | ≥0.38.0 | Framework de computação quântica |
| **numpy** | ≥1.24.0 | Computação numérica |
| **pandas** | ≥2.0.0 | Manipulação de dados |
| **scikit-learn** | ≥1.3.0 | Algoritmos clássicos de ML |
| **scipy** | ≥1.11.0 | Análises estatísticas |
| **statsmodels** | ≥0.14.0 | ANOVA e testes estatísticos |
| **plotly** | ≥5.17.0 | Visualizações interativas |
| **optuna** | ≥3.3.0 | Otimização Bayesiana (opcional) |

---


## ✅ Teste de Funcionamento

### Teste Rápido (2-3 minutos)

```bash

# Windows
$env:VQC_QUICK="1"; $env:VQC_BAYESIAN="1"; python framework_investigativo_completo.py

# Linux/Mac
VQC_QUICK=1 VQC_BAYESIAN=1 python framework_investigativo_completo.py

```text

**O que esperar**:
- Criação de pasta `resultados_YYYY-MM-DD_HH-MM-SS/`
- Mensagens de progresso com emojis
- Geração de CSVs e visualizações
- Consolidação automática ao final


---


## 🐛 Troubleshooting

### Erro: "ModuleNotFoundError: No module named 'pennylane'"
**Solução**:

```bash
pip install pennylane

```text

### Erro: "No module named 'optuna'" (modo Bayesiano)
**Solução** (opcional, apenas para modo Bayesiano):

```bash
pip install optuna

```text

### Erro: "ImportError: DLL load failed" (Windows)
**Solução**: Instale Microsoft Visual C++ Redistributable
- Download: <https://aka.ms/vs/17/release/vc_redist.x64.exe>


### Aviso: "UserWarning: PennyLane device..."
**Não é um erro**: Avisos podem ser ignorados, não afetam execução.


---


## 🔧 Instalação Avançada

### Com CUDA (GPU - Opcional)
Se você tem GPU NVIDIA e quer acelerar simulações:

```bash
pip install pennylane-lightning-gpu

```text

### Com Qiskit (Opcional)
Para integração com IBM Quantum:

```bash
pip install pennylane-qiskit

```text

### Ambiente Conda (Alternativa)

```bash
conda create -n vqc python=3.9
conda activate vqc
pip install -r requirements.txt

```text

---


## 📂 Estrutura do Projeto

Após instalação, você terá:

```

Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers/
├── .gitignore                              # Arquivos ignorados pelo Git
├── .venv/                                  # Ambiente virtual (ignorado)
├── LICENSE                                 # Licença MIT
├── README.md                               # Documentação principal
├── INSTALL.md                              # Este arquivo
├── requirements.txt                        # Dependências Python
├── framework_investigativo_completo.py     # Framework principal
├── docs/                                   # Documentação detalhada
│   ├── AUTOMACAO_FRAMEWORK.md
│   ├── CHANGELOG_v7.2.md
│   ├── GUIA_RAPIDO_v7.2.md
│   └── RESUMO_EXECUTIVO_v7.2.md
├── examples/                               # Exemplos de uso
│   └── exemplo_uso_programatico.py
└── tools/                                  # Scripts auxiliares (obsoletos)
    ├── consolidate_results.py              # (funcionalidade integrada no framework)
    └── orchestrate_framework.py            # (funcionalidade integrada no framework)

```

---


## 🎯 Próximos Passos

1. ✅ **Leia a documentação**: [README.md](README.md)
2. ✅ **Execute o teste rápido** (comando acima)
3. ✅ **Explore os modos de execução**: [GUIA_RAPIDO_v7.2.md](docs/GUIA_RAPIDO_v7.2.md)
4. ✅ **Execute o framework completo** quando estiver pronto


---


## 📞 Suporte

- **Issues**: [GitHub Issues](https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers/issues)
- **Documentação**: Consulte `docs/` para guias detalhados
- **Email**: (adicione seu email de contato)


---


## 📝 Notas de Versão

### v7.2 (Atual)
- Consolidação automática integrada
- Comparação de baselines automática
- Metadados automáticos
- Documentação completa


### Requisitos Mínimos
- Python 3.9+
- 8GB RAM
- 5GB disco (mínimo)


---


**✅ Instalação concluída! Execute `python framework_investigativo_completo.py` para começar.**

