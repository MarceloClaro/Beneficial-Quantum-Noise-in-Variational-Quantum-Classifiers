#!/usr/bin/env python3
"""
Script para gerar notebook 01_introducao_vqc.ipynb completo
Introdução didática aos VQCs mantendo rigor QUALIS A1
"""

import json
from pathlib import Path

def create_intro_notebook():
    """Cria notebook de introdução aos VQCs."""
    
    notebook = {
        "cells": [],
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3"
            },
            "language_info": {
                "codemirror_mode": {"name": "ipython", "version": 3},
                "file_extension": ".py",
                "mimetype": "text/x-python",
                "name": "python",
                "nbconvert_exporter": "python",
                "pygments_lexer": "ipython3",
                "version": "3.9.0"
            }
        },
        "nbformat": 4,
        "nbformat_minor": 4
    }
    
    # Célula 1: Título
    notebook["cells"].append({
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "# Introdução aos Variational Quantum Classifiers (VQCs)\\n",
            "\\n",
            "[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers/blob/main/notebooks/01_introducao_vqc.ipynb)\\n",
            "\\n",
            "---\\n",
            "\\n",
            "## 🎯 Objetivos deste Notebook\\n",
            "\\n",
            "Este notebook fornece uma introdução prática e didática aos Classificadores Quânticos Variacionais (VQCs),\\n",
            "preparando você para entender os conceitos avançados do framework completo.\\n",
            "\\n",
            "### O que você aprenderá:\\n",
            "\\n",
            "1. 🔬 **O que é um VQC** - conceito básico e motivação\\n",
            "2. 🧮 **Como funciona** - arquitetura e treinamento\\n",
            "3. 💻 **Implementação prática** - código passo a passo\\n",
            "4. 📊 **Primeiro experimento** - treinar um VQC em dados sintéticos\\n",
            "5. 🎨 **Visualização** - entender o processo de aprendizado\\n",
            "\\n",
            "---"
        ]
    })
    
    # Célula 2: O que é um VQC
    notebook["cells"].append({
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 1. O que é um Variational Quantum Classifier?\\n",
            "\\n",
            "### 💡 Analogia para Iniciantes\\n",
            "\\n",
            "Imagine que você tem um filtro especial de café que pode ser ajustado (girar manivelas, mudar filtros).\\n",
            "Cada ajuste muda como o café é filtrado. Um VQC é como esse filtro ajustável, mas para dados:\\n",
            "\\n",
            "- **Entrada**: dados (como números)\\n",
            "- **Filtro**: circuito quântico com parâmetros ajustáveis\\n",
            "- **Saída**: classificação (\"gato\" ou \"cachorro\", por exemplo)\\n",
            "- **Treinamento**: ajustar os parâmetros para melhorar a precisão\\n",
            "\\n",
            "### 🎓 Definição Técnica\\n",
            "\\n",
            "Um **Variational Quantum Classifier (VQC)** é um algoritmo híbrido quântico-clássico que:\\n",
            "\\n",
            "1. **Codifica** dados clássicos em estados quânticos (feature map)\\n",
            "2. **Processa** através de um circuito quântico parametrizado (ansatz variacional)\\n",
            "3. **Mede** observáveis quânticos para obter predições\\n",
            "4. **Otimiza** parâmetros classicamente via gradient descent\\n",
            "\\n",
            "**Equação fundamental**:\\n",
            "\\n",
            "$$f(x; \\\\theta) = \\\\langle \\\\psi(x, \\\\theta) | \\\\hat{O} | \\\\psi(x, \\\\theta) \\\\rangle$$\\n",
            "\\n",
            "Onde:\\n",
            "- $x$ = dados de entrada\\n",
            "- $\\\\theta$ = parâmetros variacionais\\n",
            "- $|\\\\psi\\\\rangle$ = estado quântico preparado\\n",
            "- $\\\\hat{O}$ = observável medido (ex: $\\\\sigma_z$)\\n",
            "\\n",
            "---"
        ]
    })
    
    # Célula 3: Instalação
    notebook["cells"].append({
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 2. Instalação e Setup\\n",
            "\\n",
            "Vamos instalar as bibliotecas necessárias:"
        ]
    })
    
    notebook["cells"].append({
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "%%capture\\n",
            "!pip install pennylane numpy matplotlib scikit-learn\\n",
            "\\n",
            "print('✓ Bibliotecas instaladas!')"
        ]
    })
    
    # Célula 4: Imports
    notebook["cells"].append({
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "import numpy as np\\n",
            "import matplotlib.pyplot as plt\\n",
            "from sklearn.datasets import make_moons\\n",
            "from sklearn.model_selection import train_test_split\\n",
            "from sklearn.preprocessing import StandardScaler\\n",
            "\\n",
            "import pennylane as qml\\n",
            "from pennylane import numpy as pnp\\n",
            "\\n",
            "print(f'PennyLane version: {qml.__version__}')\\n",
            "print('✓ Imports concluídos!')"
        ]
    })
    
    # Célula 5: Criar dados sintéticos
    notebook["cells"].append({
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 3. Preparar Dados Sintéticos\\n",
            "\\n",
            "Vamos criar um dataset sintético chamado \"two moons\" - dois crescentes entrelaçados.\\n",
            "Este é um problema de classificação não-linear clássico."
        ]
    })
    
    notebook["cells"].append({
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# Gerar dados\\n",
            "np.random.seed(42)\\n",
            "X, y = make_moons(n_samples=200, noise=0.1)\\n",
            "\\n",
            "# Normalizar\\n",
            "scaler = StandardScaler()\\n",
            "X = scaler.fit_transform(X)\\n",
            "\\n",
            "# Converter para +1/-1\\n",
            "y = 2 * y - 1\\n",
            "\\n",
            "# Dividir treino/teste\\n",
            "X_train, X_test, y_train, y_test = train_test_split(\\n",
            "    X, y, test_size=0.2, random_state=42\\n",
            ")\\n",
            "\\n",
            "# Visualizar\\n",
            "plt.figure(figsize=(8, 6))\\n",
            "plt.scatter(X_train[y_train==1, 0], X_train[y_train==1, 1], \\n",
            "            c='blue', label='Classe +1', alpha=0.6, edgecolors='k')\\n",
            "plt.scatter(X_train[y_train==-1, 0], X_train[y_train==-1, 1], \\n",
            "            c='red', label='Classe -1', alpha=0.6, edgecolors='k')\\n",
            "plt.xlabel('Feature 1')\\n",
            "plt.ylabel('Feature 2')\\n",
            "plt.title('Dataset Two Moons (Dados de Treinamento)')\\n",
            "plt.legend()\\n",
            "plt.grid(True, alpha=0.3)\\n",
            "plt.show()\\n",
            "\\n",
            "print(f'Treino: {len(X_train)} amostras')\\n",
            "print(f'Teste: {len(X_test)} amostras')"
        ]
    })
    
    # Célula 6: Definir circuito
    notebook["cells"].append({
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 4. Definir o Circuito Quântico\\n",
            "\\n",
            "### 💡 Para Iniciantes\\n",
            "\\n",
            "O circuito quântico é como uma \"receita\" de operações nos qubits:\\n",
            "1. **Codificar** os dados (feature map)\\n",
            "2. **Processar** com portas parametrizadas (ansatz)\\n",
            "3. **Medir** o resultado\\n",
            "\\n",
            "### 🎓 Para Especialistas\\n",
            "\\n",
            "Implementamos um ansatz simples com:\\n",
            "- Feature map: ângulos baseados nas features\\n",
            "- Ansatz variacional: RY(θ) + CNOT para emaranhamento\\n",
            "- Medição: Pauli-Z no primeiro qubit"
        ]
    })
    
    notebook["cells"].append({
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# Número de qubits\\n",
            "n_qubits = 2\\n",
            "\\n",
            "# Dispositivo quântico (simulador)\\n",
            "dev = qml.device('default.qubit', wires=n_qubits)\\n",
            "\\n",
            "@qml.qnode(dev)\\n",
            "def quantum_circuit(weights, x):\\n",
            "    \\\"\\\"\\\"\\n",
            "    Circuito quântico variacional simples.\\n",
            "    \\n",
            "    Args:\\n",
            "        weights: parâmetros treináveis\\n",
            "        x: vetor de features (2D)\\n",
            "    \\n",
            "    Returns:\\n",
            "        Expectation value de Pauli-Z\\n",
            "    \\\"\\\"\\\"\\n",
            "    # Feature map: codificar dados\\n",
            "    qml.RY(x[0], wires=0)\\n",
            "    qml.RY(x[1], wires=1)\\n",
            "    \\n",
            "    # Ansatz variacional (camada 1)\\n",
            "    qml.RY(weights[0], wires=0)\\n",
            "    qml.RY(weights[1], wires=1)\\n",
            "    qml.CNOT(wires=[0, 1])\\n",
            "    \\n",
            "    # Ansatz variacional (camada 2)\\n",
            "    qml.RY(weights[2], wires=0)\\n",
            "    qml.RY(weights[3], wires=1)\\n",
            "    qml.CNOT(wires=[1, 0])\\n",
            "    \\n",
            "    # Medição\\n",
            "    return qml.expval(qml.PauliZ(0))\\n",
            "\\n",
            "# Testar circuito\\n",
            "test_weights = pnp.array([0.1, 0.2, 0.3, 0.4], requires_grad=True)\\n",
            "test_x = pnp.array([0.5, -0.3])\\n",
            "output = quantum_circuit(test_weights, test_x)\\n",
            "print(f'Saída do circuito (teste): {output:.4f}')\\n",
            "\\n",
            "# Desenhar circuito\\n",
            "print('\\\\nEstrutura do circuito:')\\n",
            "print(qml.draw(quantum_circuit)(test_weights, test_x))"
        ]
    })
    
    # Célula 7: Função de custo
    notebook["cells"].append({
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 5. Definir Função de Custo\\n",
            "\\n",
            "A função de custo mede o quão errada é a predição do modelo.\\n",
            "Queremos minimizar esse erro durante o treinamento."
        ]
    })
    
    notebook["cells"].append({
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "def cost_function(weights, X, y):\\n",
            "    \\\"\\\"\\\"\\n",
            "    Calcula o erro quadrático médio (MSE).\\n",
            "    \\n",
            "    Args:\\n",
            "        weights: parâmetros do circuito\\n",
            "        X: dados de entrada\\n",
            "        y: rótulos verdadeiros (+1 ou -1)\\n",
            "    \\n",
            "    Returns:\\n",
            "        MSE loss\\n",
            "    \\\"\\\"\\\"\\n",
            "    predictions = [quantum_circuit(weights, x) for x in X]\\n",
            "    loss = np.mean((predictions - y)**2)\\n",
            "    return loss\\n",
            "\\n",
            "def accuracy(weights, X, y):\\n",
            "    \\\"\\\"\\\"Calcula acurácia de classificação.\\\"\\\"\\\"\\n",
            "    predictions = [np.sign(quantum_circuit(weights, x)) for x in X]\\n",
            "    return np.mean(predictions == y)\\n",
            "\\n",
            "print('✓ Funções de custo definidas!')"
        ]
    })
    
    # Célula 8: Treinar
    notebook["cells"].append({
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 6. Treinar o VQC\\n",
            "\\n",
            "### 💡 Para Iniciantes\\n",
            "\\n",
            "Treinamento = ajustar os parâmetros do circuito gradualmente para reduzir o erro.\\n",
            "É como afinar um instrumento musical: pequenos ajustes até soar perfeito!\\n",
            "\\n",
            "### 🎓 Para Especialistas\\n",
            "\\n",
            "Usamos **Gradient Descent** com **Parameter Shift Rule** para calcular gradientes:\\n",
            "\\n",
            "$$\\\\frac{\\\\partial}{\\\\partial \\\\theta_i} f(\\\\theta) = \\\\frac{1}{2}[f(\\\\theta + \\\\pi/2 \\\\hat{e}_i) - f(\\\\theta - \\\\pi/2 \\\\hat{e}_i)]$$"
        ]
    })
    
    notebook["cells"].append({
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# Inicializar pesos\\n",
            "np.random.seed(42)\\n",
            "weights = pnp.array(np.random.randn(4) * 0.1, requires_grad=True)\\n",
            "\\n",
            "# Otimizador\\n",
            "opt = qml.GradientDescentOptimizer(stepsize=0.1)\\n",
            "\\n",
            "# Treinamento\\n",
            "n_epochs = 50\\n",
            "batch_size = 10\\n",
            "\\n",
            "history = {'loss': [], 'acc_train': [], 'acc_test': []}\\n",
            "\\n",
            "print('Iniciando treinamento...')\\n",
            "print(f'Épocas: {n_epochs}, Batch size: {batch_size}')\\n",
            "print('-' * 60)\\n",
            "\\n",
            "for epoch in range(n_epochs):\\n",
            "    # Mini-batch training\\n",
            "    indices = np.random.choice(len(X_train), batch_size, replace=False)\\n",
            "    X_batch = X_train[indices]\\n",
            "    y_batch = y_train[indices]\\n",
            "    \\n",
            "    # Atualizar pesos\\n",
            "    weights = opt.step(lambda w: cost_function(w, X_batch, y_batch), weights)\\n",
            "    \\n",
            "    # Calcular métricas (a cada 5 épocas)\\n",
            "    if (epoch + 1) % 5 == 0:\\n",
            "        loss = cost_function(weights, X_train, y_train)\\n",
            "        acc_train = accuracy(weights, X_train, y_train)\\n",
            "        acc_test = accuracy(weights, X_test, y_test)\\n",
            "        \\n",
            "        history['loss'].append(loss)\\n",
            "        history['acc_train'].append(acc_train)\\n",
            "        history['acc_test'].append(acc_test)\\n",
            "        \\n",
            "        print(f'Época {epoch+1:3d} | Loss: {loss:.4f} | '\n",
            "              f'Acc Treino: {acc_train:.2%} | Acc Teste: {acc_test:.2%}')\\n",
            "\\n",
            "print('-' * 60)\\n",
            "print('✓ Treinamento concluído!')\\n",
            "print(f'\\\\nAcurácia final no teste: {history[\"acc_test\"][-1]:.2%}')"
        ]
    })
    
    # Célula 9: Visualizar resultados
    notebook["cells"].append({
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 7. Visualizar Resultados\\n",
            "\\n",
            "Vamos ver como o VQC aprendeu ao longo do tempo!"
        ]
    })
    
    notebook["cells"].append({
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# Gráfico de aprendizado\\n",
            "fig, axes = plt.subplots(1, 2, figsize=(14, 5))\\n",
            "\\n",
            "# Loss\\n",
            "axes[0].plot(range(5, n_epochs+1, 5), history['loss'], 'o-', linewidth=2)\\n",
            "axes[0].set_xlabel('Época')\\n",
            "axes[0].set_ylabel('Loss (MSE)')\\n",
            "axes[0].set_title('Evolução da Função de Custo')\\n",
            "axes[0].grid(True, alpha=0.3)\\n",
            "\\n",
            "# Acurácia\\n",
            "axes[1].plot(range(5, n_epochs+1, 5), history['acc_train'], \\n",
            "             'o-', label='Treino', linewidth=2)\\n",
            "axes[1].plot(range(5, n_epochs+1, 5), history['acc_test'], \\n",
            "             's-', label='Teste', linewidth=2)\\n",
            "axes[1].set_xlabel('Época')\\n",
            "axes[1].set_ylabel('Acurácia')\\n",
            "axes[1].set_title('Evolução da Acurácia')\\n",
            "axes[1].legend()\\n",
            "axes[1].grid(True, alpha=0.3)\\n",
            "\\n",
            "plt.tight_layout()\\n",
            "plt.show()"
        ]
    })
    
    # Célula 10: Fronteira de decisão
    notebook["cells"].append({
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "## 8. Visualizar Fronteira de Decisão\\n",
            "\\n",
            "Vamos ver como o VQC separa as duas classes no espaço de features."
        ]
    })
    
    notebook["cells"].append({
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [
            "# Criar grid para visualização\\n",
            "x_min, x_max = X[:, 0].min() - 0.5, X[:, 0].max() + 0.5\\n",
            "y_min, y_max = X[:, 1].min() - 0.5, X[:, 1].max() + 0.5\\n",
            "xx, yy = np.meshgrid(np.linspace(x_min, x_max, 50),\\n",
            "                     np.linspace(y_min, y_max, 50))\\n",
            "\\n",
            "# Predizer para cada ponto do grid\\n",
            "print('Calculando fronteira de decisão...')\\n",
            "Z = np.array([quantum_circuit(weights, np.array([x, y])) \\n",
            "              for x, y in zip(xx.ravel(), yy.ravel())])\\n",
            "Z = Z.reshape(xx.shape)\\n",
            "\\n",
            "# Plotar\\n",
            "plt.figure(figsize=(10, 8))\\n",
            "plt.contourf(xx, yy, Z, levels=20, cmap='RdBu', alpha=0.6)\\n",
            "plt.colorbar(label='Saída do VQC')\\n",
            "plt.contour(xx, yy, Z, levels=[0], colors='black', linewidths=2)\\n",
            "\\n",
            "# Plotar pontos de teste\\n",
            "plt.scatter(X_test[y_test==1, 0], X_test[y_test==1, 1], \\n",
            "            c='blue', label='Classe +1 (teste)', \\n",
            "            edgecolors='k', s=100, alpha=0.8)\\n",
            "plt.scatter(X_test[y_test==-1, 0], X_test[y_test==-1, 1], \\n",
            "            c='red', label='Classe -1 (teste)', \\n",
            "            edgecolors='k', s=100, alpha=0.8)\\n",
            "\\n",
            "plt.xlabel('Feature 1')\\n",
            "plt.ylabel('Feature 2')\\n",
            "plt.title(f'Fronteira de Decisão do VQC (Acurácia: {history[\"acc_test\"][-1]:.2%})')\\n",
            "plt.legend()\\n",
            "plt.grid(True, alpha=0.3)\\n",
            "plt.show()\\n",
            "\\n",
            "print('✓ Visualização concluída!')"
        ]
    })
    
    # Célula 11: Conclusão
    notebook["cells"].append({
        "cell_type": "markdown",
        "metadata": {},
        "source": [
            "---\\n",
            "\\n",
            "## 9. Conclusões e Próximos Passos\\n",
            "\\n",
            "### 🎉 Parabéns!\\n",
            "\\n",
            "Você acabou de treinar seu primeiro Variational Quantum Classifier!\\n",
            "\\n",
            "### O que aprendemos:\\n",
            "\\n",
            "✓ VQCs são algoritmos híbridos quântico-clássicos\\n",
            "✓ Codificam dados em estados quânticos\\n",
            "✓ Usam circuitos parametrizados para processamento\\n",
            "✓ Treinam via gradient descent (parameter shift rule)\\n",
            "✓ Podem resolver problemas de classificação não-linear\\n",
            "\\n",
            "### 🚀 Próximos Passos\\n",
            "\\n",
            "Agora que você entende o básico, explore:\\n",
            "\\n",
            "1. **Notebook 02**: [Demonstração de Ruído Benéfico](02_beneficial_noise_demo.ipynb)\\n",
            "   - Como ruído quântico pode **melhorar** o desempenho\\n",
            "   \\n",
            "2. **Notebook 03**: [Framework Completo](03_reproducao_experimentos.ipynb)\\n",
            "   - Experimentos rigorosos com 8,280 configurações\\n",
            "   - Análise estatística completa (QUALIS A1)\\n",
            "   - Múltiplas arquiteturas e tipos de ruído\\n",
            "\\n",
            "### 📚 Referências para Aprofundamento\\n",
            "\\n",
            "- **Schuld, M. et al.** (2020): *Supervised learning with quantum computers*, Springer\\n",
            "- **Cerezo, M. et al.** (2021): *Variational quantum algorithms*, Nature Reviews Physics\\n",
            "- **PennyLane Documentation**: https://pennylane.ai/\\n",
            "\\n",
            "---\\n",
            "\\n",
            "### 💬 Feedback\\n",
            "\\n",
            "Este notebook foi útil? Tem sugestões de melhoria?\\n",
            "Abra uma issue no [repositório GitHub](https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers)!\\n",
            "\\n",
            "---"
        ]
    })
    
    # Salvar notebook
    output_path = Path(__file__).parent / "notebooks" / "01_introducao_vqc.ipynb"
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(notebook, f, indent=1, ensure_ascii=False)
    
    print(f"✓ Notebook criado com sucesso: {output_path}")
    print(f"  - {len(notebook['cells'])} células")

if __name__ == "__main__":
    create_intro_notebook()
