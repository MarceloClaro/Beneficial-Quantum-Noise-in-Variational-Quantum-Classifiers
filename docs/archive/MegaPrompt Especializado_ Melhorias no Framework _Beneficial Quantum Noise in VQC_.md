# MegaPrompt Especializado: Melhorias no Framework "Beneficial Quantum Noise in VQC"

## 🎯 OBJETIVO GERAL

Refatorar e expandir o framework "Beneficial Quantum Noise in VQC" para alcançar o mais alto padrão de rigor matemático, reprodutibilidade e auditabilidade, garantindo conformidade com os critérios de periódicos Qualis A1 como Nature, Quantum e Physical Review.

---


## 📖 ÍNDICE

### PARTE I: Configuração e Planejamento
1. [Seção 0: Configuração do Projeto de Melhoria](#secao-0)
2. [Glossário de Melhorias](#glossario)


### PARTE II: Execução das Melhorias (10 tarefas)
3. [Tarefa 1: Documentação Matemática Formal](#tarefa-1)
4. [Tarefa 2: Validação Matemática dos Operadores de Kraus](#tarefa-2)
5. [Tarefa 3: Derivação Formal do QNG](#tarefa-3)
6. [Tarefa 4: Centralização e Documentação de Seeds](#tarefa-4)
7. [Tarefa 5: Geração de Manifesto de Execução](#tarefa-5)
8. [Tarefa 6: Correção de Bonferroni nos Testes Post-Hoc](#tarefa-6)
9. [Tarefa 7: Análise de Poder Estatístico](#tarefa-7)
10. [Tarefa 8: Geração de Tabela Código→Método](#tarefa-8)
11. [Tarefa 9: Integração com Cirq e Qiskit](#tarefa-9)
12. [Tarefa 10: Geração de Diagramas de Circuitos](#tarefa-10)


### PARTE III: Validação e Entrega
13. [Checklist de Conformidade Qualis A1](#checklist)
14. [Entrega Final (Pull Request)](#entrega)


---


<a name="secao-0"></a>

## PARTE I: Configuração e Planejamento

### Seção 0: Configuração do Projeto de Melhoria

**Instrução:** Clone o repositório e crie um novo branch para as melhorias.


```bash

# 1. Clone o repositório
git clone <https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers.git>
cd Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers

# 2. Crie um novo branch para as melhorias
git checkout -b feature/qualis-a1-improvements

# 3. Crie um diretório para os novos módulos
mkdir qualis_a1_modules

```text

**Arquivo de Configuração `qai_config.json`:**

Crie um arquivo `qai_config.json` na raiz do projeto para centralizar as configurações das melhorias.

```json
{
  "version": "8.0-QAI",
  "strict_mode": true,
  "default_seed": 42,
  "bonferroni_alpha": 0.05,
  "power_analysis_target": 0.80,
  "output_dir": "./resultados_qualis_a1",
  "enabled_frameworks": ["pennylane", "cirq", "qiskit"]
}

```text

---


<a name="glossario"></a>

### Glossário de Melhorias

- **Documentação Matemática Formal:** Inclusão de equações LaTeX nos docstrings para descrever a teoria por trás de cada componente.
- **Validação Matemática:** Adição de testes unitários que verificam propriedades matemáticas fundamentais (ex: completude dos operadores de Kraus).
- **Manifesto de Execução:** Arquivo JSON que registra todas as configurações, versões de bibliotecas e comandos de uma execução para garantir reprodutibilidade.
- **Correção de Bonferroni:** Método estatístico para ajustar p-valores em comparações múltiplas, controlando a taxa de erro tipo I.
- **Análise de Poder Estatístico:** Cálculo da probabilidade de detectar um efeito real, garantindo que o estudo tem tamanho amostral suficiente.
- **Tabela Código→Método:** Mapeamento explícito entre seções do artigo e linhas de código, para auditoria.


---


<a name="tarefa-1"></a>

## PARTE II: Execução das Melhorias (10 tarefas)

### Tarefa 1: Documentação Matemática Formal

**Objetivo:** Adicionar docstrings com equações LaTeX a todas as 11 classes de ruído.


**Instrução:** Para cada classe de ruído em `framework_investigativo_completo.py` (ex: `RuidoDepolarizante`, `RuidoAmplitudeDamping`), adicione um docstring detalhado com a descrição matemática, os operadores de Kraus e as referências (Nielsen & Chuang).


**Exemplo (para `RuidoDepolarizante`):**

```python
class RuidoDepolarizante(ModeloRuido):
    """
    Modelo de ruído despolarizante (depolarizing noise).
    
    Descrição Matemática:

    ---------------------

    O canal despolarizante é definido por:
    
    $$\\mathcal{E}(\\rho) = (1-p)\\rho + \\frac{p}{3}(X\\rho X + Y\\rho Y + Z\\rho Z)$$
    
    Operadores de Kraus:

    -------------------

    $$K_0 = \\sqrt{1-p} \\mathbb{I}, K_1 = \\sqrt{p/3} X, K_2 = \\sqrt{p/3} Y, K_3 = \\sqrt{p/3} Z$$
    
    Referências:

    -----------

    Nielsen, M. A., & Chuang, I. L. (2010). Quantum Computation and Quantum Information.
    """

```text

**Critério de Aceitação:** Todas as 11 classes de ruído possuem docstrings com equações LaTeX.


---


<a name="tarefa-2"></a>

### Tarefa 2: Validação Matemática dos Operadores de Kraus

**Objetivo:** Adicionar um método de validação para os operadores de Kraus.


**Instrução:** Crie um novo módulo `qualis_a1_modules/validation.py` e adicione a função `validar_operadores_kraus`.


**`qualis_a1_modules/validation.py`:**

```python
import numpy as np

def validar_operadores_kraus(operadores: list, tol: float = 1e-10) -> bool:
    """
    Valida se os operadores de Kraus satisfazem a condição de completude:
    $$\\sum_{k} K_k^\\dagger K_k = \\mathbb{I}$$
    """
    soma = sum(np.conj(K).T @ K for K in operadores)
    identidade = np.eye(operadores[0].shape[0])
    erro = np.linalg.norm(soma - identidade)
    if erro > tol:
        raise ValueError(f"Operadores de Kraus não satisfazem completude: erro = {erro:.2e}")
    return True

```text

**Integração:** Chame esta função dentro do método `aplicar_ruido` de `LindbladNoiseModel`.


**Critério de Aceitação:** A validação é executada para cada modelo de ruído.


---


<a name="tarefa-3"></a>

### Tarefa 3: Derivação Formal do QNG

**Objetivo:** Adicionar documentação matemática à classe `QNG`.


**Instrução:** Adicione um docstring detalhado à classe `QNG` em `framework_investigativo_completo.py` com a derivação do Quantum Natural Gradient, a definição do tensor métrico de Fubini-Study e as referências (Stokes et al. 2020, Yamamoto 2019).


**Exemplo:**

```python
class QNG:
    """
    Quantum Natural Gradient (QNG) optimizer.
    
    Fundamentação Matemática:

    ------------------------

    O QNG utiliza a métrica de Fubini-Study para definir a direção de descida:
    
    $$\\theta_{t+1} = \\theta_t - \\eta g^{-1}(\\theta_t) \\nabla_{\\theta} L(\\theta_t)$$
    
    onde $g(\\theta)$ é o tensor métrico quântico (Quantum Fisher Information Matrix).
    
    Referências:

    -----------

    Stokes, J., et al. (2020). Quantum Natural Gradient. Quantum, 4, 269.
    """

```text

**Critério de Aceitação:** A classe `QNG` possui docstring com derivação matemática.


---


<a name="tarefa-4"></a>

### Tarefa 4: Centralização e Documentação de Seeds

**Objetivo:** Criar uma função centralizada para configuração de seeds.


**Instrução:** Crie um novo módulo `qualis_a1_modules/reproducibility.py` e adicione a função `configurar_seeds_reprodutiveis`.


**`qualis_a1_modules/reproducibility.py`:**

```python
def configurar_seeds_reprodutiveis(seed: int = 42, verbose: bool = True):
    """
    Configura seeds para NumPy, Python random, PennyLane, Optuna e scikit-learn.
    
    Escolha da Seed:

    ---------------

    Seed 42 é utilizada por convenção (referência a "The Hitchhiker's Guide to the Galaxy").
    """

    # ... implementação completa ...

```text

**Integração:** Chame esta função no início do `main` em `framework_investigativo_completo.py`.


**Critério de Aceitação:** Todas as fontes de aleatoriedade são controladas por uma única função.


---


<a name="tarefa-5"></a>

### Tarefa 5: Geração de Manifesto de Execução

**Objetivo:** Criar um arquivo JSON que documente cada execução.


**Instrução:** No módulo `qualis_a1_modules/reproducibility.py`, adicione a função `gerar_manifesto_execucao`.


**`qualis_a1_modules/reproducibility.py`:**

```python
def gerar_manifesto_execucao(config: dict, pasta_resultados: str):
    """
    Gera manifesto JSON com todas as configurações, versões de bibliotecas e comando de execução.
    """

    # ... implementação completa ...

```text

**Integração:** Chame esta função no final do `main`.


**Critério de Aceitação:** Um arquivo `execution_manifest.json` é gerado para cada execução.


---


<a name="tarefa-6"></a>

### Tarefa 6: Correção de Bonferroni nos Testes Post-Hoc

**Objetivo:** Adicionar correção para múltiplas comparações.


**Instrução:** Modifique a classe `TestesEstatisticosAvancados` para incluir um método `testes_post_hoc_com_correcao` que utilize `statsmodels.stats.multitest.multipletests` com `method='bonferroni'`.


**Exemplo:**

```python
def testes_post_hoc_com_correcao(self, df, grupo_col, metrica_col, metodo_correcao=\'bonferroni\'):
    """
    Realiza testes post-hoc com correção para múltiplas comparações.
    
    Fundamentação:

    -------------

    A correção de Bonferroni ajusta o nível de significância: $\\alpha_{ajustado} = \\alpha / m$
    
    Referências:

    -----------

    Dunn, O. J. (1961). Multiple comparisons among means. JASA, 56(293), 52-64.
    """

    # ... implementação completa ...

```text

**Critério de Aceitação:** Os resultados dos testes post-hoc incluem p-valores ajustados.


---


<a name="tarefa-7"></a>

### Tarefa 7: Análise de Poder Estatístico

**Objetivo:** Adicionar cálculo de poder estatístico (1-β).


**Instrução:** Na classe `TestesEstatisticosAvancados`, adicione um método `analise_poder_estatistico` que calcule o poder com base no tamanho do efeito, tamanho da amostra e nível de significância.


**Exemplo:**

```python
def analise_poder_estatistico(self, effect_size, n_per_group, alpha=0.05):
    """
    Calcula o poder estatístico (1-β) para detectar um efeito de tamanho dado.
    
    Fundamentação:

    -------------

    Poder = $1 - \\beta = P(\\text{rejeitar } H_0 | H_0 \\text{ falsa})$
    
    Referências:

    -----------

    Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences.
    """

    # ... implementação completa ...

```text

**Critério de Aceitação:** O relatório final inclui a análise de poder para os principais efeitos.


---


<a name="tarefa-8"></a>

### Tarefa 8: Geração de Tabela Código→Método

**Objetivo:** Criar um mapeamento explícito entre o artigo e o código.


**Instrução:** Crie um novo módulo `qualis_a1_modules/auditing.py` e adicione a função `gerar_tabela_codigo_metodo`.


**`qualis_a1_modules/auditing.py`:**

```python
def gerar_tabela_codigo_metodo(pasta_resultados: str):
    """
    Gera tabela de rastreabilidade Código→Método para auditoria Qualis A1.
    
    Formato:
    | Componente do Método | Arquivo/Função/Linha | Parâmetros | Artefatos Gerados |
    """

    # ... implementação completa com mapeamento ...

```text

**Integração:** Chame esta função como parte do pipeline de geração de resultados.


**Critério de Aceitação:** Um arquivo `tabela_codigo_metodo.csv` é gerado.


---


<a name="tarefa-9"></a>

### Tarefa 9: Integração com Cirq e Qiskit

**Objetivo:** Aumentar a generalidade do framework.


**Instrução:** Crie dois novos arquivos, `framework_cirq.py` e `framework_qiskit.py`, que implementem as mesmas funcionalidades do framework original, mas utilizando as bibliotecas Cirq e Qiskit, respectivamente.


**Exemplo (`framework_cirq.py`):**

```python
import cirq

def circuito_hardware_efficient_cirq(n_qubits, depth, params):
    """
    Implementação do ansatz Hardware-Efficient em Cirq.
    Equivalente a: framework_investigativo_completo.py:circuito_hardware_efficient:L1675
    """

    # ... implementação completa ...

```text

**Critério de Aceitação:** É possível executar o mesmo experimento com PennyLane, Cirq e Qiskit.


---


<a name="tarefa-10"></a>

### Tarefa 10: Geração de Diagramas de Circuitos

**Objetivo:** Melhorar a didática do artigo com visualizações.


**Instrução:** Crie um novo módulo `qualis_a1_modules/visualization.py` e adicione a função `gerar_diagrama_circuito`.


**`qualis_a1_modules/visualization.py`:**

```python
def gerar_diagrama_circuito(ansatz_func, n_qubits, depth, pasta_resultados):
    """
    Gera diagrama visual do circuito quântico usando PennyLane, Cirq ou Qiskit.
    """

    # ... implementação completa ...

```

**Critério de Aceitação:** Diagramas de todos os 9 ansätze são gerados e salvos como PNG.


---


<a name="checklist"></a>

## PARTE III: Validação e Entrega

### Checklist de Conformidade Qualis A1

#### 1. Rigor Matemático (30 pts)
- [ ] Docstrings com equações LaTeX (10 pts)
- [ ] Validação de operadores de Kraus (10 pts)
- [ ] Derivação do QNG (10 pts)


#### 2. Reprodutibilidade (30 pts)
- [ ] Seeds centralizadas (15 pts)
- [ ] Manifesto de execução (15 pts)


#### 3. Rigor Estatístico (20 pts)
- [ ] Correção de Bonferroni (10 pts)
- [ ] Análise de poder (10 pts)


#### 4. Auditoria e Transparência (20 pts)
- [ ] Tabela Código→Método (10 pts)
- [ ] Integração Cirq/Qiskit (5 pts)
- [ ] Diagramas de circuitos (5 pts)


**Pontuação Final:** [Soma dos pontos] / 100


---


<a name="entrega"></a>

### Entrega Final (Pull Request)

1. ✅ Crie um Pull Request do branch `feature/qualis-a1-improvements` para o `main`.
2. ✅ No corpo do PR, inclua:
   - Resumo das 10 melhorias implementadas.
   - Pontuação final do Checklist de Conformidade Qualis A1.
   - Link para o novo diretório `qualis_a1_modules`.
   - Instruções de como executar o framework aprimorado.
3. ✅ Solicite revisão de pelo menos 2 coautores.
4. ✅ Após aprovação, faça o merge para o `main`.


---


**FIM DO MEGAPROMPT**

