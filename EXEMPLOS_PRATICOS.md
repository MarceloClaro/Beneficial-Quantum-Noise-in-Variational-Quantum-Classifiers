# Exemplos Práticos - Caso Beneficial Quantum Noise

Este documento fornece exemplos concretos de aplicação do MegaPrompt v2.0 ao projeto "Beneficial Quantum Noise in Variational Quantum Classifiers".

## 📊 Visão Geral do Projeto

**Objetivo**: Demonstrar que ruído quântico pode ser benéfico para VQCs  
**Configurações**: 2,688 experimentos  
**Frameworks**: PennyLane, Qiskit, Cirq  
**Datasets**: 4 (moons, circles, blobs, iris)  
**Ansätze**: 7 arquiteturas  
**Tipos de Ruído**: 6 modelos físicos  

## 🔧 Exemplo 1: Configuração do config.json

```json
{
  "output_mode": "MODE_A",
  "reference_policy": "R1",
  "editorial_profile": "PROFILE_PR_QUANTUM",
  "target_journals": {
    "primary": "Quantum",
    "secondary": ["Physical Review A", "Physical Review Research"]
  },
  "inputs": {
    "code_path": "https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers",
    "data_path": "[Gerado pelo código]",
    "artifacts_path": "resultados/"
  },
  "user_inputs": {
    "research_question": "Can quantum noise be systematically leveraged as a beneficial regularizer in Variational Quantum Classifiers?",
    "hypotheses": [
      "H₀: There exists an optimal noise level p* > 0 that improves generalization",
      "H₁: Different noise types exhibit distinct beneficial regimes",
      "H₂: Dynamic noise schedules outperform static configurations"
    ]
  }
}
```

## 📝 Exemplo 2: Tabela de Rastreabilidade

| Seção | Afirmação/Número | Evidência (Arquivo/Função/Linha) | Referência |
|-------|------------------|----------------------------------|------------|
| Abstract | "65.83% accuracy achieved" | `resultados/melhor_config.json:accuracy` | - |
| Methods | "Lindblad master equation" | `framework_investigativo_completo.py:RuidoDepolarizante:L1523-1548` | (Lindblad, 1976) |
| Results | "12.6% faster convergence" | `resultados/comparison_schedules.csv:epochs_mean` | - |
| Results | "Cosine schedule superior" | `resultados/anova_schedules.json:p_value=0.003` | - |
| Methods | "Quantum Natural Gradient" | `framework_investigativo_completo.py:ClassificadorVQC:L2341-2398` | (Stokes et al., 2020) |

## 🧮 Exemplo 3: Cálculo de Configurações

```python
# Fatores experimentais
factors = {
    'dataset': 4,        # moons, circles, blobs, iris
    'ansatz': 7,         # 7 arquiteturas diferentes
    'noise_type': 6,     # 6 tipos de ruído
    'schedule': 2,       # Constant, Linear
    'initialization': 8, # 8 estratégias de inicialização
    'seed': 5           # 5 seeds diferentes
}

# Total de configurações
total = 4 × 7 × 6 × 2 × 8 × 5 = 13,440 configurações

# Ajuste: apenas schedules aplicáveis
# (schedule não se aplica a noise_type='None')
adjusted_total = (4 × 7 × 1 × 1 × 8 × 5) +  # Sem ruído
                 (4 × 7 × 6 × 2 × 8 × 5)    # Com ruído
               = 1,120 + 13,440 = 14,560

# Simplificação para demo rápida
demo_configs = 4 × 1 × 2 × 1 × 1 × 1 = 8
```

## 🔬 Exemplo 4: Mapeamento Código→Método

| Componente do Método | Arquivo/Função/Linha | Parâmetros | Artefatos Gerados |
|---------------------|---------------------|------------|-------------------|
| **Ansatz: Hardware-Efficient** | `framework_investigativo_completo.py:circuito_hardware_efficient:L876-923` | `n_qubits=4, n_camadas=2, pesos` | PennyLane QNode |
| **Ruído: Depolarizing** | `framework_investigativo_completo.py:RuidoDepolarizante:L1523-1548` | `nivel: float ∈ [0,1]` | Operadores de Kraus: K₀, K₁, K₂, K₃ |
| **Otimizador: QNG** | `framework_investigativo_completo.py:ClassificadorVQC:L2341-2398` | `lr=0.01, approx='block-diag'` | Gradientes naturais |
| **Validação: K-Fold** | `framework_investigativo_completo.py:executar_experimento:L2789-2834` | `k=5, stratified=True, seed=42` | Métricas por fold |
| **Análise Estatística** | `qualis_a1_modules/statistical_extensions.py:testes_post_hoc_com_correcao:L152-201` | `method='bonferroni', alpha=0.05` | Tabela de p-valores ajustados |

## 📈 Exemplo 5: Tabela S1 (Amostra)

```csv
config_id,dataset,ansatz,noise_type,noise_level,schedule,initialization,seed,accuracy_train,accuracy_test,loss_final,epochs_to_converge
1,moons,BasicEntangler,None,0.0,Constant,random,42,0.9833,0.9500,0.1234,87
2,moons,BasicEntangler,Depolarizing,0.001,Constant,random,42,0.9750,0.9667,0.1156,92
3,moons,BasicEntangler,Depolarizing,0.005,Constant,random,42,0.9667,0.9750,0.1089,78
4,moons,BasicEntangler,Depolarizing,0.005,Linear,random,42,0.9667,0.9833,0.0978,65
...
2688,iris,RandomEntangling,GeneralizedAmplitudeDamping,0.1,Cosine,pi_8,1024,0.8900,0.8667,0.3456,145
```

## 📊 Exemplo 6: Análise Estatística

### ANOVA de Dois Fatores

```python
# Hipótese: Tipo de ruído e schedule interagem significativamente

# Resultado:
# F(noise_type) = 12.45, p < 0.001
# F(schedule) = 8.92, p = 0.003
# F(noise_type × schedule) = 3.67, p = 0.027

# Interpretação: Interação significativa confirma H₂
```

### Testes Post-Hoc com Correção de Bonferroni

```python
# Comparações múltiplas entre schedules
# α_original = 0.05
# Número de comparações = 3 (Constant vs Linear, Constant vs Cosine, Linear vs Cosine)
# α_ajustado = 0.05/3 = 0.0167

comparisons = [
    ('Constant', 'Linear', p=0.023, adjusted_p=0.069, significant=False),
    ('Constant', 'Cosine', p=0.004, adjusted_p=0.012, significant=True),
    ('Linear', 'Cosine', p=0.156, adjusted_p=0.468, significant=False)
]

# Resultado: Apenas Constant vs Cosine é significativo após correção
```

### Tamanho de Efeito

```python
# Cohen's d entre melhor configuração com ruído vs sem ruído
# Média com ruído: 0.9583 (DP: 0.0234)
# Média sem ruído: 0.9200 (DP: 0.0312)

d = (0.9583 - 0.9200) / sqrt((0.0234² + 0.0312²)/2)
d = 0.0383 / 0.0276
d = 1.39  # Grande efeito (> 0.8)

# Interpretação: Ruído tem efeito substancial na acurácia
```

## 🎯 Exemplo 7: Problema Formal

```latex
\textbf{Seja:}
\begin{itemize}
    \item $\mathcal{D} = \{(x_i, y_i)\}_{i=1}^{N}$ um dataset com $N$ amostras
    \item $U(\theta)$ um circuito quântico parametrizado por $\theta \in \mathbb{R}^{P}$
    \item $\mathcal{N}_p(\cdot)$ um canal quântico com parâmetro de ruído $p \in [0, p_{max}]$
    \item $L(\theta, p)$ a função de custo (cross-entropy loss)
\end{itemize}

\textbf{O problema de otimização é:}
\begin{equation}
    (\theta^*, p^*) = \arg\min_{\theta, p} \mathbb{E}_{(x,y) \sim \mathcal{D}_{test}} [L(y, f(x; \theta, p))]
\end{equation}

\textbf{Sujeito a:}
\begin{align}
    p &\in [0, 0.5] \\
    |\theta| &= 2 \times n_{qubits} \times n_{layers} \\
    T &\leq 100 \text{ epochs}
\end{align}

\textbf{Hipótese Principal (H₀):}
\begin{equation}
    \exists p^* > 0 : \mathbb{E}[\text{Acc}(f_{\theta^*, p^*})] > \mathbb{E}[\text{Acc}(f_{\theta^*_0, 0})]
\end{equation}
```

## 📝 Exemplo 8: Algorithm 1 em LaTeX

```latex
\begin{algorithm}[H]
\caption{Experimental Pipeline for Beneficial Noise Analysis}
\label{alg:pipeline}
\begin{algorithmic}[1]
\REQUIRE Datasets $\mathcal{D} = \{\mathcal{D}_1, \ldots, \mathcal{D}_4\}$
\REQUIRE Configurations $\mathcal{C} = \{\text{ansatz}, \text{noise}, \text{schedule}, \ldots\}$
\REQUIRE Seeds $S = \{42, 123, 456, 789, 1024\}$
\ENSURE Results table $R$ with accuracy, loss, epochs
\STATE Initialize $R \leftarrow \emptyset$
\FOR{each dataset $D \in \mathcal{D}$}
    \FOR{each configuration $c \in \mathcal{C}$}
        \FOR{each seed $s \in S$}
            \STATE Set random seed: \texttt{np.random.seed}$(s)$
            \STATE Initialize model $M_c$ with configuration $c$
            \STATE Split: $(D_{train}, D_{val}, D_{test}) \leftarrow$ \texttt{train\_test\_split}$(D, s)$
            \STATE $M_c^* \leftarrow$ \texttt{train}$(M_c, D_{train}, D_{val})$
            \STATE $(acc, loss) \leftarrow$ \texttt{evaluate}$(M_c^*, D_{test})$
            \STATE Append $(D, c, s, acc, loss)$ to $R$
        \ENDFOR
    \ENDFOR
\ENDFOR
\STATE Perform statistical analysis: ANOVA, post-hoc, effect sizes
\RETURN $R$
\end{algorithmic}
\end{algorithm}
```

## 🔍 Exemplo 9: Marcadores de Integridade

### Caso: Versão de biblioteca não documentada

**Texto Original:**
> "We used PennyLane for quantum circuit simulation."

**Com Marcador:**
> "We used PennyLane **[INFORMAÇÃO AUSENTE: versão não especificada]** for quantum circuit simulation. Based on installation logs, PennyLane 0.38.0 was likely used."

### Caso: Experimento não executado

**Texto Original:**
> "We also tested on a 10-qubit circuit."

**Com Marcador:**
> "We also **planned** to test on a 10-qubit circuit **[NÃO DISPONÍVEL: experimento não executado devido a restrições computacionais]**."

### Caso: Referência faltando (modo R0)

**Texto Original:**
> "Neural tangent kernel theory suggests..."

**Com Marcador:**
> "Neural tangent kernel theory **[LACUNA DE CITAÇÃO: referência não disponível em modo R0]** suggests..."

## ✅ Exemplo 10: Quality Gates

### Quality Gate F1 (Fase 1)

```markdown
## Quality Gate F1 - Checklist

- [x] Cada item tem origem (arquivo/função/linha/config/log)
  - Exemplo: RuidoDepolarizante → framework_investigativo_completo.py:L1523-1548
- [x] O total de configurações foi calculado e conferido
  - Calculado: 2,688 = 4 × 7 × 6 × 2 × 8 × 1
  - Conferido: Geração de tabela S1 produziu 2,688 linhas
- [x] Ausências explicitadas com [INFORMAÇÃO AUSENTE]/[NÃO DISPONÍVEL]
  - Versões de bibliotecas: Documentadas em requirements.txt
  - Hiperparâmetros: Todos documentados em config
  - Seeds: Fixadas e reportadas (42, 123, 456, 789, 1024)

Status: ✅ APROVADO - Pode prosseguir para Fase 2
```

### Quality Gate F4 (Fase 4)

```markdown
## Quality Gate F4 - Checklist

- [x] Sem números sem lastro
  - Todos os valores numéricos rastreáveis via tabela de rastreabilidade
  - Exemplo: 65.83% → resultados/melhor_config.json:accuracy
- [x] R0 respeitado (se aplicável)
  - Modo R1 ativo: Referências adicionadas com DOI
- [x] Methods completo (notação + equações + algoritmo + mapa)
  - Notação: Seção 2.1 (símbolos definidos)
  - Equações: 15 equações LaTeX numeradas
  - Algorithm 1: Incluído em Methods
  - Mapa código→método: Tabela 3 páginas

Status: ✅ APROVADO - Pode prosseguir para Fase 5
```

## 🎓 Exemplo 11: Resultado Final

### Pontuação Audit Checklist

```
Categoria                 Pontos    Max    %
─────────────────────────────────────────────
Reprodutibilidade           30      30   100%
Rastreabilidade             30      30   100%
Rigor Estatístico           20      20   100%
Transparência               20      20   100%
─────────────────────────────────────────────
TOTAL                      100     100   100%

Status: ✅ EXCELENTE - Pronto para submissão
```

### Índice de Consistência

```
Verificações Totais: 247
Verificações Aprovadas: 242
Problemas Encontrados: 5 (menores)
Avisos/Marcadores: 3

Índice de Consistência: 98.0%

Status: ✅ EXCELENTE - Meta atingida (≥95%)
```

---

**Nota**: Estes exemplos são baseados no projeto real "Beneficial Quantum Noise in VQC" e demonstram a aplicação prática do MegaPrompt v2.0. Adapte os valores, arquivos e configurações ao seu projeto específico.
