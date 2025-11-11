## Computational Complexity e Limitações Computacionais

Este documento resume a análise de complexidade computacional do framework e fornece recomendações para escalabilidade e reprodução científica.

### 1. Notação e premissas

- n: número de qubits do circuito
- d: profundidade (n_camadas)
- m: número de medidas/observáveis por execução
- Para estados mistos simulados (PennyLane default.mixed) o custo de memória cresce com o espaço de densidade: O(2^{2n}).

### 2. Complexidade temporal

- Evolução de um circuito parametrizado (simulador de estado misto): custo por camada ≈ O(2^{2n}) operações elementares; com d camadas → O(d · 2^{2n}).
- Treinamento (épocas E, dataset T amostras): O(E · T · d · 2^{2n}).

Exemplo: n=4 (como no trabalho), d=2, E=15, T≈300 → custo computacional compatível com execução em CPU com 16GB RAM.

### 3. Complexidade de memória

- Estados puros (vetor): O(2^{n}) memória.
- Estados mistos (matriz densidade): O(2^{2n}) memória — principal limitador para n>6.

### 4. Recomendações para escalabilidade

1. Para n>6, considerar aproximações:
   - Tensor networks / MPS (Matrix Product States) para circuitos com baixa entropia de corte.
   - Técnicas de truncamento (SVD threshold) para manter custo e memória controlados.
2. Subamostragem e curriculum training: iniciar com E baixo e aumentar apenas em trials promissores (pruning adaptativo já implementado com Optuna ajuda nisso).
3. Offload e hospedagem: usar instâncias com mais memória (64–128 GB) ou GPUs compatíveis quando possível (p.ex. PennyLane Lightning + GPU em cloud).

### 5. Medidas práticas adotadas no projeto

- Limitação consciente a 4 qubits para cobertura exaustiva (8,280 experimentos) e permitir 300+ amostras por dataset.
- Uso de `VQC_QUICK` para reduzir épocas em testes exploratórios.

### 6. O que documentar no artigo

- Declarar explicitamente: uso de estados mistos (default.mixed), limites de memória e motivos da escolha n=4.
- Incluir tabela com custo estimado para n ∈ {4,6,8} e recomendações de hardware.

---

Arquivo gerado automaticamente para submissão Qualis A1 — incluir link no README.md (seção "Computational Complexity").
