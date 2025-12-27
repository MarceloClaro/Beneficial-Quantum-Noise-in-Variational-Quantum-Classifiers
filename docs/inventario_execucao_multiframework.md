# Inventário de Execução Multiframework

**Data:** 2025-12-27  
**Status:** 🚧 EM ANDAMENTO  
**Versão:** 1.0

---

## 📋 Visão Geral

Este documento mapeia todos os componentes do repositório necessários para execução multiframework (PennyLane, Qiskit, Cirq) com rastreabilidade completa.

---

## 🗂️ Estrutura do Repositório

### Diretórios Principais

```
.
├── configs/                    # Configurações experimentais unificadas
├── docs/                       # Documentação
│   └── execution_plan/         # Planos de execução
├── logs/                       # Logs de execução por framework
│   ├── pennylane/
│   ├── qiskit/
│   └── cirq/
├── results/                    # Resultados por framework
│   ├── pennylane/
│   ├── qiskit/
│   ├── cirq/
│   └── comparisons/           # Análises comparativas
├── figures/                    # Figuras geradas
├── manifests/                  # Manifestos de execução
│   ├── pennylane/
│   ├── qiskit/
│   └── cirq/
├── tests/                      # Testes unitários e smoke tests
├── artigo_cientifico/         # Documentação do artigo
└── tools/                      # Ferramentas de validação
```

---

## 🔑 Arquivos-Chave Identificados

### Framework Implementations

| Framework | Arquivo Principal | Executor | Status |
|-----------|------------------|----------|--------|
| **PennyLane (Baseline)** | `framework_investigativo_completo.py` (5,661 linhas) | Built-in | ✅ Funcional |
| **Qiskit** | `framework_qiskit.py` | `executar_framework_qiskit.py` | ✅ Funcional |
| **Cirq** | `framework_cirq.py` | `executar_framework_cirq.py` | ✅ Funcional |

### Comparison & Analysis Tools

| Ferramenta | Arquivo | Função | Status |
|------------|---------|--------|--------|
| **Multiframework Executor** | `executar_multiframework.py` | Executa os 3 frameworks | ✅ Existe |
| **Quick Executor** | `executar_multiframework_rapido.py` | Execução rápida | ✅ Existe |
| **Comparative Analysis** | `comparacao_multiframework_completa.py` | Análise comparativa | ✅ Existe |
| **Results Generator** | `generate_comparative_results.py` | Gera tabelas comparativas | ✅ Existe |

### Support Components

| Componente | Arquivo | Função | Status |
|------------|---------|--------|--------|
| **TREX Error Mitigation** | `trex_error_mitigation.py` | Correção de erros de leitura | ✅ Implementado |
| **AUEC Framework** | `adaptive_unified_error_correction.py` | Correção adaptativa unificada | ✅ Implementado |
| **Visualization** | `visualize_results.py` | Geração de figuras | ✅ Existe |

---

## 🔬 Componentes Experimentais

### PennyLane (Baseline)

**Arquivo:** `framework_investigativo_completo.py`

**Componentes Principais:**
- **11 Modelos de Ruído:**
  1. `RuidoDepolarizante` - Depolarizing noise
  2. `RuidoAmplitudeDamping` - Amplitude damping (T₁)
  3. `RuidoPhaseDamping` - Phase damping (T₂)
  4. `RuidoBitFlip` - Bit flip errors
  5. `RuidoPhaseFlip` - Phase flip errors
  6. `RuidoGeneralizado` - Generalized noise
  7. `RuidoAmortecimentoThermal` - Thermal damping
  8. `RuidoCruzado` - Crosstalk
  9. `RuidoDrift` - Parameter drift
  10. `RuidoCombinado` - Combined noise models
  11. `RuidoDinamico` - Dynamic noise schedules

- **7 Ansätze (Arquiteturas de Circuito):**
  1. `SimplifiedTwoLocal` - 2-local simplified
  2. `HardwareEfficientAnsatz` - Hardware-efficient
  3. `StronglyEntanglingLayers` - Strong entanglement
  4. `AmplitudeEmbedding` - Amplitude encoding
  5. `AngleEmbedding` - Angle encoding
  6. `IQPEmbedding` - IQP embedding
  7. `BasicEntanglerLayers` - Basic entangling

- **4 Schedules de Ruído:**
  1. Static (constante)
  2. Linear (decrescente linear)
  3. Exponential (decrescente exponencial)
  4. Cosine (decrescente cosenoidal)

- **5 Estratégias de Inicialização:**
  1. Random
  2. Xavier/Glorot
  3. He/Kaiming
  4. Identity-preserving
  5. Low-variance

**Fatores Experimentais:**
- Datasets: 4 (Iris, Wine, Digits, Breast Cancer)
- Modelos de Ruído: 7 principais
- Ansätze: 6 testados
- Schedules: 2 (Static, Cosine)
- Seeds: 8 [42, 43]

**Total de Configurações:** 4 × 7 × 6 × 2 × 8 = 2,688 configurações

**Métricas Coletadas:**
- Accuracy (treino, validação, teste)
- Loss (treino, validação)
- Tempo de execução
- Número de portas
- Profundidade do circuito
- Gradiente médio (para detectar barren plateaus)

### Qiskit

**Arquivo:** `framework_qiskit.py`

**Simulador:** Aer (state vector + noise models)

**Equivalências com PennyLane:**
- Ansätze traduzidos usando `QuantumCircuit`
- Noise models via `NoiseModel` do Qiskit Aer
- Otimização via COBYLA/SPSA

**Diferenças Conhecidas:**
- Backend específico: `AerSimulator`
- Noise model: Usa `depolarizing_error()` do Qiskit
- Transpilation: Pode adicionar overhead de portas

**Métricas Específicas:**
- Transpiled circuit depth
- Gate count após otimização
- Simulation time (pode ser maior devido ao overhead)

### Cirq

**Arquivo:** `framework_cirq.py`

**Simulador:** DensityMatrixSimulator / Simulator

**Equivalências com PennyLane:**
- Ansätze traduzidos usando `cirq.Circuit`
- Noise via `cirq.depolarize()`, `cirq.amplitude_damp()`, etc.
- Otimização via SciPy optimizers

**Diferenças Conhecidas:**
- Qubit naming: `cirq.GridQubit` vs. wires
- Noise application: Per-gate basis
- Measurement: Projective measurements

**Métricas Específicas:**
- Circuit moment depth (diferente de gate depth)
- Fidelity tracking
- Execution time (geralmente mais rápido que Qiskit)

---

## ⚙️ Parametrização e Configuração

### Formato Atual

**PennyLane:** Constantes internas no `framework_investigativo_completo.py`
- Seeds: Linha ~1250
- Hiperparâmetros: Linha ~890-950
- Datasets: Linha ~320-480

**Qiskit:** Parâmetros em `executar_framework_qiskit.py`
- CLI arguments ou valores hard-coded

**Cirq:** Parâmetros em `executar_framework_cirq.py`
- CLI arguments ou valores hard-coded

### Formato Alvo (Unificado)

**Arquivo:** `configs/experiment_unified.yaml` (a ser criado)

```yaml
version: "1.0"
run_id: "20251227_multiframework_baseline"

# Seeds para reprodutibilidade
seeds:
  global: 42
  frameworks:
    pennylane: [42, 43]
    qiskit: [42, 43]
    cirq: [42, 43]

# Fatores experimentais
experimental_design:
  datasets: ["iris", "wine", "digits", "breast_cancer"]
  noise_models: ["depolarizing", "amplitude_damping", "phase_damping"]
  ansatze: ["simplified_two_local", "hardware_efficient", "strongly_entangling"]
  schedules: ["static", "cosine"]
  
# Parâmetros do ansatz
ansatz_params:
  n_qubits: 4
  depth: 2
  entanglement: "full"

# Noise model
noise_params:
  gamma_range: [0.0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2]
  
# Otimizador
optimizer:
  name: "Adam"
  learning_rate: 0.01
  epochs: 100
  batch_size: 32
  
# Métricas a coletar
metrics:
  - accuracy_train
  - accuracy_val
  - accuracy_test
  - loss_train
  - loss_val
  - execution_time
  - circuit_depth
  - gate_count
  - gradient_mean
  - gradient_variance

# Critérios de convergência
convergence:
  patience: 10
  min_delta: 0.001
```

---

## 📊 Métricas e Schema de Dados

### Schema de `metrics.csv`

```csv
framework,run_id,dataset,seed,config_id,noise_model,ansatz,schedule,gamma,accuracy_train,accuracy_val,accuracy_test,loss_train,loss_val,execution_time,circuit_depth,gate_count,gradient_mean,gradient_var
pennylane,20251227_001,iris,42,cfg_001,depolarizing,simplified_two_local,static,0.005,0.95,0.93,0.91,0.12,0.15,12.3,8,45,0.023,0.0012
```

### Schema de `summary.csv`

```csv
framework,dataset,noise_model,ansatz,schedule,accuracy_mean,accuracy_std,accuracy_ci_lower,accuracy_ci_upper,n_samples
pennylane,iris,depolarizing,simplified_two_local,static,0.91,0.02,0.89,0.93,8
```

### Schema de `comparative_table.csv`

```csv
condition,pennylane_mean,pennylane_std,qiskit_mean,qiskit_std,cirq_mean,cirq_std,anova_f,anova_p,effect_size_cohens_d
iris_depolarizing_static,0.91,0.02,0.89,0.03,0.90,0.02,2.34,0.12,0.45
```

---

## 🧪 Testes e Validação

### Smoke Tests (a implementar)

**Arquivo:** `tests/smoke_multiframework.py`

**Testes Necessários:**
1. ✅ Importação de cada framework sem erro
2. ✅ Execução mínima (1 configuração) por framework
3. ✅ Geração de outputs no schema correto
4. ✅ Leitura e parse de arquivos de configuração
5. ✅ Validação de seeds (determinismo)

### Testes Unitários Existentes

**Diretório:** `tests/`

**Cobertura:**
- ⚠️ Limitada - necessita expansão
- ⚠️ Sem testes específicos para equivalência entre frameworks

---

## 🔄 Fluxo de Execução

### Ordem Mandatória

1. **Preparação:**
   - Criar `configs/experiment_unified.yaml`
   - Instalar dependências (requirements.txt)
   - Validar ambiente (Python, libs)

2. **Execução PennyLane (Baseline):**
   ```bash
   python framework_investigativo_completo.py --config configs/experiment_unified.yaml --output results/pennylane/run_001
   ```

3. **Execução Qiskit:**
   ```bash
   python executar_framework_qiskit.py --config configs/experiment_unified.yaml --output results/qiskit/run_001
   ```

4. **Execução Cirq:**
   ```bash
   python executar_framework_cirq.py --config configs/experiment_unified.yaml --output results/cirq/run_001
   ```

5. **Análise Comparativa:**
   ```bash
   python generate_comparative_results.py --run-id run_001 --output results/comparisons/run_001
   ```

6. **Geração de Figuras:**
   ```bash
   python visualize_results.py --run-id run_001 --output figures/run_001
   ```

---

## 📝 Documentação a Atualizar

### Arquivos do Artigo

**Diretório:** `artigo_cientifico/fase4_secoes/`

| Arquivo | Seções Afetadas | Tipo de Atualização |
|---------|-----------------|---------------------|
| `metodologia_completa.md` | 3.2, 3.3 | Adicionar detalhes multiframework |
| `resultados_completo.md` | 4.10 | Atualizar com run_id dos resultados |
| `discussao_completa.md` | 5.8 | Interpretar comparações |

### Documentos de Rastreabilidade (a criar)

1. **`docs/equivalencias_e_limitacoes.md`**
   - Diferenças inevitáveis entre frameworks
   - Justificativas e impactos esperados

2. **`docs/melhorias_map.md`**
   - Mapeamento de melhorias do MegaPrompt
   - Status de implementação

3. **`CHANGELOG_EXECUCOES.md`**
   - Histórico de execuções (run_ids)
   - Mudanças entre versões

4. **`relatorio_consistencia.md`**
   - Auditoria código ↔ dados ↔ documentação
   - Discrepâncias e resoluções

---

## 🎯 Próximas Ações (Roadmap)

### Fase 1: Inventário ✅ (Este Documento)
- [x] Mapear estrutura do repositório
- [x] Identificar arquivos-chave
- [x] Documentar componentes experimentais
- [x] Definir schema de dados

### Fase 2: Configuração Unificada (Em Andamento)
- [ ] Criar `configs/experiment_unified.yaml`
- [ ] Adaptar scripts para ler config unificada
- [ ] Validar paridade de configurações

### Fase 3: Melhorias do MegaPrompt
- [ ] Extrair requisitos do MegaPrompt
- [ ] Criar `docs/melhorias_map.md`
- [ ] Implementar melhorias prioritárias

### Fase 4: Testes e Validação
- [ ] Criar `tests/smoke_multiframework.py`
- [ ] Validar determinismo (seeds)
- [ ] Testar schema de outputs

### Fases 5-9: Execução e Análise (Requer Ambiente Computacional)
- [ ] Executar PennyLane, Qiskit, Cirq
- [ ] Gerar análises comparativas
- [ ] Criar figuras
- [ ] Atualizar documentação

---

## 📞 Contato e Tracking

**Documento de Referência:** Este arquivo  
**Atualizações:** A cada conclusão de fase  
**Issues/Pendências:** Registrar em `docs/execution_plan/issues.md`

---

**Última Atualização:** 2025-12-27  
**Responsável:** Equipe de Desenvolvimento  
**Reviewer:** @MarceloClaro
