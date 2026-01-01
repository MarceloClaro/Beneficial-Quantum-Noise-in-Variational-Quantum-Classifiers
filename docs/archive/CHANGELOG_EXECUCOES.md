# Changelog de Execuções Multiframework

**Propósito:** Rastrear todas as execuções experimentais, mudanças de configuração e run_ids "golden" (referência).


---


## Formato de Entrada

Cada execução deve seguir este formato:

```markdown

### Run ID: YYYYMMDD_XXX_<framework/all>

**Data:** YYYY-MM-DD HH:MM  
**Frameworks:** PennyLane | Qiskit | Cirq | Todos  
**Tipo:** Baseline | Experimento | Validação | Ablation  
**Status:** ✅ Sucesso | ⚠️ Parcial | ❌ Falha


#### Configuração:
- Config file: `configs/experiment_unified.yaml` (v1.0)
- Seeds: [42, 43]
- Datasets: [iris, wine]
- Noise models: [depolarizing, phase_damping]
- Total configs: 128


#### Mudanças em Relação ao Run Anterior:
- Adicionado schedule "cosine"
- Corrigido bug em TREX calibration
- Atualizado PennyLane 0.33.0 → 0.33.1


#### Resultados Principais:
- Accuracy média: 0.87 ± 0.03
- Melhor config: iris + phase_damping + cosine
- Tempo total: 2h 15min


#### Artefatos Gerados:
- `results/<framework>/20251227_001/metrics.csv`
- `results/comparisons/20251227_001/comparative_table.csv`
- `figures/20251227_001/*.png`
- `manifests/<framework>/20251227_001/manifest.json`


#### Notas:
- Framework Qiskit 12% mais lento devido a transpilation
- Cirq density matrix mode requer 2GB RAM por config
- Algumas configs não convergiram (registrado em logs)


**Golden Run:** ❌ Não | ✅ Sim (se for referência para artigo)

```text

---


## Histórico de Execuções

### Run ID: 20251227_000_baseline

**Data:** 2025-12-27 10:00 (PLANEJADO)  
**Frameworks:** Todos (PennyLane, Qiskit, Cirq)  
**Tipo:** Baseline  
**Status:** 🔄 Planejado


#### Configuração:
- Config file: `configs/experiment_unified.yaml` (v1.0)
- Seeds: [42, 43]
- Datasets: [iris, wine, digits, breast_cancer]
- Noise models: [depolarizing, amplitude_damping, phase_damping]
- Ansätze: [simplified_two_local, hardware_efficient, strongly_entangling]
- Schedules: [static, cosine]
- Total configs: 4 × 3 × 3 × 2 × 2 = 144 configs × 3 frameworks = 432 runs


#### Mudanças em Relação ao Run Anterior:
- Primeira execução multiframework completa
- Implementação de config unificado
- Geração de manifestos habilitada


#### Resultados Principais:
- [A SER PREENCHIDO APÓS EXECUÇÃO]


#### Artefatos Gerados:
- [A SER PREENCHIDO APÓS EXECUÇÃO]


#### Notas:
- Esta é a execução baseline de referência
- Todos os 3 frameworks devem ser executados com mesma config
- Resultados servirão como benchmark para futuras melhorias


**Golden Run:** ✅ Sim (quando completado com sucesso)


---


## Run IDs de Referência

### Golden Runs (Referência para Artigo)

| Run ID | Data | Descrição | Seção do Artigo | Status |
|--------|------|-----------|-----------------|--------|
| 20251227_000_baseline | 2025-12-27 | Baseline multiframework completo | 4.10 (Resultados Multiframework) | 🔄 Planejado |
| TBD | TBD | Validação com TREX+AUEC | 4.10, 5.8.4 (Sinergia) | 🔄 Futuro |
| TBD | TBD | Ablation study schedules | 5.7 (Discussão Schedules) | 🔄 Futuro |

---


## Tracking de Mudanças de Configuração

### v1.0 → v1.1 (Planejado)

#### Mudanças:
- [ ] Adicionar TREX error mitigation habilitado
- [ ] Adicionar AUEC framework habilitado
- [ ] Aumentar seeds para [42, 43, 44]
- [ ] Adicionar dataset "synthetic" para controle


#### Impacto Esperado:
- +30% tempo de execução (overhead TREX/AUEC)
- +50% uso de memória (AUEC Kalman filter state)
- Melhoria esperada: +6-7% accuracy


**Decisão:** Pendente validação do baseline


---


## Troubleshooting e Issues

### Issue #1: Qiskit Transpilation Overhead (2025-12-27)

**Problema:** Transpilation adiciona 20-30% de portas extras  
**Impacto:** Tempo de execução +15-20%  
**Solução:** Documentado em `docs/equivalencias_e_limitacoes.md` como limitação conhecida  
**Status:** ✅ Documentado, ⚠️ Não mitigado (aceitável)


### Issue #2: Cirq Density Matrix RAM (2025-12-27)

**Problema:** Density matrix simulator usa 4× mais RAM que state vector  
**Impacto:** Limite de ~12 qubits vs. ~20 qubits  
**Solução:** Para n_qubits=4, aceitável (1-2GB RAM)  
**Status:** ✅ Validado OK para este projeto


---


## Template para Nova Execução

```markdown

### Run ID: YYYYMMDD_XXX_<descrição>

**Data:** YYYY-MM-DD HH:MM  
**Frameworks:** [Lista]  
**Tipo:** [Baseline | Experimento | Validação | Ablation]  
**Status:** [Status]


#### Configuração:
- Config file: `path/to/config.yaml` (vX.X)
- [Listar parâmetros principais]


#### Mudanças em Relação ao Run Anterior:
- [Listar mudanças]


#### Resultados Principais:
- [Após execução]


#### Artefatos Gerados:
- [Após execução]


#### Notas:
- [Observações importantes]


**Golden Run:** [Sim/Não]

```

---


**Última Atualização:** 2025-12-27  
**Mantido por:** Equipe de Desenvolvimento  
**Reviewer:** @MarceloClaro

