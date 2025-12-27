# Multiframework Execution Roadmap

**Date:** December 27, 2025  
**Status:** 🚧 IMPLEMENTATION PLAN  
**Based on:** MegaPrompt Operational Plan from @MarceloClaro

---

## 📋 Executive Summary

This document outlines the roadmap for implementing the comprehensive multiframework execution plan across PennyLane, Qiskit, and Cirq, following the 13-section MegaPrompt operational specification.

### Scope

Execute and compare implementations across three quantum frameworks, applying improvements from "MegaPrompt Especializado_ Melhorias no Framework _Beneficial Quantum Noise in VQC_.md", and update all artifacts with full traceability and reproducibility.

---

## 🎯 Key Objectives

1. **Functional Parity**: Ensure identical experimental hypotheses, factors, levels, and metrics across all three frameworks
2. **Complete Execution**: Run full pipeline on PennyLane (baseline), Qiskit, and Cirq
3. **Comparative Analysis**: Generate statistical comparisons with confidence intervals and effect sizes
4. **Documentation Update**: Refresh all .md files, tables, and figures with new results
5. **Reproducibility**: Generate execution manifests and logs for full auditability

---

## 📂 Current State Assessment

### Existing Files Identified

**Framework Implementations:**
- ✅ `framework_investigativo_completo.py` (PennyLane baseline)
- ✅ `framework_qiskit.py` (Qiskit implementation)
- ✅ `framework_cirq.py` (Cirq implementation)

**Execution Scripts:**
- ✅ `executar_multiframework.py`
- ✅ `executar_multiframework_rapido.py`
- ✅ `executar_framework_qiskit.py`
- ✅ `executar_framework_cirq.py`

**Comparison Tools:**
- ✅ `comparacao_multiframework_completa.py`
- ✅ `generate_comparative_results.py`

**Enhancement Documents:**
- ✅ `MegaPrompt Especializado_ Melhorias no Framework _Beneficial Quantum Noise in VQC_.md`
- ✅ `MegaPrompt Aprimorado v3.0_ Melhorias no Framework _Beneficial Quantum Noise in VQC_.md`
- ✅ `Planjo: Executar Multiframework e Atualizar os Resultados.md`

---

## 🗺️ Implementation Phases

### Phase 1: Repository Inventory & Alignment (2-4 hours)

**Tasks:**
- [x] Identify all framework entrypoints and dependencies
- [ ] Create `docs/inventario_execucao_multiframework.md`
- [ ] Map code structure to execution flow
- [ ] Identify parametrization mechanisms (CLI, YAML, JSON)
- [ ] Document current experimental configurations

**Deliverables:**
- Repository structure map
- Entrypoint documentation
- Current config inventory

### Phase 2: Configuration Standardization (2-3 hours)

**Tasks:**
- [ ] Create unified `configs/experiment.yaml` (or use existing)
- [ ] Define factors and levels consistently
- [ ] Standardize seeds across frameworks
- [ ] Document equivalences and limitations
- [ ] Create `docs/equivalencias_e_limitacoes.md`

**Deliverables:**
- `configs/experiment_unified.yaml`
- Equivalence mapping document
- Seed configuration documentation

### Phase 3: Improvements Implementation (4-8 hours)

**Tasks:**
- [ ] Extract requirements from MegaPrompt documents
- [ ] Create `docs/melhorias_map.md` (improvement tracking)
- [ ] Implement improvements per framework:
  - [ ] PennyLane
  - [ ] Qiskit
  - [ ] Cirq
- [ ] Add smoke tests for each framework
- [ ] Create `tests/smoke_multiframework.py`

**Deliverables:**
- Improvement tracking document
- Updated framework code
- Smoke test suite

### Phase 4: Environment Preparation (1-2 hours)

**Tasks:**
- [ ] Create isolated environments (or validate existing)
- [ ] Freeze dependencies: `requirements-lock.txt`
- [ ] Document Python + library versions
- [ ] Configure global random seeds
- [ ] Test environment reproducibility

**Deliverables:**
- `requirements-lock.txt`
- Environment documentation
- Seed configuration scripts

### Phase 5: Execution - PennyLane Baseline (2-4 hours)

**Tasks:**
- [ ] Run complete pipeline with unified config
- [ ] Generate `logs/pennylane/<run_id>/`
- [ ] Generate `results/pennylane/<run_id>/metrics.csv`
- [ ] Create execution manifest
- [ ] Validate outputs

**Deliverables:**
- PennyLane results (CSV, logs, manifest)
- Baseline metrics for comparison

### Phase 6: Execution - Qiskit (2-4 hours)

**Tasks:**
- [ ] Run with equivalent configuration
- [ ] Generate `logs/qiskit/<run_id>/`
- [ ] Generate `results/qiskit/<run_id>/metrics.csv`
- [ ] Create execution manifest
- [ ] Document any necessary adaptations

**Deliverables:**
- Qiskit results (CSV, logs, manifest)
- Adaptation notes if needed

### Phase 7: Execution - Cirq (2-4 hours)

**Tasks:**
- [ ] Run with equivalent configuration
- [ ] Generate `logs/cirq/<run_id>/`
- [ ] Generate `results/cirq/<run_id>/metrics.csv`
- [ ] Create execution manifest
- [ ] Document any necessary adaptations

**Deliverables:**
- Cirq results (CSV, logs, manifest)
- Adaptation notes if needed

### Phase 8: Comparative Analysis (3-5 hours)

**Tasks:**
- [ ] Standardize CSV schemas across frameworks
- [ ] Run `generate_comparative_results.py`
- [ ] Generate `results/comparisons/<run_id>/comparative_table.csv`
- [ ] Calculate statistical metrics:
  - [ ] Confidence intervals (95%)
  - [ ] Effect sizes
  - [ ] Multiple comparison corrections
- [ ] Create `results/comparisons/<run_id>/stats_report.json`

**Deliverables:**
- Comparative tables
- Statistical analysis report
- Framework performance comparison

### Phase 9: Figure Generation (2-3 hours)

**Tasks:**
- [ ] Generate mandatory figures:
  - [ ] Performance curves by framework
  - [ ] Noise sensitivity analysis
  - [ ] Robustness distributions (by seed)
  - [ ] Cost vs. performance (Pareto)
- [ ] Save to `figures/<run_id>/`
- [ ] Create `figures/<run_id>/figure_manifest.json`

**Deliverables:**
- Figure set (PNG/PDF/SVG)
- Figure generation manifest

### Phase 10: Documentation Update (4-6 hours)

**Tasks:**
- [ ] Update all .md files in `artigo_cientifico/`
- [ ] Embed new tables referencing run_id
- [ ] Update figure links
- [ ] Refresh methodology notes
- [ ] Create `CHANGELOG_EXECUCOES.md`

**Deliverables:**
- Updated article sections
- Changelog documenting changes
- Table/figure updates

### Phase 11: Consistency Verification (2-3 hours)

**Tasks:**
- [ ] Run consistency checks:
  - [ ] Code vs. docs alignment
  - [ ] Config vs. docs alignment
  - [ ] Metrics implementation vs. docs
  - [ ] Reported numbers vs. CSV/JSON
  - [ ] Figures vs. data sources
- [ ] Create `relatorio_consistencia.md`
- [ ] Address any discrepancies

**Deliverables:**
- Consistency report
- Discrepancy resolution log

### Phase 12: Quality Gates (1-2 hours)

**Tasks:**
- [ ] Verify all frameworks executed successfully
- [ ] Confirm comparative results generated
- [ ] Validate markdown references to run_id
- [ ] Check manifests and logs exist
- [ ] Review equivalence documentation

**Deliverables:**
- Quality gate checklist (completed)
- Sign-off for production use

---

## 📊 Deliverables Summary

### Logs & Manifests
```
logs/
├── pennylane/<run_id>/
│   ├── stdout.log
│   ├── stderr.log
│   └── manifest_execucao.json
├── qiskit/<run_id>/
│   ├── stdout.log
│   ├── stderr.log
│   └── manifest_execucao.json
└── cirq/<run_id>/
    ├── stdout.log
    ├── stderr.log
    └── manifest_execucao.json
```

### Results
```
results/
├── pennylane/<run_id>/
│   ├── metrics.csv
│   └── summary.csv
├── qiskit/<run_id>/
│   ├── metrics.csv
│   └── summary.csv
├── cirq/<run_id>/
│   ├── metrics.csv
│   └── summary.csv
└── comparisons/<run_id>/
    ├── comparative_table.csv
    └── stats_report.json
```

### Figures
```
figures/<run_id>/
├── performance_curves.png
├── noise_sensitivity.png
├── robustness_distributions.png
├── cost_vs_performance.png
└── figure_manifest.json
```

### Documentation
```
docs/
├── inventario_execucao_multiframework.md
├── equivalencias_e_limitacoes.md
├── melhorias_map.md
└── inventario_docs.md

CHANGELOG_EXECUCOES.md
relatorio_consistencia.md
```

---

## ⏱️ Estimated Timeline

| Phase | Duration | Dependencies |
|-------|----------|--------------|
| 1. Inventory | 2-4h | None |
| 2. Config Standardization | 2-3h | Phase 1 |
| 3. Improvements | 4-8h | Phase 2 |
| 4. Environment Prep | 1-2h | Phase 3 |
| 5. PennyLane Execution | 2-4h | Phase 4 |
| 6. Qiskit Execution | 2-4h | Phase 4 |
| 7. Cirq Execution | 2-4h | Phase 4 |
| 8. Comparative Analysis | 3-5h | Phases 5-7 |
| 9. Figure Generation | 2-3h | Phase 8 |
| 10. Documentation Update | 4-6h | Phases 8-9 |
| 11. Consistency Verification | 2-3h | Phase 10 |
| 12. Quality Gates | 1-2h | Phase 11 |
| **TOTAL** | **27-48 hours** | Sequential execution |

**Note:** Phases 5-7 can partially overlap if running on different machines/environments.

---

## 🚦 Definition of Done

The implementation is complete when:

- [ ] Golden run_id exists with complete results for all 3 frameworks
- [ ] Comparative statistics generated and coherent
- [ ] All markdown and figures reference the golden run_id
- [ ] Manifests + logs enable full reproduction
- [ ] Zero inconsistencies (or documented with mitigation)
- [ ] All quality gates passed
- [ ] `CHANGELOG_EXECUCOES.md` documents the golden run

---

## 🔧 Next Immediate Steps

### 1. Create Inventory (Start Here)

```bash
# Create docs directory structure
mkdir -p docs/execution_plan
mkdir -p logs/{pennylane,qiskit,cirq}
mkdir -p results/{pennylane,qiskit,cirq,comparisons}
mkdir -p figures

# Start inventory
python -c "
import os
import json

inventory = {
    'frameworks': {
        'pennylane': 'framework_investigativo_completo.py',
        'qiskit': 'framework_qiskit.py',
        'cirq': 'framework_cirq.py'
    },
    'execution_scripts': {
        'multiframework': 'executar_multiframework.py',
        'multiframework_rapid': 'executar_multiframework_rapido.py',
        'qiskit': 'executar_framework_qiskit.py',
        'cirq': 'executar_framework_cirq.py'
    },
    'comparison': {
        'complete': 'comparacao_multiframework_completa.py',
        'generator': 'generate_comparative_results.py'
    }
}

with open('docs/execution_plan/inventory.json', 'w') as f:
    json.dump(inventory, f, indent=2)

print('Inventory created at docs/execution_plan/inventory.json')
"
```

### 2. Review MegaPrompt Documents

```bash
# List all improvement documents
ls -lh "MegaPrompt"* "Planjo"*

# Priority reading order:
# 1. Planjo: Executar Multiframework e Atualizar os Resultados.md
# 2. MegaPrompt Especializado_ Melhorias no Framework _Beneficial Quantum Noise in VQC_.md
# 3. MegaPrompt Aprimorado v3.0_ Melhorias no Framework _Beneficial Quantum Noise in VQC_.md
```

### 3. Identify Existing Configs

```bash
# Find configuration files
find . -name "*.yaml" -o -name "*.json" -o -name "*config*" | grep -v node_modules | grep -v ".git"
```

---

## 📞 Contact & Coordination

For questions or coordination on this execution plan:
- Reference this document: `MULTIFRAMEWORK_EXECUTION_ROADMAP.md`
- Track progress in: `docs/execution_plan/progress.md` (to be created)
- Log issues in: `docs/execution_plan/issues.md` (to be created)

---

## 🎯 Success Metrics

| Metric | Target | Measurement |
|--------|--------|-------------|
| **Framework Parity** | 100% | Identical configs across 3 frameworks |
| **Execution Success Rate** | 100% | All frameworks complete without critical errors |
| **Result Completeness** | 100% | All expected metrics generated |
| **Documentation Currency** | 100% | All .md files reference current run_id |
| **Reproducibility** | 100% | External user can reproduce with manifests |
| **Consistency** | ≥95% | Code↔data↔docs alignment |

---

**Status:** 🚧 READY FOR PHASE 1 - INVENTORY

**Next Action:** Begin Phase 1 - Repository Inventory & Alignment

**Owner:** TBD  
**Reviewer:** @MarceloClaro  
**Estimated Start:** Ready to begin  
**Estimated Completion:** 27-48 hours of focused work
