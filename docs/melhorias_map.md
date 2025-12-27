# Mapeamento de Melhorias - MegaPrompt para Implementação

**Data:** 2025-12-27  
**Versão:** 1.0  
**Fonte:** `MegaPrompt Especializado_ Melhorias no Framework _Beneficial Quantum Noise in VQC_.md`

---

## 📋 Visão Geral

Este documento mapeia cada melhoria proposta no MegaPrompt para sua implementação concreta, status atual, e evidências de conclusão.

---

## 🎯 10 Tarefas Principais

### Tarefa 1: Documentação Matemática Formal

**ID:** IMP-001  
**Prioridade:** Alta  
**Categoria:** Documentação

**Descrição:**
Adicionar docstrings com equações LaTeX a todas as 11 classes de ruído em `framework_investigativo_completo.py`.

**Arquivos Afetados:**
- `framework_investigativo_completo.py` (linhas ~150-750)

**Critérios de Aceite:**
- [ ] Cada classe de ruído tem docstring com equação matemática em LaTeX
- [ ] Operadores de Kraus documentados
- [ ] Referências a Nielsen & Chuang incluídas
- [ ] Exemplos de uso fornecidos

**Status:** 🔴 Não Implementado

**Evidência:** N/A

**Impacto Esperado:**
- Melhora reprodutibilidade matemática
- Facilita auditoria por revisores
- Aumenta pontuação QUALIS A1 (rigor matemático)

**Notas de Implementação:**
```python
# Exemplo para RuidoDepolarizante
class RuidoDepolarizante(ModeloRuido):
    """
    Modelo de ruído despolarizante (depolarizing noise).
    
    Descrição Matemática:
    $$\\mathcal{E}(\\rho) = (1-p)\\rho + \\frac{p}{3}(X\\rho X + Y\\rho Y + Z\\rho Z)$$
    
    Operadores de Kraus:
    $$K_0 = \\sqrt{1-p}I, K_1 = \\sqrt{p/3}X, K_2 = \\sqrt{p/3}Y, K_3 = \\sqrt{p/3}Z$$
    
    Referência: Nielsen & Chuang, "Quantum Computation and Quantum Information", Cap. 8
    """
```

---

### Tarefa 2: Validação Matemática dos Operadores de Kraus

**ID:** IMP-002  
**Prioridade:** Alta  
**Categoria:** Testes

**Descrição:**
Criar módulo `qualis_a1_modules/validador_kraus.py` com testes que verificam completude (∑ K†K = I) e trace-preserving para todos os modelos de ruído.

**Arquivos Afetados:**
- `qualis_a1_modules/validador_kraus.py` (novo)
- `tests/test_kraus_validation.py` (novo)

**Critérios de Aceite:**
- [ ] Função `validar_completude_kraus()` implementada
- [ ] Função `validar_trace_preserving()` implementada
- [ ] Testes unitários passando para 11 modelos de ruído
- [ ] Tolerância numérica: 1e-10

**Status:** 🔴 Não Implementado

**Evidência:** N/A

**Impacto Esperado:**
- Garante correção matemática dos modelos
- Aumenta confiança dos revisores
- Detecta bugs em implementações futuras

---

### Tarefa 3: Derivação Formal do QNG (Quantum Natural Gradient)

**ID:** IMP-003  
**Prioridade:** Média  
**Categoria:** Documentação + Código

**Descrição:**
Documentar derivação formal do QNG e implementar cálculo da métrica de Fisher quântica em `qualis_a1_modules/qng_formal.py`.

**Arquivos Afetados:**
- `qualis_a1_modules/qng_formal.py` (novo)
- `docs/qng_derivation.md` (novo)

**Critérios de Aceite:**
- [ ] Derivação matemática completa em `docs/qng_derivation.md`
- [ ] Implementação de `calcular_metric_fisher_quantica()`
- [ ] Comparação QNG vs. Adam em benchmark
- [ ] Documentação de trade-offs (velocidade vs. precisão)

**Status:** 🔴 Não Implementado

**Evidência:** N/A

**Impacto Esperado:**
- Justifica escolha de otimizador
- Adiciona profundidade metodológica
- Possível melhoria de performance

---

### Tarefa 4: Centralização e Documentação de Seeds

**ID:** IMP-004  
**Prioridade:** **Crítica**  
**Categoria:** Reprodutibilidade

**Descrição:**
Centralizar todos os seeds em `qai_config.json` (agora `configs/experiment_unified.yaml`) e documentar controle de aleatoriedade.

**Arquivos Afetados:**
- `configs/experiment_unified.yaml` (criado ✅)
- `framework_investigativo_completo.py` (modificar para ler config)
- `framework_qiskit.py` (modificar para ler config)
- `framework_cirq.py` (modificar para ler config)
- `docs/reproducibility_protocol.md` (novo)

**Critérios de Aceite:**
- [x] Configuração unificada criada
- [ ] Seeds centralizados (global, numpy, random, torch, frameworks)
- [ ] Scripts adaptados para ler do config
- [ ] Documentação de protocolo de reprodutibilidade
- [ ] Validação: mesma seed → mesmos resultados (determinismo)

**Status:** 🟡 Parcialmente Implementado (config criado, scripts não adaptados)

**Evidência:** `configs/experiment_unified.yaml` criado

**Impacto Esperado:**
- **Crítico** para reprodutibilidade
- Facilita auditoria
- Atende requisito mandatório QUALIS A1

**Próximos Passos:**
1. Adaptar scripts para ler `seeds` do YAML
2. Implementar função `set_all_seeds(config)` em módulo comum
3. Validar determinismo com testes

---

### Tarefa 5: Geração de Manifesto de Execução

**ID:** IMP-005  
**Prioridade:** **Crítica**  
**Categoria:** Reprodutibilidade

**Descrição:**
Implementar `qualis_a1_modules/manifesto_gerador.py` que registra tudo sobre uma execução (commit hash, versões de libs, comando executado, seeds, tempo, status).

**Arquivos Afetados:**
- `qualis_a1_modules/manifesto_gerador.py` (novo)
- `manifests/<framework>/<run_id>/manifest_execucao.json` (output)

**Critérios de Aceite:**
- [ ] Função `gerar_manifesto()` implementada
- [ ] Captura: commit hash, SO, Python version, pip freeze
- [ ] Captura: comando executado, seeds, config
- [ ] Captura: timestamps (início, fim), tempo total
- [ ] Captura: status de saída (sucesso/erro)
- [ ] JSON válido e versionado

**Status:** 🔴 Não Implementado

**Evidência:** N/A

**Impacto Esperado:**
- **Crítico** para reprodutibilidade completa
- Permite auditoria externa
- Facilita debugging de experimentos

**Exemplo de Manifesto:**
```json
{
  "manifest_version": "1.0",
  "run_id": "20251227_001_pennylane",
  "timestamp_start": "2025-12-27T10:00:00Z",
  "timestamp_end": "2025-12-27T12:30:00Z",
  "duration_seconds": 9000,
  "git_commit": "d03683e",
  "git_branch": "copilot/generate-scientific-article-yet-again",
  "environment": {
    "os": "Ubuntu 22.04",
    "python_version": "3.10.12",
    "cpu": "Intel Xeon 8-core",
    "ram_gb": 32
  },
  "dependencies": {
    "pennylane": "0.33.1",
    "numpy": "1.24.3",
    "scipy": "1.11.4"
  },
  "command": "python framework_investigativo_completo.py --config configs/experiment_unified.yaml",
  "config_file": "configs/experiment_unified.yaml",
  "seeds": {
    "global": 42,
    "numpy": 42,
    "random": 42
  },
  "status": "success",
  "outputs": {
    "metrics_csv": "results/pennylane/20251227_001/metrics.csv",
    "summary_csv": "results/pennylane/20251227_001/summary.csv"
  }
}
```

---

### Tarefa 6: Correção de Bonferroni nos Testes Post-Hoc

**ID:** IMP-006  
**Prioridade:** Alta  
**Categoria:** Análise Estatística

**Descrição:**
Implementar correção de Bonferroni em `qualis_a1_modules/stats_aprimorada.py` para controlar erro tipo I em comparações múltiplas.

**Arquivos Afetados:**
- `qualis_a1_modules/stats_aprimorada.py` (novo)
- Scripts de análise estatística (modificar)

**Critérios de Aceite:**
- [ ] Função `corrigir_bonferroni(p_values, alpha=0.05)` implementada
- [ ] Integração com testes post-hoc (Tukey, Scheffé)
- [ ] Relatório de p-valores ajustados vs. originais
- [ ] Documentação de quando usar vs. FDR (Benjamini-Hochberg)

**Status:** 🔴 Não Implementado

**Evidência:** N/A

**Impacto Esperado:**
- Corrige análise estatística (atualmente sem correção)
- Atende requisito QUALIS A1 de rigor estatístico
- Pode mudar significância de alguns resultados

---

### Tarefa 7: Análise de Poder Estatístico

**ID:** IMP-007  
**Prioridade:** Média  
**Categoria:** Análise Estatística

**Descrição:**
Implementar análise de poder estatístico (power analysis) para validar que o tamanho amostral é suficiente para detectar efeitos.

**Arquivos Afetados:**
- `qualis_a1_modules/power_analysis.py` (novo)
- Relatórios de resultados (adicionar seção de poder)

**Critérios de Aceite:**
- [ ] Função `calcular_poder(effect_size, n, alpha=0.05)` implementada
- [ ] Cálculo para cada comparação principal
- [ ] Target: poder ≥ 0.80
- [ ] Recomendação de n se poder < 0.80

**Status:** 🔴 Não Implementado

**Evidência:** N/A

**Impacto Esperado:**
- Valida robustez estatística
- Responde preocupação comum de revisores
- Pode justificar tamanho amostral atual

---

### Tarefa 8: Geração de Tabela Código→Método

**ID:** IMP-008  
**Prioridade:** Média  
**Categoria:** Auditoria

**Descrição:**
Gerar tabela que mapeia seções do artigo para linhas de código, facilitando auditoria.

**Arquivos Afetados:**
- `qualis_a1_modules/auditoria_mapeador.py` (novo)
- `docs/codigo_metodo_mapping.md` (output)

**Critérios de Aceite:**
- [ ] Tabela com colunas: Seção Artigo | Método/Técnica | Arquivo | Linhas | Comentário
- [ ] Cobertura: 100% dos métodos mencionados no artigo
- [ ] Formato: Markdown + CSV

**Status:** 🔴 Não Implementado

**Evidência:** N/A

**Impacto Esperado:**
- Facilita revisão por pares
- Aumenta transparência
- Acelera processo de publicação

**Exemplo:**
| Seção | Método | Arquivo | Linhas | Comentário |
|-------|--------|---------|--------|------------|
| 3.2.6 | TREX Error Mitigation | trex_error_mitigation.py | 45-128 | Inversão de matriz confusão |
| 3.2.7 | AUEC Framework | adaptive_unified_error_correction.py | 67-245 | Filtro de Kalman adaptativo |

---

### Tarefa 9: Integração com Cirq e Qiskit

**ID:** IMP-009  
**Prioridade:** **Crítica** (parte do plano multiframework)  
**Categoria:** Multiframework

**Descrição:**
Garantir paridade funcional entre PennyLane, Qiskit e Cirq, documentando equivalências e limitações.

**Arquivos Afetados:**
- `framework_qiskit.py` (validar/adaptar)
- `framework_cirq.py` (validar/adaptar)
- `docs/equivalencias_e_limitacoes.md` (novo)

**Critérios de Aceite:**
- [ ] Mesma configuração experimental nos 3 frameworks
- [ ] Equivalências documentadas (ansätze, noise models)
- [ ] Limitações documentadas (diferenças inevitáveis)
- [ ] Testes de equivalência (mesmo seed → resultados similares)

**Status:** 🟡 Parcialmente Implementado (frameworks existem, paridade não validada)

**Evidência:** `framework_qiskit.py` e `framework_cirq.py` existem

**Impacto Esperado:**
- Valida portabilidade do fenômeno
- Primeira triangulação multiframework na literatura
- Fortalece argumentação científica

**Próximos Passos:**
- Criar `docs/equivalencias_e_limitacoes.md`
- Validar que mesma config gera outputs comparáveis
- Documentar diferenças (e.g., transpilation overhead)

---

### Tarefa 10: Geração de Diagramas de Circuitos

**ID:** IMP-010  
**Prioridade:** Baixa  
**Categoria:** Visualização

**Descrição:**
Gerar diagramas visuais dos circuitos quânticos para cada ansatz usando `circuit.draw()` de cada framework.

**Arquivos Afetados:**
- `qualis_a1_modules/circuit_visualizer.py` (novo)
- `figures/circuits/` (outputs)

**Critérios de Aceite:**
- [ ] Diagramas para 7 ansätze principais
- [ ] Formato: PNG + SVG (alta resolução)
- [ ] Anotações de: qubits, portas, profundidade
- [ ] Inclusão no material suplementar do artigo

**Status:** 🔴 Não Implementado

**Evidência:** N/A

**Impacto Esperado:**
- Melhora clareza metodológica
- Facilita compreensão para não-especialistas
- Valorizado por revisores e leitores

---

## 📊 Resumo de Status

| ID | Tarefa | Prioridade | Status | Evidência |
|----|--------|------------|--------|-----------|
| IMP-001 | Documentação Matemática | Alta | 🔴 Não Impl. | - |
| IMP-002 | Validação Kraus | Alta | 🔴 Não Impl. | - |
| IMP-003 | QNG Derivação | Média | 🔴 Não Impl. | - |
| **IMP-004** | **Seeds Centralização** | **Crítica** | 🟡 Parcial | config criado |
| **IMP-005** | **Manifesto Execução** | **Crítica** | 🔴 Não Impl. | - |
| IMP-006 | Correção Bonferroni | Alta | 🔴 Não Impl. | - |
| IMP-007 | Análise Poder | Média | 🔴 Não Impl. | - |
| IMP-008 | Tabela Código→Método | Média | 🔴 Não Impl. | - |
| **IMP-009** | **Multiframework Paridade** | **Crítica** | 🟡 Parcial | frameworks existem |
| IMP-010 | Diagramas Circuitos | Baixa | 🔴 Não Impl. | - |

**Legenda:**
- 🔴 Não Implementado
- 🟡 Parcialmente Implementado
- 🟢 Implementado
- ✅ Verificado

---

## 🎯 Priorização para Implementação Imediata

### Fase Atual (Infraestrutura):

1. **IMP-004 (Crítica):** Completar centralização de seeds
   - Adaptar scripts para ler config
   - Implementar `set_all_seeds()`
   - Validar determinismo

2. **IMP-005 (Crítica):** Implementar gerador de manifestos
   - Criar `manifesto_gerador.py`
   - Integrar em todos os scripts de execução

3. **IMP-009 (Crítica):** Documentar equivalências multiframework
   - Criar `docs/equivalencias_e_limitacoes.md`
   - Validar paridade de configurações

### Fase Posterior (Análise):

4. **IMP-006 (Alta):** Implementar correção Bonferroni
5. **IMP-002 (Alta):** Validação matemática Kraus
6. **IMP-001 (Alta):** Documentação matemática
7. **IMP-007 (Média):** Análise de poder
8. **IMP-008 (Média):** Tabela código→método
9. **IMP-003 (Média):** QNG derivação
10. **IMP-010 (Baixa):** Diagramas circuitos

---

## 📝 Tracking de Progresso

**Atualizado:** 2025-12-27  
**Responsável:** Equipe de Desenvolvimento  
**Reviewer:** @MarceloClaro

**Próxima Revisão:** Após conclusão de IMP-004, IMP-005, IMP-009
