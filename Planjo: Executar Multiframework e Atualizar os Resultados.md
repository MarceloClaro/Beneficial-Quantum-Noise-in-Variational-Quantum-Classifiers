# Plano Operacional (MegaPrompt) — Execução Multiframework e Atualização Integral de Resultados
**Escopo:** Executar e comparar implementações em **PennyLane, Qiskit e Cirq**, aplicando as melhorias descritas em `MegaPrompt Especializado_ Melhorias no Framework Beneficial Quantum Noise in VQC.md`, e atualizar **todos** os artefatos (dados, tabelas, figuras e Markdown) com rastreabilidade e reprodutibilidade.

---

## 1) Compreensão da tarefa (declaração operacional)
Você deve:
1) **Aplicar melhorias** (refactor, correções, novas métricas/rotinas, ajustes estatísticos, etc.) aos três frameworks, garantindo **paridade funcional** (mesma hipótese experimental, mesmos fatores, mesmos níveis, mesmas métricas).
2) **Executar o pipeline completo** em cada framework (PennyLane como baseline, Qiskit e Cirq como replicações).
3) **Gerar resultados comparativos** entre frameworks (tabelas, estatísticas, figuras).
4) **Atualizar toda a documentação** (arquivos `.md`, tabelas e figuras) para refletir os novos resultados.
5) **Gerar manifestos de execução** e logs para reprodutibilidade/auditoria.
6) **Validar consistência**: código ↔ dados ↔ texto ↔ figuras.

**Regras de integridade (não negociáveis):**
- Não inventar números. Se algo não foi gerado pelo pipeline, marcar como **[NÃO DISPONÍVEL]**.
- Não inferir configuração experimental ausente. Se faltar, marcar **[INFORMAÇÃO AUSENTE]** e registrar pendência.
- Toda atualização deve ser rastreável: *qual commit/arquivo/execução gerou qual resultado*.

---

## 2) Entradas obrigatórias
1) Repositório do projeto (código completo).
2) Documento de melhorias: `MegaPrompt Especializado_ Melhorias no Framework Beneficial Quantum Noise in VQC.md`.
3) Dados (ou instruções oficiais de obtenção) e licenças.
4) Scripts de execução existentes (se houver) ou entrypoints definidos.

---

## 3) Saídas obrigatórias (artefatos a entregar)
### 3.1 Logs e manifestos (reprodutibilidade)
- `logs/<framework>/<run_id>/stdout.log`
- `logs/<framework>/<run_id>/stderr.log`
- `manifests/<framework>/<run_id>/manifest_execucao.json` contendo:
  - hash/commit do código
  - ambiente (SO, CPU/GPU, RAM)
  - versões exatas (Python + libs)
  - seeds
  - parâmetros e configs
  - comandos executados
  - tempo de execução e status

### 3.2 Resultados tabulares e métricas (dados “fonte da verdade”)
- `results/<framework>/<run_id>/metrics.csv` (por configuração/seed)
- `results/<framework>/<run_id>/summary.csv` (agregado por condição)
- `results/comparisons/<run_id>/comparative_table.csv`
- `results/comparisons/<run_id>/stats_report.json` (testes, IC, effect sizes)

### 3.3 Figuras
- `figures/<run_id>/...` (PNG/PDF/SVG conforme padrão do repositório)
- `figures/<run_id>/figure_manifest.json` (fonte, script e parâmetros)

### 3.4 Documentação atualizada
- Atualização de todos os `.md` de resultados (mantendo coerência com tabelas/figuras).
- `CHANGELOG_EXECUCOES.md` (o que mudou, por quê, e qual run_id valida).

---

## 4) Inventário e alinhamento do repositório (antes de executar)
### 4.1 Analisar estrutura e entrypoints
- Identificar:
  - diretórios de código, scripts, configs, resultados e docs
  - como o pipeline é parametrizado (CLI, YAML/JSON, constantes)
  - onde métricas são calculadas e persistidas
- Listar “arquivos-chave” com função e dependências.

### 4.2 Padronizar configurações experimentais (paridade entre frameworks)
Definir, em um único arquivo de configuração (se ainda não existir), por exemplo:
- `configs/experiment.yaml` (ou `.json`) com:
  - fatores e níveis
  - número de seeds e lista de seeds
  - número de repetições
  - parâmetros do ansatz (n_qubits, depth L, etc.)
  - noise model e parâmetros
  - otimizador e hiperparâmetros
  - datasets e splits
  - métricas a coletar
  - critérios de convergência

**Regra:** Qiskit e Cirq devem ler exatamente a mesma config e produzir métricas no mesmo schema.

---

## 5) Aplicação das melhorias (controle de mudanças)
### 5.1 Extrair requisitos do arquivo de melhorias
- Converter cada melhoria em um item rastreável:
  - `ID`, descrição, arquivos afetados, critérios de aceite, impacto esperado.
- Criar `docs/melhorias_map.md` com:
  - lista de melhorias
  - status (implementado/não implementado)
  - evidência (commit/arquivo/teste)

### 5.2 Implementar melhorias com testes mínimos
- Para cada framework:
  - implementar a mudança
  - adicionar/verificar testes (unitários ou “smoke tests” de execução curta)
- Criar `tests/smoke_multiframework.py` (ou equivalente) para validar:
  - importação do framework
  - execução mínima de 1 configuração
  - escrita de outputs no schema esperado

---

## 6) Preparação do ambiente (determinística)
### 6.1 Gerenciamento de dependências
- Usar ambiente isolado por framework ou um ambiente único com constraints.
- Congelar versões:
  - `requirements.txt` + `requirements-lock.txt` (ou `conda env export`)
- Validar:
  - `python --version`
  - `pip freeze` (persistir no manifesto)

### 6.2 Controle de aleatoriedade
- Definir seeds globais e por biblioteca:
  - Python `random`, NumPy, e seeds dos simuladores/backends.
- Persistir seeds no `manifest_execucao.json`.

---

## 7) Execução do pipeline por framework (ordem e critérios)
### 7.1 Ordem mandatória
1) PennyLane (baseline)
2) Qiskit
3) Cirq
4) Comparativo (script agregador)

### 7.2 Execução PennyLane (baseline)
- Rodar o pipeline completo conforme config.
- Gerar:
  - métricas por seed/configuração
  - agregados
  - logs e manifesto

### 7.3 Execução Qiskit
- Rodar em simulador definido (ex.: Aer) e backend consistente.
- Garantir equivalência de:
  - ansatz (ou tradução formal equivalente)
  - noise model (ou aproximação documentada)
  - métrica e função objetivo

### 7.4 Execução Cirq
- Rodar em simulador definido e pipeline equivalente.
- Mesmas garantias de equivalência do item 7.3.

**Regra de equivalência:** Se houver divergência inevitável (ex.: modelagem de ruído), registrar em:
- `docs/equivalencias_e_limitacoes.md` com:
  - diferença
  - justificativa
  - impacto esperado
  - como comparar de forma justa

---

## 8) Comparação entre frameworks (comparative science)
### 8.1 Padronização de schema
Antes de comparar, validar que `metrics.csv` de cada framework contém colunas obrigatórias, por exemplo:
- `framework`, `dataset`, `seed`, `config_id`
- fatores (colunas por fator)
- métricas (accuracy, loss, etc.)
- custo (tempo, portas, profundidade, memória se disponível)

### 8.2 Estatística comparativa
- Produzir:
  - tabelas por condição
  - IC 95% (por condição e por framework)
  - tamanhos de efeito (quando aplicável)
  - correção para múltiplas comparações (se aplicável)
- Salvar tudo em `results/comparisons/<run_id>/`.

### 8.3 Figuras obrigatórias (mínimo)
- Curvas de desempenho por framework (mesma escala).
- Sensibilidade ao ruído / parâmetro-chave (por framework).
- Robustez por seeds (distribuições).
- Custo vs desempenho (fronteira de Pareto, se houver dados).

---

## 9) Atualização de Markdown, tabelas e figuras (controle editorial)
### 9.1 Fonte da verdade
- Markdown nunca “inventa” valores: ele deve puxar de:
  - `summary.csv`
  - `comparative_table.csv`
  - `stats_report.json`
  - figuras versionadas por `run_id`

### 9.2 Mecanismo de atualização
- Centralizar referências a resultados por `run_id`.
- Atualizar:
  - seções de resultados
  - tabelas embutidas
  - links para figuras
  - notas metodológicas afetadas por melhorias

### 9.3 Registro de alterações
- Atualizar `CHANGELOG_EXECUCOES.md` com:
  - o que mudou (melhoria aplicada)
  - quais resultados mudaram (e por quê)
  - qual run_id é o “golden run”

---

## 10) Verificações finais (auditoria de consistência)
### 10.1 Checklist de consistência código ↔ resultados ↔ documentação
- Arquiteturas/modelos: código = docs
- Fatores e níveis: config = docs
- Total de configurações: enumerador = docs
- Métricas: implementação = docs
- Números reportados no texto: existem em CSV/JSON
- Figuras: geradas a partir dos dados do run_id

### 10.2 Quality Gate final (não finalizar sem cumprir)
- [ ] Todos os frameworks executaram sem erro (ou erros documentados com causa e workaround).
- [ ] Resultados comparativos gerados e versionados.
- [ ] Markdown atualizado aponta para run_id válido.
- [ ] Manifestos e logs existem para cada framework.
- [ ] Diferenças de equivalência (se houver) documentadas.

---

## 11) Arquivos-chave (a localizar e confirmar no repositório)
**Linha base (PennyLane):**
- `framework_investigativo_completo.py`

**Qiskit:**
- `framework_qiskit.py`
- `executar_framework_qiskit.py`

**Cirq:**
- `framework_cirq.py`
- `executar_framework_cirq.py`

**Comparação:**
- `generate_comparative_results.py`

**Documentação:**
- múltiplos `.md` de resultados (inventariar e listar em `docs/inventario_docs.md`)

---

## 12) Entregáveis de documentação (obrigatórios)
1) `docs/inventario_execucao_multiframework.md` (o plano + status)
2) `docs/equivalencias_e_limitacoes.md` (paridade e desvios)
3) `docs/melhorias_map.md` (melhoria → implementação → evidência)
4) `CHANGELOG_EXECUCOES.md` (run_id “golden” e mudanças)
5) `relatorio_consistencia.md` (discrepâncias e correções)

---

## 13) Definição de “pronto” (Definition of Done)
O trabalho está pronto quando:
- Existe um run_id final com resultados completos para os 3 frameworks;
- Comparativos gerados e coerentes;
- Todos os Markdown e figuras atualizados referenciam esse run_id;
- Manifestos + logs permitem reprodução;
- Inconsistências são zero, ou documentadas com mitigação clara.
