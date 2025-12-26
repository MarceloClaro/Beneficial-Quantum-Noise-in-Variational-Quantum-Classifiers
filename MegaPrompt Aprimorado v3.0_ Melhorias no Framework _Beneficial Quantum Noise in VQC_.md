# MegaPrompt Aprimorado v3.0: Melhorias no Framework "Beneficial Quantum Noise in VQC"

## 🎯 OBJETIVO GERAL

Refatorar e expandir o framework "Beneficial Quantum Noise in VQC" para alcançar o mais alto padrão de rigor matemático, reprodutibilidade e auditabilidade, garantindo conformidade com os critérios de periódicos Qualis A1 como Nature, Quantum e Physical Review. Este prompt integra as melhores práticas do MegaPrompt v2.0 com o bundle multiframework, criando um sistema de execução e auditoria de classe mundial.

---

## 📖 ÍNDICE

### PARTE I: Princípios e Configuração
1. [Princípios Não Negociáveis](#principios)
2. [Seção 0: Configuração do Projeto de Melhoria](#secao-0)
3. [Glossário de Melhorias](#glossario)

### PARTE II: Automação e Quality Gates
4. [Automação com Task](#automacao)
5. [Quality Gates (E1-E4 + FINAL)](#quality-gates)

### PARTE III: Execução das Melhorias (12 tarefas)
6. [Tarefa 1: Documentação Matemática Formal](#tarefa-1)
7. [Tarefa 2: Validação Matemática dos Operadores de Kraus](#tarefa-2)
8. [Tarefa 3: Derivação Formal do QNG](#tarefa-3)
9. [Tarefa 4: Centralização e Documentação de Seeds](#tarefa-4)
10. [Tarefa 5: Geração de Manifesto de Execução](#tarefa-5)
11. [Tarefa 6: Correção de Bonferroni nos Testes Post-Hoc](#tarefa-6)
12. [Tarefa 7: Análise de Poder Estatístico](#tarefa-7)
13. [Tarefa 8: Geração de Tabela Código→Método](#tarefa-8)
14. [Tarefa 9: Integração com Cirq e Qiskit](#tarefa-9)
15. [Tarefa 10: Geração de Diagramas de Circuitos](#tarefa-10)
16. [Tarefa 11: Validação com Schemas JSON](#tarefa-11)
17. [Tarefa 12: Documento de Equivalências e Limitações](#tarefa-12)

### PARTE IV: Validação e Entrega
18. [Checklist de Conformidade Qualis A1](#checklist)
19. [Entrega Final (Pull Request)](#entrega)

---

<a name="principios"></a>
## PARTE I: Princípios e Configuração

### Princípios Não Negociáveis

1. **Sem Números Inventados:** Todo número deve existir em CSV/JSON gerado pelo pipeline; caso contrário marque **[NÃO DISPONÍVEL]**.
2. **Sem Inferências de Configuração:** Fatores/níveis/seeds devem vir de `configs/experiment.yaml`; caso contrário **[INFORMAÇÃO AUSENTE]**.
3. **Paridade entre Frameworks:** Diferenças inevitáveis devem ser registradas em `docs/equivalencias_e_limitacoes.md`.
4. **Audit Trail Completo:** Toda seção/figura/tabela deve apontar fonte (arquivo/config/log) e `run_id`.

---

<a name="secao-0"></a>
### Seção 0: Configuração do Projeto de Melhoria

**Instrução:** Clone o repositório e crie a estrutura de diretórios.

```bash
# 1. Clone o repositório
git clone https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers.git
cd Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers

# 2. Crie um novo branch para as melhorias
git checkout -b feature/qualis-a1-v3

# 3. Crie a estrutura de diretórios
mkdir -p configs schemas scripts logs results figures manifests docs
```

**Arquivo de Configuração `configs/experiment.yaml`:**
Crie o arquivo `configs/experiment.yaml` como fonte única de verdade.

```yaml
# configs/experiment.yaml
run:
  run_id: "AUTO"              # se AUTO, gerar timestamp + hash curto
  output_root: "results"

reproducibility:
  seeds: [1, 2, 3, 4, 5]

model:
  n_qubits: 4
  depth_L: 3
  ansatz: "StronglyEntangling"

noise:
  enabled: true
  model: "depolarizing"
  params:
    p: [0.0, 0.01, 0.02]

frameworks:
  pennylane:
    backend: "default.qubit"
  qiskit:
    backend: "aer_simulator"
  cirq:
    backend: "cirq.Simulator"
```

---

<a name="glossario"></a>
### Glossário de Melhorias

- **run_id:** Identificador único que rastreia todos os artefatos de uma execução.
- **Manifesto de Execução:** Arquivo JSON que registra todas as configurações, versões de bibliotecas e comandos de uma execução.
- **Paridade Experimental:** Garantia de que os experimentos são executados com as mesmas configurações em diferentes frameworks.
- **Quality Gates:** Checkpoints obrigatórios que garantem a qualidade e a consistência do pipeline.

---

<a name="automacao"></a>
## PARTE II: Automação e Quality Gates

### Automação com Task

**Instrução:** Crie um `Taskfile.yml` na raiz do projeto para automatizar todo o pipeline.

```yaml
# Taskfile.yml
version: '3'

tasks:
  default:
    cmds:
      - task: all

  setup:
    desc: "Instalar dependências"
    cmds:
      - pip install -r requirements.txt

  smoke:
    desc: "Smoke test (1 config, 1 seed)"
    cmds:
      - python scripts/run_pennylane.py --config configs/experiment.yaml --smoke
      - python scripts/run_qiskit.py --config configs/experiment.yaml --smoke
      - python scripts/run_cirq.py --config configs/experiment.yaml --smoke

  all:
    desc: "Pipeline completo"
    cmds:
      - task: smoke
      - task: run:pl
      - task: run:qiskit
      - task: run:cirq
      - task: compare
      - task: update-docs
      - task: audit
```

---

<a name="quality-gates"></a>
### Quality Gates (E1-E4 + FINAL)

**Instrução:** Implemente um script `scripts/audit_consistency.py` que verifique os Quality Gates.

- **Gate E1:** Smoke test passa nos 3 frameworks.
- **Gate E2:** Execução completa gera `metrics.csv` válido.
- **Gate E3:** Comparativo gera `comparative_table.csv` e `stats_report.json`.
- **Gate E4:** Docs atualizados referenciam `run_id` e valores existem em artefatos.
- **Gate FINAL:** Auditoria de consistência sem discrepâncias não explicadas.

---

<a name="tarefa-1"></a>
## PARTE III: Execução das Melhorias (12 tarefas)

### Tarefa 1: Documentação Matemática Formal
**Objetivo:** Adicionar docstrings com equações LaTeX a todas as 11 classes de ruído.

**Instrução:** Para cada classe de ruído, adicione um docstring detalhado com a descrição matemática, os operadores de Kraus e as referências.

---

### Tarefa 2: Validação Matemática dos Operadores de Kraus
**Objetivo:** Adicionar um método de validação para os operadores de Kraus.

**Instrução:** Crie um módulo `scripts/validation.py` com a função `validar_operadores_kraus` e integre-o ao pipeline.

---

### Tarefa 3: Derivação Formal do QNG
**Objetivo:** Adicionar documentação matemática à classe `QNG`.

**Instrução:** Adicione um docstring detalhado à classe `QNG` com a derivação do Quantum Natural Gradient e as referências.

---

### Tarefa 4: Centralização e Documentação de Seeds
**Objetivo:** Criar uma função centralizada para configuração de seeds.

**Instrução:** Crie um módulo `scripts/reproducibility.py` com a função `configurar_seeds_reprodutiveis` e chame-a no início de cada script de execução.

---

### Tarefa 5: Geração de Manifesto de Execução
**Objetivo:** Criar um arquivo JSON que documente cada execução.

**Instrução:** No módulo `scripts/reproducibility.py`, adicione a função `gerar_manifesto_execucao` com os campos do bundle (commit, host, dependências, etc.).

---

### Tarefa 6: Correção de Bonferroni nos Testes Post-Hoc
**Objetivo:** Adicionar correção para múltiplas comparações.

**Instrução:** Modifique a classe `TestesEstatisticosAvancados` para incluir um método `testes_post_hoc_com_correcao` com `method='bonferroni'`.

---

### Tarefa 7: Análise de Poder Estatístico
**Objetivo:** Adicionar cálculo de poder estatístico (1-β).

**Instrução:** Na classe `TestesEstatisticosAvancados`, adicione um método `analise_poder_estatistico`.

---

### Tarefa 8: Geração de Tabela Código→Método
**Objetivo:** Criar um mapeamento explícito entre o artigo e o código.

**Instrução:** Crie um módulo `scripts/auditing.py` com a função `gerar_tabela_codigo_metodo`.

---

### Tarefa 9: Integração com Cirq e Qiskit
**Objetivo:** Aumentar a generalidade do framework.

**Instrução:** Crie os scripts `scripts/run_cirq.py` e `scripts/run_qiskit.py` que leem `configs/experiment.yaml` e produzem `metrics.csv` no mesmo schema.

---

### Tarefa 10: Geração de Diagramas de Circuitos
**Objetivo:** Melhorar a didática do artigo com visualizações.

**Instrução:** Crie um módulo `scripts/visualization.py` com a função `gerar_diagrama_circuito`.

---

### Tarefa 11: Validação com Schemas JSON
**Objetivo:** Garantir a consistência dos dados gerados.

**Instrução:** Crie um script `scripts/validate_schemas.py` que valide `metrics.csv` e `stats_report.json` contra os schemas em `schemas/`.

---

### Tarefa 12: Documento de Equivalências e Limitações
**Objetivo:** Documentar diferenças inevitáveis entre frameworks.

**Instrução:** Crie o arquivo `docs/equivalencias_e_limitacoes.md` e documente as diferenças na implementação de ruído, otimizadores, etc.

---

<a name="checklist"></a>
## PARTE IV: Validação e Entrega

### Checklist de Conformidade Qualis A1

**1. Rigor Matemático (30 pts)**
- [ ] Docstrings com equações LaTeX (10 pts)
- [ ] Validação de operadores de Kraus (10 pts)
- [ ] Derivação do QNG (10 pts)

**2. Reprodutibilidade (30 pts)**
- [ ] Seeds centralizadas (10 pts)
- [ ] Manifesto de execução completo (10 pts)
- [ ] Configuração única (`experiment.yaml`) (10 pts)

**3. Rigor Estatístico (20 pts)**
- [ ] Correção de Bonferroni (10 pts)
- [ ] Análise de poder (10 pts)

**4. Auditoria e Transparência (20 pts)**
- [ ] Tabela Código→Método (5 pts)
- [ ] Integração Cirq/Qiskit (5 pts)
- [ ] Diagramas de circuitos (5 pts)
- [ ] Validação com schemas (5 pts)

**Pontuação Final:** [Soma dos pontos] / 100

---

<a name="entrega"></a>
### Entrega Final (Pull Request)

1. ✅ Crie um Pull Request do branch `feature/qualis-a1-v3` para o `main`.
2. ✅ No corpo do PR, inclua:
   - Resumo das 12 melhorias implementadas.
   - Pontuação final do Checklist de Conformidade Qualis A1.
   - Link para o `run_id` da execução de validação.
3. ✅ Solicite revisão de pelo menos 2 coautores.
4. ✅ Após aprovação, faça o merge para o `main`.

---

**FIM DO MEGAPROMPT**
GAPROMPT**
