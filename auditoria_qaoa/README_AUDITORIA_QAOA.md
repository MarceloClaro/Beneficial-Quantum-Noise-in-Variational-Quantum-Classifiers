# Auditoria QAOA — QUALIS A1

Este diretório contém artefatos de auditoria e consolidação de resultados para experimentos QAOA com análise de ruído benéfico.

## Artefatos

- `auditoria_qaoa_master.csv`: tabela mestre consolidada de todos os runs (CSV/JSON), incluindo colunas enriquecidas
- `auditoria_qaoa_energia.png` e `auditoria_qaoa_energia.html`: gráficos de energia por experimento
- `auditoria_qaoa_tempo.png` e `auditoria_qaoa_tempo.html`: gráficos de tempo por experimento
- `ambiente_execucao.json`: snapshot do ambiente (Python, SO, versões de pacotes)
- `manifest_codigo.json`: hashes SHA-256 dos scripts para rastreabilidade de versão

## Metodologia

1. **Execução dos experimentos**: geração de resultados por run (CSV/JSON)
2. **Enriquecimento dos CSVs** com metadados e métricas derivadas:
   - Metadados: Run Timestamp, Run ID, Qubits, P-layers, Shots, Max Iter, Seed
   - Energia Normalizada (%) = (Energia Final / Energia Máxima) × 100
   - Melhora vs Sem Ruído (%) = ((Energia Final − Energia Sem Ruído) / Energia Sem Ruído) × 100
   - Classificação: "Ruído benéfico" (> +1%), "prejudicial" (< −1%), "negligenciável" (−1% a +1%)
   - AUE (heurística): max(0, Melhora %) × (Energia Normalizada / 100)
   - TREX: Energia Final / Tempo (s)
3. **Consolidação multi-run** em `auditoria_qaoa_master.csv`
4. **Geração de gráficos** (energia e tempo)
5. **Rastreabilidade de versão** via SHA-256 dos scripts principais

## Rastreabilidade de Versão

### Hashes SHA-256

Cada run captura um manifest dos hashes SHA-256 dos scripts principais em `manifest_codigo.json`:

| Script | SHA-256 |
| --- | --- |
| experimento_qaoa_otimizado.py | [hash calculado] |
| enriquecer_resultados_qaoa.py | [hash calculado] |
| auditoria_qaoa_resultados.py | [hash calculado] |
| validar_auditoria_qaoa.py | [hash calculado] |
| framework_qaoa_100qubits.py | [hash calculado] |

### Verificação de Integridade

Para verificar se os scripts foram alterados desde a última auditoria:

```powershell
& ".\.venv\Scripts\python.exe" "calculador_hashes_qaoa.py"
```

Compare os hashes exibidos com os valores em `manifest_codigo.json`.

## Reprodutibilidade

- **Run ID**: baseado em timestamp (ex.: `run_20251228_145335`)
- **Metadados completos**: timestamp, ambiente, versões de pacotes
- **Hashes de código**: garantem integridade e versão do código utilizado
- **Scripts rastreados**:
  - `experimento_qaoa_otimizado.py`: geração de experimentos
  - `enriquecer_resultados_qaoa.py`: enriquecimento de CSVs
  - `auditoria_qaoa_resultados.py`: consolidação e gráficos
  - `validar_auditoria_qaoa.py`: validação QUALIS A1
  - `framework_qaoa_100qubits.py`: framework QAOA

## Como atualizar

### 1. Rodar experimentos

```powershell
& ".\.venv\Scripts\python.exe" "experimento_qaoa_otimizado.py"
```

Gera arquivos em `resultados_qaoa_otimizado/resultados_YYYYMMDD_HHMMSS.{csv,json}`

### 2. Enriquecer CSVs

```powershell
& ".\.venv\Scripts\python.exe" "enriquecer_resultados_qaoa.py"
```

Adiciona colunas de metadados e métricas; insere linha de resumo com hashes.

### 3. Consolidar e gerar gráficos

```powershell
& ".\.venv\Scripts\python.exe" "auditoria_qaoa_resultados.py"
```

Produz:
- `auditoria_qaoa_master.csv`
- `auditoria_qaoa_energia.{png,html}`
- `auditoria_qaoa_tempo.{png,html}`
- `ambiente_execucao.json`
- `manifest_codigo.json`

### 4. Validar QUALIS A1

```powershell
& ".\.venv\Scripts\python.exe" "validar_auditoria_qaoa.py"
```

Retorna `✅ Validação QUALIS A1: PASS` ou erros específicos.

## Observações

- CSVs errôneos (vazios) são mantidos com registro de erro para transparência
- Hashes SHA-256 inclusos no JSON de resumo de cada run
- Para experimentos maiores, ajuste `n_qubits`, `shots`, `max_iter` preservando tempo de execução
- Simulator Qiskit tem limite de ~30 qubits; QAOA atual usa 6 qubits (comprovadamente estável)
