# Kit de Reprodutibilidade e Validação — Banca CNPq / Qualis A1

Este diretório organiza o projeto em um fluxo **profissional, auditável e didático** para banca e submissão Qualis A1.

## Objetivo

Garantir que qualquer avaliador consiga:
1. Reproduzir o experimento sem ambiguidade.
2. Validar os artefatos estatísticos e de rastreabilidade.
3. Verificar critérios objetivos de aceite metodológico.

## Estrutura

- `01_GUIA_RAPIDO.md`: guia de execução em 15–30 min.
- `02_PROCEDIMENTO_COMPLETO.md`: protocolo detalhado para auditoria.
- `03_CHECKLIST_BANCA.md`: checklist formal de validação.
- `04_MODELO_RELATORIO_BANCA.md`: template de relatório técnico.
- `scripts/run_reproducibility.sh`: execução padronizada ponta a ponta.
- `scripts/validate_artifacts.py`: verificação automática de artefatos e quality gates.
- `templates/ATA_VALIDACAO.md`: modelo de ata para evidência institucional.
- `examples/`: exemplos de chamadas e saídas esperadas.

## Fluxo recomendado

1. Seguir `01_GUIA_RAPIDO.md`.
2. Rodar `scripts/run_reproducibility.sh`.
3. Rodar `scripts/validate_artifacts.py` no diretório de resultados.
4. Preencher `03_CHECKLIST_BANCA.md` e `templates/ATA_VALIDACAO.md`.
5. Consolidar no `04_MODELO_RELATORIO_BANCA.md`.

## Resultado esperado

A banca terá evidência de:
- Execução reproduzível (seed, ambiente, protocolo).
- Qualidade estatística (p-valor, IC, permutação, estabilidade).
- Integridade de arquivos (checksums).
- Documentação clara para publicação A1.
