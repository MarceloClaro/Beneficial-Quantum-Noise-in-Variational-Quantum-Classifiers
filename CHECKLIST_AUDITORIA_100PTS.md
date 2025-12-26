# Checklist de Auditoria para Artigos Científicos Qualis A1
## Sistema de Pontuação: 0-100 pontos

---

## CATEGORIA 1: REPRODUTIBILIDADE (30 pontos)

### 1.1 Ambiente Computacional Documentado (10 pontos)

- [ ] **[3 pts]** Sistema operacional especificado (ex: Ubuntu 22.04 LTS)
- [ ] **[2 pts]** Versão do Python documentada (ex: Python 3.9.7)
- [ ] **[2 pts]** Hardware descrito (CPU, GPU, memória RAM)
- [ ] **[2 pts]** Arquivo `requirements.txt` ou `environment.yml` presente e completo
- [ ] **[1 pt]** Versões exatas de todas as bibliotecas críticas (com `==`, não `>=`)

**Subtotal 1.1**: _____ / 10

---

### 1.2 Seeds Fixas e Reportadas (10 pontos)

- [ ] **[3 pts]** Seeds aleatórias documentadas no código
- [ ] **[2 pts]** Seeds reportadas no Methods (texto do artigo)
- [ ] **[2 pts]** Função de fixação de seed implementada corretamente (numpy, random, torch, etc.)
- [ ] **[2 pts]** Múltiplas seeds usadas para robustez (mínimo 3-5)
- [ ] **[1 pt]** Seeds justificadas (ex: valores "padrão" como 42, ou valores aleatórios documentados)

**Subtotal 1.2**: _____ / 10

---

### 1.3 Pipeline Executável (10 pontos)

- [ ] **[4 pts]** Script principal executa sem erros (testado em ambiente limpo)
- [ ] **[2 pts]** Comandos de execução documentados (README ou Methods)
- [ ] **[2 pts]** Tempo de execução estimado fornecido
- [ ] **[1 pt]** Logs de execução incluídos ou geráveis
- [ ] **[1 pt]** Dockerfile ou ambiente containerizado disponível (opcional, mas recomendado)

**Subtotal 1.3**: _____ / 10

---

**TOTAL CATEGORIA 1 (Reprodutibilidade)**: _____ / 30

---

## CATEGORIA 2: RASTREABILIDADE (30 pontos)

### 2.1 Tabela de Rastreabilidade Completa (15 pontos)

- [ ] **[5 pts]** Todas as afirmações quantitativas têm evidência rastreável
- [ ] **[4 pts]** Tabela de rastreabilidade `Seção → Evidência → Origem` preenchida
- [ ] **[3 pts]** Evidências são verificáveis (arquivos/linhas existem e são corretos)
- [ ] **[2 pts]** Nenhum número inventado ou sem lastro
- [ ] **[1 pt]** Marcadores [INFORMAÇÃO AUSENTE]/[NÃO DISPONÍVEL] usados quando apropriado

**Subtotal 2.1**: _____ / 15

---

### 2.2 Mapa Código→Método Completo (15 pontos)

- [ ] **[5 pts]** Tabela `Componente Método → Arquivo:Função:Linha` preenchida
- [ ] **[4 pts]** Todos os componentes metodológicos principais mapeados
- [ ] **[3 pts]** Parâmetros documentados (valores e justificativas)
- [ ] **[2 pts]** Artefatos gerados listados (figuras, tabelas, CSVs)
- [ ] **[1 pt]** Dependências de bibliotecas documentadas com versões

**Subtotal 2.2**: _____ / 15

---

**TOTAL CATEGORIA 2 (Rastreabilidade)**: _____ / 30

---

## CATEGORIA 3: RIGOR ESTATÍSTICO (20 pontos)

### 3.1 Testes Apropriados (5 pontos)

- [ ] **[2 pts]** Testes estatísticos escolhidos são adequados para os dados/hipóteses
- [ ] **[1 pt]** Pressupostos dos testes verificados (normalidade, homoscedasticidade)
- [ ] **[1 pt]** Testes paramétricos e não-paramétricos quando apropriado
- [ ] **[1 pt]** Justificativa para escolha dos testes fornecida

**Subtotal 3.1**: _____ / 5

---

### 3.2 Correção para Múltiplas Comparações (5 pontos)

- [ ] **[3 pts]** Correção aplicada quando há múltiplas comparações (Bonferroni, FDR, etc.)
- [ ] **[1 pt]** Tipo de correção documentado e justificado
- [ ] **[1 pt]** p-values ajustados reportados (não apenas p-values brutos)

**Subtotal 3.2**: _____ / 5

---

### 3.3 Intervalos de Confiança (5 pontos)

- [ ] **[3 pts]** Intervalos de confiança de 95% (ou outro nível justificado) reportados
- [ ] **[1 pt]** Método de cálculo documentado (bootstrap, t-distribution, etc.)
- [ ] **[1 pt]** ICs visualizados em figuras (barras de erro, ribbons)

**Subtotal 3.3**: _____ / 5

---

### 3.4 Tamanhos de Efeito (5 pontos)

- [ ] **[2 pts]** Tamanhos de efeito calculados (Cohen's d, η², r, etc.)
- [ ] **[2 pts]** Effect sizes reportados junto com p-values
- [ ] **[1 pt]** Interpretação dos tamanhos de efeito (pequeno/médio/grande)

**Subtotal 3.4**: _____ / 5

---

**TOTAL CATEGORIA 3 (Rigor Estatístico)**: _____ / 20

---

## CATEGORIA 4: TRANSPARÊNCIA (20 pontos)

### 4.1 Código Disponível Publicamente (10 pontos)

- [ ] **[5 pts]** Código completo disponível em repositório público (GitHub, GitLab, etc.)
- [ ] **[2 pts]** Repositório bem documentado (README detalhado)
- [ ] **[2 pts]** Licença de código aberto especificada (MIT, Apache, GPL)
- [ ] **[1 pt]** DOI ou identificador persistente para o código (Zenodo, figshare)

**Subtotal 4.1**: _____ / 10

---

### 4.2 Dados Disponíveis Publicamente (5 pontos)

- [ ] **[3 pts]** Dados brutos ou processados disponíveis publicamente
- [ ] **[1 pt]** Formato de dados bem documentado (schemas, codebooks)
- [ ] **[1 pt]** Licença de dados especificada (CC-BY, CC0, Open Data Commons)

*Nota: Se dados são gerados sinteticamente, documentar o script de geração*

**Subtotal 4.2**: _____ / 5

---

### 4.3 Limitações e Ameaças à Validade Discutidas (5 pontos)

- [ ] **[2 pts]** Seção "Threats to Validity" ou "Limitations" presente
- [ ] **[1 pt]** Ameaças à validade interna discutidas
- [ ] **[1 pt]** Ameaças à validade externa discutidas
- [ ] **[1 pt]** Scope conditions claramente especificadas

**Subtotal 4.3**: _____ / 5

---

**TOTAL CATEGORIA 4 (Transparência)**: _____ / 20

---

## PONTUAÇÃO FINAL

| Categoria | Pontos Obtidos | Pontos Máximos |
|-----------|----------------|----------------|
| 1. Reprodutibilidade | _____ | 30 |
| 2. Rastreabilidade | _____ | 30 |
| 3. Rigor Estatístico | _____ | 20 |
| 4. Transparência | _____ | 20 |
| **TOTAL** | **_____** | **100** |

---

## INTERPRETAÇÃO DA PONTUAÇÃO

### Classificação Qualis A1

| Pontuação | Classificação | Interpretação |
|-----------|---------------|---------------|
| 90-100 | 🥇 **Excelente** | Pronto para submissão a Nature, Science, Physical Review, Quantum |
| 75-89 | 🥈 **Muito Bom** | Pronto para submissão a periódicos Qualis A1/A2 com revisões menores |
| 60-74 | 🥉 **Bom** | Requer melhorias antes de submissão a Qualis A1 |
| 40-59 | ⚠️ **Insuficiente** | Requer melhorias substanciais |
| 0-39 | ❌ **Inadequado** | Não está pronto para submissão |

---

## AÇÕES RECOMENDADAS POR FAIXA

### 90-100 pontos (Excelente)
✅ **Pronto para submissão!**
- Revisar uma última vez com coautores
- Preparar cover letter destacando rigor metodológico
- Considerar submissão a periódicos de máximo impacto

### 75-89 pontos (Muito Bom)
⚠️ **Pequenos ajustes necessários**
- Revisar seções que perderam pontos
- Completar [INFORMAÇÃO AUSENTE] se houver
- Adicionar elementos faltantes (IC, effect sizes, etc.)
- Estimar 1-2 semanas para melhorias

### 60-74 pontos (Bom)
🔧 **Melhorias moderadas necessárias**
- Priorizar rastreabilidade e reprodutibilidade
- Adicionar análises estatísticas faltantes
- Melhorar documentação de código/dados
- Estimar 2-4 semanas para melhorias

### 40-59 pontos (Insuficiente)
🚧 **Trabalho substancial necessário**
- Revisar fundamentalmente a metodologia
- Adicionar seeds, logs, documentação
- Refazer análises estatísticas com rigor
- Estimar 1-2 meses para melhorias

### 0-39 pontos (Inadequado)
🛑 **Repensar abordagem**
- Revisar se o estudo está maduro para publicação
- Considerar experimentos adicionais
- Buscar mentoria/consultoria metodológica
- Estimar 3+ meses para melhorias

---

## REQUISITOS MÍNIMOS POR PERIÓDICO

### Nature, Science (Cell Press)
- **Mínimo**: 90 pontos
- **Ênfase**: Reprodutibilidade (30/30), Transparência (20/20)
- **Extras**: Código depositado em repositório persistente, Statement de Disponibilidade de Dados

### Physical Review (APS), Quantum
- **Mínimo**: 85 pontos
- **Ênfase**: Rigor Estatístico (18/20), Rastreabilidade (25/30)
- **Extras**: Equações correspondem ao código, notação precisa

### Periódicos Brasileiros Qualis A1
- **Mínimo**: 75 pontos
- **Ênfase**: Reprodutibilidade (25/30), Normas ABNT
- **Extras**: Resumo em português e inglês, palavras-chave de tesauro

---

## CHECKLIST ADICIONAL DE CONFORMIDADE EDITORIAL

### Para MODO A (Inglês/Internacional)

- [ ] Texto em inglês revisado por native speaker
- [ ] Referências no formato do periódico (APS, IEEE, Nature style)
- [ ] Figuras em alta resolução (≥300 DPI)
- [ ] Equações numeradas e referenciadas corretamente
- [ ] Unidades SI usadas consistentemente

### Para MODO B (Português/ABNT)

- [ ] Texto em português formal (não coloquial)
- [ ] Citações autor-data conforme NBR 10520
- [ ] Referências completas conforme NBR 6023
- [ ] Resumo ≤250 palavras
- [ ] Abstract alinhado com Resumo

---

## SCRIPT DE AUTO-AVALIAÇÃO

```python
#!/usr/bin/env python3
"""
auto_avaliacao_qualis.py
Calcula automaticamente pontuação em itens verificáveis.
"""

import os
import re
from pathlib import Path

def avaliar_repositorio(repo_path):
    """Avalia automaticamente itens objetivos."""
    pontos = 0
    
    # 1.1: Ambiente documentado
    if (Path(repo_path) / 'requirements.txt').exists():
        pontos += 2
    if (Path(repo_path) / 'README.md').exists():
        pontos += 2
    
    # 1.2: Seeds no código
    for py_file in Path(repo_path).rglob('*.py'):
        with open(py_file) as f:
            if 'seed' in f.read().lower():
                pontos += 3
                break
    
    # 4.1: Código público (verificar se está em Git)
    if (Path(repo_path) / '.git').exists():
        pontos += 5
    
    return pontos

if __name__ == '__main__':
    repo = '.'
    pts = avaliar_repositorio(repo)
    print(f"Pontuação automática (parcial): {pts}/100")
    print("Complete a avaliação manual para pontuação total.")
```

---

## REGISTRO DE AUDITORIA

**Data da Auditoria**: _____________  
**Auditor**: _____________________  
**Versão do Artigo**: _____________  
**Pontuação Total**: _____ / 100  
**Classificação**: ________________  

**Observações**:
```
[Espaço para notas do auditor]




```

**Assinatura**: _____________________  
**Data**: ___________________________

---

**Template Version**: 1.0  
**Última Atualização**: 26/12/2025  
**Compatível com**: Qualis 2021-2024, APS, Nature, IEEE standards
