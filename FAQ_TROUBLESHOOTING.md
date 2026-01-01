# FAQ e Troubleshooting - Geração de Artigos Científicos

## Perguntas Frequentes (FAQ)

### Configuração e Planejamento

**P: Qual modo devo escolher: MODE_A ou MODE_B?**


#### R:
- **MODE_A**: Para submissão a periódicos internacionais (Nature, Science, Physical Review, Quantum, IEEE). Artigo em INGLÊS com formatação LaTeX.
- **MODE_B**: Para submissão a periódicos brasileiros ou lusófonos. Artigo em PORTUGUÊS com normas ABNT (NBR 10520/6023).


Escolha com base no periódico-alvo primário listado em `config.json`.

---


**P: Quando usar R0 vs R1?**


#### R:
- **R0 (Referências Travadas)**: Use quando:
  - Há restrições institucionais sobre referências
  - Trabalha com um conjunto pré-aprovado de citações
  - Tem acesso limitado a bases de dados acadêmicas
  - Quer garantir que apenas referências específicas sejam usadas

  
- **R1 (Referências Expandidas)**: Use quando:
  - Tem acesso completo a bases de dados (Web of Science, Scopus, arXiv)
  - Precisa encontrar referências para técnicas recentes
  - Quer maximizar a fundamentação bibliográfica
  - Aceita adicionar referências durante o processo


**Recomendação**: R1 é geralmente preferível para artigos de alto impacto.


---


**P: O que significa PROFILE_PR_QUANTUM vs PROFILE_GENERAL?**


#### R:
- **PROFILE_PR_QUANTUM**: Tom técnico, matemática rigorosa, notação avançada. Para Physical Review, Quantum, Nature Physics, periódicos de física/computação quântica.
- **PROFILE_GENERAL**: Tom mais acessível, explicações didáticas, ênfase em contexto e impacto. Para periódicos multidisciplinares ou de audiência mais ampla.


---


### Problemas Comuns e Soluções

**P: O código não tem seeds fixas. O que fazer?**


**R:**
1. **Documentar**: Marque no `analise_codigo_inicial.md` como `[INFORMAÇÃO AUSENTE: Seeds aleatórias não fixadas]`
2. **Adicionar às Threats to Validity**: Na seção de Discussão, inclua:

   ```text

   **Threats to Validity - Reprodutibilidade Estocástica**: O código atual não

   utiliza seeds fixas para geradores de números aleatórios, o que limita a
   reprodutibilidade bit-a-bit dos resultados.
   ```

3. **Executar múltiplas vezes**: Execute o pipeline 10-30 vezes e reporte:
   - Média e desvio padrão das métricas
   - Intervalos de confiança de 95%
   - Distribuições dos resultados
4. **Trabalho Futuro**: Mencione na Conclusão que versões futuras devem incluir seeds fixas


**Script auxiliar**:

```bash

# Execute 10 vezes e agregue resultados
for i in {1..10}; do
    python seu_pipeline.py --output results_$i.json
done
python agregar_resultados.py results_*.json --output resultados_agregados.csv

```text

---


**P: Como proceder se não há logs de execução?**


**R:**


**Opção 1 (Preferível)**: Executar e gerar logs

```bash

# Execute pipeline com logging completo
python seu_pipeline.py --verbose --log-file execucao.log 2>&1 | tee console.log

```text

**Opção 2**: Se execução for inviável (muito tempo/recursos):
1. Marque todos os números como `[NÃO DISPONÍVEL]`
2. Foque o artigo na **metodologia** e **design experimental**
3. Apresente resultados esperados ou simulados como **plano experimental**
4. Seja explícito: "Os resultados apresentados são projeções baseadas em execuções parciais"


**Opção 3**: Use resultados de execuções anteriores:

```bash

# Procure logs/resultados em outras pastas
find . -name "*.log" -o -name "results*.json" -o -name "*.csv"

```text

---


**P: Quando usar [INFORMAÇÃO AUSENTE] vs [NÃO DISPONÍVEL] vs [LACUNA DE CITAÇÃO]?**


**R:**


| Marcador | Uso | Exemplo | Ação |
|----------|-----|---------|------|
| `[INFORMAÇÃO AUSENTE]` | Info deveria existir mas não foi encontrada | Versão de biblioteca não especificada | Investigar código/docs e adicionar |
| `[NÃO DISPONÍVEL]` | Info não pode ser gerada/obtida | Resultados de pipeline que não executa | Executar pipeline ou remover claim |
| `[LACUNA DE CITAÇÃO]` | Falta referência (apenas R0) | Técnica sem citação em modo R0 | Mudar para R1 ou reformular texto |

**Regra de ouro**: Nunca invente informação. Se não sabe, marque explicitamente.


---


**P: O total de configurações experimentais não bate. Como calcular corretamente?**


**R:**


**Fórmula**: Total = Fator1 × Fator2 × ... × FatorN × Repetições


**Exemplo Beneficial Quantum Noise**:

```

Total = 4 datasets × 7 ansätze × 6 tipos de ruído × 2 schedules × 8 inicializações × 5 seeds
      = 4 × 7 × 6 × 2 × 8 × 5
      = 13.440 configurações

```text

**Validação no código**:

```python
import itertools

datasets = ['moons', 'circles', 'blobs', 'iris']
ansatze = ['BasicEntangler', 'StronglyEntangling', ...]
noise_types = ['Depolarizing', 'AmplitudeDamping', ...]
schedules = ['Constant', 'Linear']
inits = ['random', 'zeros', ...]
seeds = [42, 123, 456, 789, 1024]

total = len(list(itertools.product(datasets, ansatze, noise_types, schedules, inits, seeds)))
print(f"Total configurações: {total}")

```text

---


**P: Como adaptar para áreas fora de Computação Quântica?**


**R:**


O framework é **agnóstico ao domínio**. Substitua:

| Quantum Computing | Aprendizado de Máquina | Bioinformática | Ciências Sociais |
|-------------------|------------------------|----------------|------------------|
| Ansätze | Arquiteturas de rede | Algoritmos de alinhamento | Modelos estatísticos |
| Qubits | Neurônios/camadas | Sequências | Variáveis/fatores |
| Canais de ruído | Dropout/regularização | Mutações/erros | Viés/confounders |
| Circuito quântico | Pipeline de ML | Pipeline bioinformático | Modelo causal |

**Exemplo adaptado (ML)**:

```json
{
  "models": ["ResNet50", "VGG16", "EfficientNet"],
  "datasets": ["CIFAR-10", "ImageNet", "MNIST"],
  "regularization": ["Dropout", "L2", "BatchNorm"],
  "optimizers": ["Adam", "SGD", "RMSprop"]
}

```text

---


### Problemas Técnicos

**P: O `gerador_artigo_completo.py` falha ao analisar o código. O que fazer?**


**R:**


**Diagnóstico**:

```bash
python gerador_artigo_completo.py --repositorio . --output teste_debug --verbose

```text

**Problemas comuns**:


1. **Erro ao ler `requirements.txt`**:

   ```bash

   # Verifique se existe e está bem formatado
   cat requirements.txt

   # Formato esperado: biblioteca==versão (uma por linha)
   ```text

2. **Erro ao encontrar arquivos Python**:

   ```bash

   # Verifique estrutura
   find . -name "*.py" | head -20
   ```text

3. **Encoding issues**:

   ```python

   # Adicione ao código
   with open(arquivo, 'r', encoding='utf-8') as f:
       conteudo = f.read()
   ```text

---


**P: A Fase X falhou no Quality Gate. Como corrigir?**


**R:**


**Quality Gate F1 (Auditoria)**:
- ✅ Toda listagem tem origem clara (arquivo:linha)
- ✅ Total de configurações calculado
- ✅ Nenhuma afirmação sem [MARCADOR] quando aplicável


**Correção**: Revise `analise_codigo_inicial.md` e adicione origens faltantes.


---


**Quality Gate F2 (Enquadramento)**:
- ✅ Pergunta de pesquisa é clara e testável
- ✅ Lacuna é operacionalizável (não vaga)


**Correção**: Revise `linha_de_pesquisa.md` e seja mais específico.


---


**Quality Gate F3 (Bibliografia)**:
- ✅ Cada referência tem DOI (quando disponível)
- ✅ Toda técnica central tem citação
- ✅ Contrapontos incluídos (não só papers favoráveis)


**Correção**:

```bash

# Busque DOIs faltantes
grep "\[DOI:" referencias_compiladas.md | grep "N/A"

# Para cada um, busque em <https://doi.org> ou CrossRef

```text

---


**Quality Gate F4 (Redação)**:
- ✅ Nenhum número sem lastro
- ✅ Tom consistente (MODE_A ou MODE_B)
- ✅ Referências formatadas corretamente


**Correção**: Use a ferramenta de consistência:

```bash
python tools/verificar_consistencia.py metodologia_completa.md --origem analise_codigo_inicial.md

```text

---


**Quality Gate F5 (Suplementar)**:
- ✅ Tabela S1 tem total correto de linhas
- ✅ Todas as figuras têm descrição detalhada


**Correção**: Conte linhas e compare:

```bash
wc -l tabela_s1_configuracoes.csv

# Deve ser: total_configs + 1 (header)

```text

---


**Quality Gate Final**:
- ✅ Consistência ≥ 95%
- ✅ Citações ↔ Referências 100%
- ✅ Ambiente reprodutível documentado


**Correção**: Execute auditoria completa:

```bash
python tools/auditoria_final.py --artigo fase6_consolidacao/artigo_abnt_final.md

```text

---


### Performance e Otimização

**P: A geração está muito lenta. Como acelerar?**


**R:**


1. **Use modo paralelo** (se disponível):

   ```bash
   python gerador_artigo_completo.py --repositorio . --output artigo --jobs 4
   ```text

2. **Pule fases já completas**:

   ```bash
   python gerador_artigo_completo.py --skip-fase 1,2 --start-from 3
   ```text

3. **Desabilite análise profunda**:

   ```bash
   python gerador_artigo_completo.py --shallow-analysis
   ```text

---


**P: Meu repositório é enorme (>10GB). Como lidar?**


**R:**


1. **Ignore pastas grandes** (adicione ao `.gitignore` ou ao config):

   ```json
   {
     "exclude_patterns": [
       "data/raw/*",
       "checkpoints/*",
       "*.h5",
       "*.pkl"
     ]
   }
   ```text

2. **Analise apenas arquivos relevantes**:

   ```json
   {
     "include_patterns": [
       "src/**/*.py",
       "experiments/*.py",
       "models/*.py"
     ]
   }
   ```text

---


### Formatação e Estilo

**P: Como garantir que ABNT está correto?**


**R:**


**Validação automática**:

```bash

# Instale validador
pip install abntex2-validator

# Valide referências
abntex2-validator agradecimentos_referencias.md

```text

**Checklist manual** (NBR 6023):
- [ ] Autor em CAIXA ALTA
- [ ] Título em **negrito**
- [ ] Ano entre parênteses ou após vírgula
- [ ] DOI no final


**Exemplo correto**:

```

NIELSEN, M. A.; CHUANG, I. L. **Quantum Computation and Quantum Information**.
10th Anniversary ed. Cambridge: Cambridge University Press, 2010.
DOI: 10.1017/CBO9780511976667

```text

---


**P: Como converter de MODE_B (ABNT/Markdown) para MODE_A (LaTeX)?**


**R:**


**Conversão automática** (use o conversor incluído):

```bash
python tools/converter_abnt_para_latex.py \

    --input artigo_abnt_final.md \
    --output artigo_latex.tex \
    --template ieee  # ou: physical_review, nature

```text

**Conversão manual**: Use Pandoc

```bash
pandoc artigo_abnt_final.md -o artigo.tex --template=template_physrev.tex

```text

---


### Integração e Workflows

**P: Como integrar com Overleaf?**


**R:**


1. **Gere LaTeX** (MODE_A):

   ```bash
   python gerador_artigo_completo.py --mode A --output artigo_latex
   ```text

2. **Crie projeto Overleaf**:
   - Novo projeto → Upload Project
   - Faça upload de `artigo_latex/`


3. **Sincronize com GitHub**:
   - Overleaf → Menu → GitHub → Link repo
   - Push/pull automáticos


---


**P: Como automatizar a geração em CI/CD?**


**R:**


**GitHub Actions** (`.github/workflows/gerar_artigo.yml`):

```yaml
name: Gerar Artigo
on:
  push:
    branches: [main]
jobs:
  gerar:
    runs-on: ubuntu-latest
    steps:

      - uses: actions/checkout@v3
      - name: Setup Python

        uses: actions/setup-python@v4
        with:
          python-version: '3.9'

      - name: Install dependencies

        run: pip install -r requirements.txt

      - name: Gerar artigo

        run: python gerador_artigo_completo.py --repositorio . --output artigo_gerado

      - name: Upload artifact

        uses: actions/upload-artifact@v3
        with:
          name: artigo-cientifico
          path: artigo_gerado/

```text

---


## Solução de Problemas Avançados

### Problema: "ImportError: No module named X"

**Causa**: Dependência faltando ou versão incompatível.


**Solução**:

```bash

# Reinstale dependências
pip install -r requirements.txt --upgrade

# Se persistir, crie ambiente limpo
python -m venv venv_artigo
source venv_artigo/bin/activate  # Windows: venv_artigo\Scripts\activate
pip install -r requirements.txt

```text

---


### Problema: "Memory Error ao processar grandes logs"

**Solução**: Processe em chunks

```python

# No código
def processar_log_grande(arquivo_log, chunk_size=10000):
    with open(arquivo_log, 'r') as f:
        while True:
            linhas = list(itertools.islice(f, chunk_size))
            if not linhas:
                break
            processar_chunk(linhas)

```text

---


### Problema: "Encoding error ao ler arquivo"

**Solução**:

```python

# Tente diferentes encodings
for encoding in ['utf-8', 'latin-1', 'cp1252']:
    try:
        with open(arquivo, 'r', encoding=encoding) as f:
            conteudo = f.read()
        break
    except UnicodeDecodeError:
        continue

```

---


## Contato e Suporte

**Issues GitHub**: Para bugs e solicitações de features, abra issue em:

<https://github.com/MarceloClaro/Beneficial-Quantum-Noise-in-Variational-Quantum-Classifiers/issues>

**Documentação**: Consulte sempre a documentação principal:
- README.md
- GERADOR_ARTIGO_README.md
- GLOSSARIO.md


---


**Última atualização**: 26/12/2025

