# Glossário de Termos - Geração de Artigos Científicos Qualis A1

## A

**Auditoria Técnica**: Processo sistemático de inventariar e verificar todos os componentes de um projeto de pesquisa (código, dados, configurações, logs, artefatos) para garantir que o que foi implementado corresponde exatamente ao que será reportado no artigo científico. Inclui contagem de linhas, identificação de dependências, mapeamento de funções, extração de parâmetros e validação de resultados.


## C

**Código→Método (Mapa)**: Tabela de rastreabilidade que conecta cada componente metodológico descrito no artigo (ex: definição do ansatz, canal de ruído, otimizador) com sua implementação exata no código-fonte (arquivo, função, linha), permitindo que revisores verifiquem a correspondência entre texto e implementação.


**Conivência Código-Texto**: Percentual de correspondência entre afirmações/números no texto do artigo e suas evidências no código/dados/logs. Meta: ≥95% (idealmente 100%). Medida através de verificação automatizada que cruza seções do artigo com arquivos de origem.


## E

**Effect Sizes (Tamanhos de Efeito)**: Métricas estatísticas que quantificam a magnitude de um efeito ou diferença entre grupos, independentemente do tamanho da amostra. Principais medidas:
- **Cohen's d**: Diferença padronizada entre médias (pequeno: 0.2, médio: 0.5, grande: 0.8)
- **Glass's Δ**: Similar a Cohen's d, mas usa desvio padrão do grupo controle
- **Hedges' g**: Versão corrigida de Cohen's d para amostras pequenas
- **η² (Eta quadrado)**: Proporção de variância explicada em ANOVA
- **r (Correlação de Pearson)**: Força de associação linear


## F

**[INFORMAÇÃO AUSENTE]**: Marcador utilizado quando uma informação que deveria existir não foi encontrada no código, dados ou documentação. Exemplo: versão de uma biblioteca não especificada, hiperparâmetro não documentado. Diferente de [NÃO DISPONÍVEL], indica que houve falha na documentação.


## I

**IMRAD**: Estrutura clássica de artigos científicos experimentais:
- **I**ntroduction (Introdução)
- **M**ethods (Metodologia/Métodos)
- **R**esults (Resultados)
- **A**nd (E)
- **D**iscussion (Discussão)


## L

**[LACUNA DE CITAÇÃO]**: Marcador utilizado (apenas em modo R0) quando falta uma referência bibliográfica para suportar uma afirmação, mas não é permitido adicionar novas referências. Indica necessidade de revisão ou de mudança para modo R1.


## M

**MODE_A vs MODE_B**:
- **MODE_A**: Artigo em INGLÊS, formatação LaTeX, estilo internacional (IEEE, Physical Review, Nature)
- **MODE_B**: Artigo em PORTUGUÊS, normas ABNT (NBR 10520 citações, NBR 6023 referências), para periódicos brasileiros


## N

**[NÃO DISPONÍVEL]**: Marcador utilizado quando uma informação não pode ser gerada ou obtida. Exemplo: resultados de um pipeline que não executa devido a erro, dados de experimento que não foi realizado. Diferente de [INFORMAÇÃO AUSENTE], indica impossibilidade de obtenção.


## P

**PROFILE_PR_QUANTUM vs PROFILE_GENERAL**:
- **PROFILE_PR_QUANTUM**: Tom técnico, matemático rigoroso, para Physical Review, Quantum, Nature Physics
- **PROFILE_GENERAL**: Tom mais narrativo, didático, para periódicos multidisciplinares ou de divulgação


## Q

**Quality Gate**: Ponto de verificação ao final de cada fase do processo de geração do artigo para garantir que critérios de qualidade foram atendidos antes de prosseguir para a próxima fase. Cada fase tem seus próprios critérios:
- **F1**: Rastreabilidade completa dos componentes técnicos
- **F2**: Lacuna operacionalizável e falsificável
- **F3**: Referências com DOI e contrapontos incluídos
- **F4**: Números com lastro verificável
- **F5**: Consistência experimental
- **Final**: Consistência ≥95% + reprodutibilidade completa


**Qualis A1**: Estrato superior do sistema Qualis de classificação de periódicos científicos no Brasil (CAPES). Periódicos A1 são considerados de excelência internacional, com alto fator de impacto e rigor na revisão por pares.


## R

**R0 vs R1 (Políticas de Referências)**:
- **R0** (Referências Travadas): Conjunto de referências é fixo e pré-aprovado. Não é permitido adicionar novas referências durante a geração. Se faltar uma citação, marca-se como [LACUNA DE CITAÇÃO]. Usado quando há restrições institucionais ou limitações de acesso a bases de dados.
- **R1** (Referências Expandidas): Permite buscar e adicionar novas referências durante a geração, desde que com DOI e justificativa. Busca organizada em 7 categorias: fundacionais, estado da arte, metodológicas, estatísticas, frameworks, críticas, aplicações.


**Rastreabilidade**: Capacidade de traçar cada afirmação, número ou resultado em um artigo científico de volta à sua origem exata no código, dados ou logs de execução. Implementada através de tabelas que mapeiam:
- Seção do artigo → Afirmação/número → Arquivo/função/linha → Referência bibliográfica


**Reprodutibilidade**: Capacidade de um terceiro replicar exatamente os resultados de um estudo. Requer documentação completa de:
- Ambiente computacional (OS, hardware, versões de software)
- Dependências (bibliotecas e versões exatas)
- Seeds aleatórias fixas
- Configurações de hiperparâmetros
- Pipeline de execução (comandos, ordem, flags)
- Dados de entrada (versão, pré-processamento)


## S

**Scope Conditions (Condições de Escopo)**: Condições sob as quais os resultados de um estudo são válidos e generalizáveis. Exemplo: "Resultados válidos para VQCs com ≤8 qubits e datasets sintéticos 2D". Essencial para delimitar claims e evitar overgeneralization.


## T

**Threats to Validity (Ameaças à Validade)**: Fatores que podem comprometer a validade das conclusões de um estudo. Categorias:
- **Validade Interna**: Ameaças à inferência causal (ex: confounders não controlados, seleção enviesada)
- **Validade Externa**: Ameaças à generalização (ex: características únicas da amostra, contexto específico)
- **Validade de Construto**: Ameaças à medição (ex: operacionalização inadequada de conceitos, instrumentos não validados)
- **Validade Estatística**: Ameaças à inferência estatística (ex: baixo poder, violação de pressupostos, múltiplas comparações sem correção)


## Siglas e Acrônimos

- **ABNT**: Associação Brasileira de Normas Técnicas
- **CARS**: Create A Research Space (modelo de introdução científica)
- **CI**: Confidence Interval (Intervalo de Confiança)
- **DOI**: Digital Object Identifier
- **IC**: Intervalo de Confiança (equivalente português de CI)
- **NBR**: Norma Brasileira
- **VQC**: Variational Quantum Classifier (Classificador Variacional Quântico)


## Símbolos Matemáticos Comuns

- **θ (theta)**: Parâmetros variacionais do circuito quântico
- **ψ (psi)**: Estado quântico
- **ρ (rho)**: Matriz densidade (estado misto)
- **U(θ)**: Operador unitário parametrizado
- **ℋ**: Espaço de Hilbert
- **𝒩**: Canal quântico (mapa de ruído)
- **L(θ)**: Função de custo (loss)
- **𝒟**: Dataset


## Referências de Formatação

### ABNT (MODE_B)
- **Citação autor-data**: (NIELSEN; CHUANG, 2010)
- **Referência completa**: NIELSEN, M. A.; CHUANG, I. L. **Quantum Computation and Quantum Information**. 10th Anniversary ed. Cambridge: Cambridge University Press, 2010.


### IEEE/APS (MODE_A)
- **Citação numérica**: [1]
- **Referência completa**: M. A. Nielsen and I. L. Chuang, *Quantum Computation and Quantum Information*, 10th Anniversary ed. Cambridge, UK: Cambridge University Press, 2010.


---


**Nota**: Este glossário é vivo e será atualizado conforme novos termos forem introduzidos no processo de geração de artigos científicos.

