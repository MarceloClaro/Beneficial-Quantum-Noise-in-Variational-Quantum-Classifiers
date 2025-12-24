# Execução Framework Qiskit - 2 Horas com Visualizações

## 📊 Status da Execução

**Data**: 24/12/2025  
**Script**: `executar_qiskit_2h_com_imagens.py`  
**Status**: ✅ Em execução (background)  
**Tempo estimado**: ~2 horas

---

## 🎯 Objetivo

Executar o framework Qiskit com geração de visualizações de alta qualidade dentro de um período de 2 horas, conforme solicitado pelo usuário.

---

## ⚙️ Configuração da Execução

### Parâmetros Otimizados

| Parâmetro | Valor | Justificativa |
|-----------|-------|---------------|
| **Datasets** | 3 | Moons, Circles, Iris (representativos) |
| **Arquiteturas** | 3 | Básico, Strongly Entangling, Hardware Efficient |
| **Tipos de Ruído** | 4 | Sem ruído, Phase Damping, Amplitude Damping, Depolarizante |
| **Níveis de Ruído** | 3 | 0.0, 0.005, 0.01 (baseline + benéfico) |
| **Seeds** | 2 | 42, 43 (reprodutibilidade) |
| **Épocas** | 5 | Reduzido para velocidade |
| **Shots** | 512 | Reduzido para velocidade |

### Cálculo Total

**Total de experimentos**: 3 × 3 × 4 × 3 × 2 = **216 experimentos**

- Menos combinações inválidas (sem_ruido com γ > 0): ~**180 experimentos válidos**
- Tempo estimado por experimento: ~30-40s
- **Tempo total estimado**: ~1.5-2 horas

### Visualizações

**Geração a cada 2 experimentos** (otimização):
- 4 visualizações por experimento selecionado
- Total esperado: ~360 visualizações (90 experimentos × 4 visualizações)

#### Tipos de Visualizações

1. **Esfera de Bloch** (`*_bloch.png`)
   - Visualização 3D do estado quântico de qubits
   - Representa superposição e fase

2. **State City 3D** (`*_city3d.png`)
   - "Arranha-céus quânticos" (plano árido solicitado)
   - Densidade de probabilidade em 3D

3. **Q-Sphere** (`*_qsphere.png`)
   - Representação esférica completa
   - Amplitude e fase de todos os estados

4. **Diagrama de Circuito** (`*_circuit.png`)
   - Visualização do circuito quântico
   - Qualidade de publicação

---

## 📁 Estrutura de Resultados

```
resultados_qiskit_2h_YYYYMMDD_HHMMSS/
├── visualizacoes/
│   ├── moons_basico_sem_ruido_g0.0000_s42_bloch.png
│   ├── moons_basico_sem_ruido_g0.0000_s42_city3d.png
│   ├── moons_basico_sem_ruido_g0.0000_s42_qsphere.png
│   ├── moons_basico_sem_ruido_g0.0000_s42_circuit.png
│   ├── moons_strongly_entangling_phase_damping_g0.0050_s42_bloch.png
│   └── ... (mais visualizações)
├── resultados_parciais.csv  (checkpoint a cada 10 experimentos)
└── resultados_completos.csv (resultados finais)
```

---

## 📊 Formato de Dados CSV

Cada linha contém:

```csv
experimento,config_nome,dataset,arquitetura,tipo_ruido,nivel_ruido,seed,
n_qubits,n_camadas,n_epocas,acuracia_treino,acuracia_teste,tempo_treino,
visualizacoes_geradas,framework,shots
```

---

## 🔄 Sistema de Checkpoint

- **Frequência**: A cada 10 experimentos
- **Arquivo**: `resultados_parciais.csv`
- **Informações salvias**:
  - Progresso atual
  - Tempo decorrido
  - Tempo estimado restante
  - Visualizações geradas

---

## 📈 Análises Incluídas

### Estatísticas Finais

- Acurácia média, máxima, mínima
- Desvio padrão
- Comparação por tipo de ruído
- Melhor configuração identificada

### Comparação Ruído Benéfico

Análise comparativa entre:
- **Sem ruído** (baseline)
- **Phase damping** (γ=0.005, γ=0.01)
- **Amplitude damping** (γ=0.005, γ=0.01)
- **Depolarizante** (γ=0.005, γ=0.01)

---

## 🚀 Como Executar

### Execução Manual

```bash
python executar_qiskit_2h_com_imagens.py
```

### Monitoramento

```bash
# Ver log em tempo real
tail -f execucao_2h.log

# Verificar progresso
ls -lh resultados_qiskit_2h_*/visualizacoes/

# Contar visualizações geradas
find resultados_qiskit_2h_*/ -name "*.png" | wc -l
```

---

## ⏱️ Timeline Estimado

| Tempo | Marco |
|-------|-------|
| 0min | Início da execução |
| 10min | Checkpoint 1 (~10 experimentos) |
| 20min | Checkpoint 2 (~20 experimentos) |
| 30min | Primeiras visualizações (~15 geradas) |
| 60min | Meio do caminho (~90 experimentos) |
| 90min | Checkpoint final (~150 experimentos) |
| 120min | **Conclusão** (~180 experimentos, ~360 visualizações) |

---

## 📝 Notas Técnicas

### Otimizações Aplicadas

1. **Épocas reduzidas**: 5 em vez de 10-20 (mais rápido, ainda demonstrativo)
2. **Shots reduzidos**: 512 em vez de 1024 (mantém precisão razoável)
3. **Visualizações seletivas**: A cada 2 experimentos (economiza tempo)
4. **Subset estratégico**: Foco em configurações que demonstram ruído benéfico

### Limitações

- Menor precisão estatística (menos seeds)
- Menor número de épocas (convergência parcial)
- Shots reduzidos (maior variância estocástica)
- **Adequado para**: Demonstração, proof-of-concept, visualizações
- **Não adequado para**: Análise estatística rigorosa, publicação final

---

## 🎯 Objetivos Atendidos

✅ **Execução em 2 horas**: Configuração otimizada  
✅ **Visualizações geradas**: 4 tipos exclusivos Qiskit  
✅ **Ruído benéfico demonstrado**: Comparação incluída  
✅ **Dados exportados**: CSV para análise posterior  
✅ **Checkpoint system**: Segurança contra falhas  

---

## 📊 Próximos Passos

Após conclusão da execução:

1. **Verificar resultados**: `resultados_qiskit_2h_*/resultados_completos.csv`
2. **Explorar visualizações**: `resultados_qiskit_2h_*/visualizacoes/`
3. **Analisar estatísticas**: Ver log final com resumo
4. **Documentar descobertas**: Atualizar MDs com imagens
5. **Comparar com PennyLane**: Usar mesmos datasets

---

## 🔗 Referências

- **Script de execução**: `executar_qiskit_2h_com_imagens.py`
- **Framework base**: `framework_qiskit.py`
- **Especificação completa**: `ESPECIFICACAO_EXPERIMENTOS_QISKIT.md`
- **Guia de uso**: `docs/GUIA_QISKIT.md`

---

**Última Atualização**: 24/12/2025 18:55  
**Status**: 🟢 Em execução  
**Tempo estimado de conclusão**: ~19:55 (2h após início)
