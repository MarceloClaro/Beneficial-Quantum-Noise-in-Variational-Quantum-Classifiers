# Melhorias para Análise e Busca de Melhor Performance

**Data:** 24 de dezembro de 2025  
**Framework Version:** 7.2  
**Objetivo:** Otimizar a busca de hiperparâmetros e melhorar análise de performance


---


## 🎯 Melhorias Disponíveis

### 1. Otimização Bayesiana Avançada (Já Implementado)

O framework já usa **Tree-structured Parzen Estimator (TPE)** da biblioteca Optuna, que é estado-da-arte em otimização Bayesiana.

**Configuração Atual:**

```python
sampler = TPESampler(
    seed=42,
    n_startup_trials=10  # Trials aleatórios iniciais para exploração
)
pruner = MedianPruner(
    n_startup_trials=5,
    n_warmup_steps=3  # Interrompe trials ruins após 3 épocas
)

```text

**Melhorias Possíveis:**


#### 1.1. Aumentar Trials Iniciais de Exploração

```bash

# Modificar no código (linha 3959):
n_startup_trials=20  # Mais exploração inicial (padrão: 10)

```text

**Benefício:** Melhor cobertura do espaço de hiperparâmetros antes de focar em regiões promissoras.


#### 1.2. Ajustar Estratégia de Pruning

```bash

# Modificar no código (linha 3960):
n_warmup_steps=5  # Aguardar mais épocas antes de podar (padrão: 3)

```text

**Benefício:** Evita descartar configurações que começam mal mas melhoram depois.


---


### 2. Usar Samplers Alternativos

#### 2.1. CMA-ES Sampler (Melhor para Espaços Contínuos)

```python
from optuna.samplers import CmaEsSampler

sampler = CmaEsSampler(
    seed=42,
    sigma0=0.1  # Desvio padrão inicial
)

```text

**Quando usar:** Ideal quando os hiperparâmetros são principalmente contínuos (learning rate, noise level).


#### Vantagens:
- ✅ Melhor para otimização de parâmetros contínuos
- ✅ Convergência mais rápida em espaços suaves
- ✅ Robusto a ruído


#### 2.2. NSGAII Sampler (Multi-Objetivo)

```python
from optuna.samplers import NSGAIISampler

sampler = NSGAIISampler(
    population_size=50,
    mutation_prob=0.1,
    crossover_prob=0.9
)

```text

**Quando usar:** Otimizar múltiplos objetivos simultaneamente (e.g., acurácia + tempo de treinamento).


**Exemplo de uso:**

```python
def objective_multi(trial):

    # ... configurar VQC ...
    accuracy = vqc.score(X_test, y_test)
    training_time = vqc.tempo_treinamento
    
    # Retornar múltiplos objetivos
    return accuracy, -training_time  # Maximizar acurácia, minimizar tempo

```text

---


### 3. Hyperband Pruning (Mais Agressivo)

```python
from optuna.pruners import HyperbandPruner

pruner = HyperbandPruner(
    min_resource=1,        # Mínimo de épocas
    max_resource=20,       # Máximo de épocas
    reduction_factor=3     # Fator de redução
)

```text

**Benefício:** Descarta configurações ruins muito mais rápido, economizando tempo computacional.


---


### 4. Análise de Importância Aprimorada

#### 4.1. Importância via fANOVA (Já Implementado)

```python
importances = optuna.importance.get_param_importances(
    study,
    evaluator=optuna.importance.FanovaImportanceEvaluator()
)

```text

#### 4.2. Adicionar Importância via Permutação

```python
from optuna.importance import MeanDecreaseImpurityImportanceEvaluator

importances_mdi = optuna.importance.get_param_importances(
    study,
    evaluator=MeanDecreaseImpurityImportanceEvaluator()
)

```text

**Benefício:** Múltiplas métricas de importância para validação cruzada.


---


### 5. Busca em Duas Fases (Recomendado)

**Fase 1: Exploração Global (Bayesian)**

```bash

# 1. Encontrar região promissora (1-2 horas)
export VQC_BAYESIAN=1
export VQC_QUICK=1
python framework_investigativo_completo.py --bayes --trials 200 --dataset-bayes all

```text

**Fase 2: Refinamento Local (Grid Search Focado)**

```python

# 2. Grid search refinado na região ótima encontrada
# Modificar ranges baseado nos melhores trials:

# Se melhor foi: depolarizante, nivel=0.001, strongly_entangling, quantico
# Criar grid focado:
noise_levels_refined = [0.0005, 0.001, 0.0015, 0.002, 0.0025]
architectures_refined = ['strongly_entangling', 'hardware_efficient']
init_strategies_refined = ['quantico', 'matematico']

```text

**Benefício:** Combina eficiência da busca Bayesiana com rigor do grid search.


---


### 6. Early Stopping Inteligente

```python

# Modificar ClassificadorVQC para incluir:
early_stopping_patience = 10  # Parar se não melhorar por 10 épocas
early_stopping_delta = 0.001   # Melhoria mínima considerada significativa

def fit(self, X, y, X_val=None, y_val=None):
    best_loss = float('inf')
    patience_counter = 0
    
    for epoca in range(n_epocas):

        # ... treinamento ...
        
        if X_val is not None:
            val_loss = self.compute_loss(X_val, y_val)
            
            if val_loss < best_loss - early_stopping_delta:
                best_loss = val_loss
                patience_counter = 0
            else:
                patience_counter += 1
                
            if patience_counter >= early_stopping_patience:
                logger.info(f"Early stopping na época {epoca}")
                break

```text

**Benefício:** Economiza tempo de treinamento sem perder qualidade.


---


### 7. Cross-Validation K-Fold

```python
from sklearn.model_selection import KFold

def evaluate_with_cv(config, X, y, n_folds=5):
    """Avalia configuração com validação cruzada"""
    kf = KFold(n_splits=n_folds, shuffle=True, random_state=42)
    scores = []
    
    for fold, (train_idx, val_idx) in enumerate(kf.split(X)):
        X_train, X_val = X[train_idx], X[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]
        
        vqc = ClassificadorVQC(**config)
        vqc.fit(X_train, y_train)
        score = vqc.score(X_val, y_val)
        scores.append(score)
    
    return np.mean(scores), np.std(scores)

```text

**Benefício:** Estimativa mais robusta da performance real, reduz overfitting na validação.


---


### 8. Ensemble de Modelos

```python
def ensemble_predict(models, X):
    """Predição por votação majoritária"""
    predictions = np.array([model.predict(X) for model in models])
    
    # Votação majoritária
    from scipy.stats import mode
    ensemble_pred = mode(predictions, axis=0)[0].flatten()
    
    return ensemble_pred

# Treinar múltiplos modelos com melhores configurações
best_configs = study.best_trials[:5]  # Top 5 configurações
models = []

for config in best_configs:
    vqc = ClassificadorVQC(**config.params)
    vqc.fit(X_train, y_train)
    models.append(vqc)

# Predição ensemble
y_pred_ensemble = ensemble_predict(models, X_test)
accuracy_ensemble = accuracy_score(y_test, y_pred_ensemble)

```text

**Benefício:** Melhora robustez e pode superar modelos individuais.


---


### 9. Análise de Performance Detalhada

#### 9.1. Curvas de Aprendizado

```python
def plot_learning_curves(vqc, X_train, y_train, X_test, y_test):
    """Plota curvas de aprendizado para diagnóstico"""
    train_sizes = np.linspace(0.1, 1.0, 10)
    train_scores = []
    test_scores = []
    
    for size in train_sizes:
        n_samples = int(len(X_train) * size)
        X_subset = X_train[:n_samples]
        y_subset = y_train[:n_samples]
        
        vqc_temp = ClassificadorVQC(**vqc.get_params())
        vqc_temp.fit(X_subset, y_subset)
        
        train_scores.append(vqc_temp.score(X_subset, y_subset))
        test_scores.append(vqc_temp.score(X_test, y_test))
    
    plt.figure(figsize=(10, 6))
    plt.plot(train_sizes * len(X_train), train_scores, label='Treino')
    plt.plot(train_sizes * len(X_train), test_scores, label='Teste')
    plt.xlabel('Número de Amostras de Treino')
    plt.ylabel('Acurácia')
    plt.title('Curvas de Aprendizado')
    plt.legend()
    plt.grid(True)
    plt.savefig('learning_curves.png', dpi=300)

```text

#### 9.2. Matriz de Confusão Detalhada

```python
from sklearn.metrics import confusion_matrix, classification_report
import seaborn as sns

def detailed_performance_analysis(y_true, y_pred):
    """Análise detalhada de performance"""

    # Matriz de confusão
    cm = confusion_matrix(y_true, y_pred)
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
    plt.title('Matriz de Confusão')
    plt.ylabel('Verdadeiro')
    plt.xlabel('Predito')
    plt.savefig('confusion_matrix.png', dpi=300)
    
    # Relatório de classificação
    report = classification_report(y_true, y_pred)
    print(report)
    
    return cm, report

```text

---


### 10. Paralelização de Trials

```python

# No código de otimização Bayesiana:
study.optimize(
    objective,
    n_trials=200,
    n_jobs=4,  # NOVO: Executar 4 trials em paralelo
    show_progress_bar=True
)

```text

**Benefício:** Reduz tempo total de execução em 4x (com 4 cores).


**Requisito:** Processador multi-core disponível.


---


## 🚀 Recomendações Práticas

### Para Pesquisa Exploratória (1-2 horas)

```bash

# Usar Bayesian otimizado com mais trials
export VQC_BAYESIAN=1
export VQC_QUICK=1
python framework_investigativo_completo.py \

    --bayes \
    --trials 500 \
    --dataset-bayes all

```text

### Para Máxima Performance (4-6 horas)

```bash

# Fase 1: Bayesian (200 trials, 2 horas)
python framework_investigativo_completo.py \

    --bayes \
    --trials 200 \
    --dataset-bayes all


# Fase 2: Grid Search Refinado na região ótima (2-4 horas)
# Ajustar ranges baseado nos resultados da Fase 1
python framework_investigativo_completo.py \

    --grid-focused \
    --focus-region "best_from_phase1"

```text

### Para Publicação Científica (15-20 horas)

```bash

# Grid Search Completo
python framework_investigativo_completo.py

```text

---


## 📊 Métricas de Performance Adicionais

### 1. AUC-ROC (Area Under Curve)

```python
from sklearn.metrics import roc_auc_score, roc_curve

# Para classificação binária
y_proba = vqc.predict_proba(X_test)
auc = roc_auc_score(y_test, y_proba)

```text

### 2. F1-Score (Balanceamento Precision/Recall)

```python
from sklearn.metrics import f1_score

f1 = f1_score(y_test, y_pred, average='weighted')

```text

### 3. Matthews Correlation Coefficient (MCC)

```python
from sklearn.metrics import matthews_corrcoef

mcc = matthews_corrcoef(y_test, y_pred)

```text

---


## 💡 Sugestões Específicas para Este Framework

### 1. Adicionar Análise de Sensibilidade

```python
def sensitivity_analysis(base_config, param_name, param_range):
    """Analisa sensibilidade a um hiperparâmetro específico"""
    results = []
    
    for value in param_range:
        config = base_config.copy()
        config[param_name] = value
        
        vqc = ClassificadorVQC(**config)
        vqc.fit(X_train, y_train)
        score = vqc.score(X_test, y_test)
        
        results.append({'value': value, 'score': score})
    
    return pd.DataFrame(results)

# Exemplo de uso:
noise_sensitivity = sensitivity_analysis(
    base_config=best_config,
    param_name='nivel_ruido',
    param_range=np.linspace(0, 0.02, 20)
)

```text

### 2. Implementar AutoML com AutoGluon

```python

# Alternativa: Usar AutoGluon para busca automática
from autogluon.tabular import TabularPredictor

# Criar dataset com features extraídas do quantum circuit
X_quantum_features = extract_quantum_features(X)

predictor = TabularPredictor(label='target').fit(
    train_data=X_quantum_features,
    time_limit=3600  # 1 hora
)

```text

---


## 📈 Benchmarks Esperados

| Método | Trials/Experimentos | Tempo | Acurácia Esperada | Uso |
|--------|---------------------|-------|-------------------|-----|
| Quick Bayesian (atual) | 5 | 7 min | 80-82% | Validação rápida |
| Bayesian Otimizado | 200 | 1-2h | 85-90% | Pesquisa exploratória |
| Bayesian Avançado | 500 | 4-6h | 88-92% | Máxima performance |
| Grid Search Completo | 8,280 | 15-20h | 90-95% | Publicação científica |
| Híbrido (2 fases) | 200+50 | 3-4h | 90-93% | Melhor custo-benefício |

---


## ✅ Checklist de Implementação

Para implementar as melhorias:

- [ ] Aumentar `n_trials` para 200-500
- [ ] Ajustar `n_startup_trials` para 20
- [ ] Implementar ensemble de top-5 modelos
- [ ] Adicionar early stopping inteligente
- [ ] Implementar cross-validation k-fold
- [ ] Adicionar análise de sensibilidade
- [ ] Gerar curvas de aprendizado
- [ ] Calcular métricas adicionais (AUC, F1, MCC)
- [ ] Paralelizar trials (n_jobs=4)
- [ ] Implementar busca em duas fases


---


## 🔧 Modificações no Código

### Exemplo: Melhorar Sampler TPE

```python

# No arquivo framework_investigativo_completo.py, linha ~3959:

# ANTES:
sampler = TPESampler(seed=42, n_startup_trials=10)

# DEPOIS:
sampler = TPESampler(
    seed=42,
    n_startup_trials=20,        # Mais exploração inicial
    n_ei_candidates=24,         # Mais candidatos para Expected Improvement
    multivariate=True,          # Considerar correlações entre hiperparâmetros
    warn_independent_sampling=True
)

```text

### Exemplo: Melhorar Pruner

```python

# ANTES:
pruner = MedianPruner(n_startup_trials=5, n_warmup_steps=3)

# DEPOIS:
from optuna.pruners import SuccessiveHalvingPruner

pruner = SuccessiveHalvingPruner(
    min_resource=1,
    reduction_factor=3,
    min_early_stopping_rate=0
)

```

---


## 📚 Referências

1. **Optuna Documentation:** <https://optuna.readthedocs.io/>
2. **Bergstra et al. (2011):** "Algorithms for Hyper-Parameter Optimization"
3. **Li et al. (2018):** "Hyperband: A Novel Bandit-Based Approach to Hyperparameter Optimization"
4. **Akiba et al. (2019):** "Optuna: A Next-generation Hyperparameter Optimization Framework"


---


**Conclusão:** O framework já possui uma base sólida de otimização Bayesiana. As melhorias sugeridas podem aumentar a acurácia em 5-15% e reduzir o tempo de busca em até 50%, dependendo da configuração escolhida.

