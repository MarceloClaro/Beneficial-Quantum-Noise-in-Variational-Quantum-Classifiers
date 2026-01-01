# CPU Optimization Guide for Multiframework Execution

## 🎯 Overview

All three frameworks (PennyLane, Qiskit, Cirq) are **fully CPU-compatible** and do not require GPU. This guide optimizes execution for CPU environments to complete Phase 5-9 experimental runs efficiently.

## ✅ Framework Analysis Results

### CPU Compatibility
- ✅ **PennyLane**: Uses `default.mixed` device (CPU-only, no GPU dependencies)
- ✅ **Qiskit**: Uses Aer simulators (`qasm_simulator`, `statevector_simulator`) - CPU-optimized
- ✅ **Cirq**: Uses density matrix simulator - pure NumPy/SciPy (CPU-only)


### No GPU Dependencies Found
- No `torch.cuda`, `tensorflow.gpu`, or CUDA calls detected
- All computations use NumPy/SciPy linear algebra (CPU-optimized)
- Simulators are designed for CPU execution


## ⚙️ CPU Optimization Strategy

### 1. Reduce Computational Load (While Maintaining Scientific Validity)

**Original Config (configs/experiment_unified.yaml):**

```yaml
training:
  epochs: 100          # Too many for initial validation
  batch_size: 32
  
pennylane:
  shots: null          # Analytic (good)
  
qiskit:
  shots: 1024          # High for validation
  
cirq:
  shots: 1024          # High for validation

```text

**Optimized for CPU (Phase 5-7 Validation):**

```yaml
training:
  epochs: 30           # Sufficient for convergence detection
  batch_size: 16       # Reduced memory footprint
  
pennylane:
  shots: null          # Keep analytic (fastest)
  
qiskit:
  shots: 512           # Reduced for speed, still statistically valid
  
cirq:
  shots: 512           # Reduced for speed

```text

### 2. Subset Experimental Grid

#### Full Grid (2,688 configurations):
- 4 datasets × 5 noise models × 6 ansätze × 4 schedules × 3 frameworks × 2 seeds


#### Validation Grid (432 configurations, ~2-3h runtime):
- 4 datasets × 3 noise models × 3 ansätze × 2 schedules × 3 frameworks × 2 seeds


#### Minimal Grid (144 configurations, ~45-60 min runtime):
- 2 datasets × 3 noise models × 3 ansätze × 2 schedules × 3 frameworks × 1 seed


### 3. Parallel Execution (Optional)

```bash

# Run frameworks in parallel using GNU parallel (if available)
parallel -j 3 ::: \
  "python framework_investigativo_completo.py --config configs/experiment_cpu_optimized.yaml" \
  "python executar_framework_qiskit.py --config configs/experiment_cpu_optimized.yaml" \
  "python executar_framework_cirq.py --config configs/experiment_cpu_optimized.yaml"

```text

## 📋 Execution Plan

### Phase 5: PennyLane Baseline (Est. ~40-60 min)

```bash
python framework_investigativo_completo.py \

  --config configs/experiment_cpu_optimized.yaml \
  --output results/pennylane/20251227_001 \
  --mode validation

```text

#### Optimizations:
- Analytic mode (no shots) = ~2x faster
- 30 epochs sufficient for convergence
- Early stopping if validation loss plateaus


### Phase 6: Qiskit Execution (Est. ~60-90 min)

```bash
python executar_framework_qiskit.py \

  --config configs/experiment_cpu_optimized.yaml \
  --output results/qiskit/20251227_001 \
  --mode validation

```text

#### Optimizations:
- 512 shots (vs 1024) = ~2x faster
- Transpilation caching enabled
- Single-thread Aer backend (avoid overhead)


### Phase 7: Cirq Execution (Est. ~50-70 min)

```bash
python executar_framework_cirq.py \

  --config configs/experiment_cpu_optimized.yaml \
  --output results/cirq/20251227_001 \
  --mode validation

```text

#### Optimizations:
- 512 shots = ~2x faster
- Density matrix only when needed (state vector preferred)
- Optimized gate decomposition


### Phase 8: Comparative Analysis (Est. ~10-15 min)

```bash
python generate_comparative_results.py \

  --run-id 20251227_001 \
  --mode validation

```text

### Phase 9: Figure Generation (Est. ~5-10 min)

```bash
python tools/generate_all_figures.py \

  --run-id 20251227_001

```text

## 🚀 Quick Start Commands

### Option 1: Validation Run (Recommended - ~3-4h total)

```bash

# Install dependencies
pip install -r requirements.txt

# Run validation grid (432 configs)
bash scripts/run_multiframework_validation.sh

```text

### Option 2: Minimal Run (Testing - ~1h total)

```bash

# Run minimal grid (144 configs)
bash scripts/run_multiframework_minimal.sh

```text

### Option 3: Full Run (Production - ~8-12h total)

```bash

# Run full grid (2,688 configs) - for final publication
bash scripts/run_multiframework_full.sh

```text

## 📊 Expected Performance

| Configuration | Total Configs | Est. Time (CPU) | Output Size |
|--------------|---------------|-----------------|-------------|
| **Minimal** | 144 | 45-60 min | ~50 MB |
| **Validation** | 432 | 2-3 hours | ~150 MB |
| **Full** | 2,688 | 8-12 hours | ~1 GB |

#### Hardware Requirements:
- **CPU**: 4+ cores recommended (2 cores minimum)
- **RAM**: 8 GB minimum, 16 GB recommended
- **Disk**: 2 GB free space (validation), 5 GB (full)


## 🔧 CPU-Specific Optimizations Applied

### NumPy/SciPy Threading

```python
import os

# Optimize for CPU (avoid oversubscription)
os.environ["OMP_NUM_THREADS"] = "4"
os.environ["OPENBLAS_NUM_THREADS"] = "4"
os.environ["MKL_NUM_THREADS"] = "4"
os.environ["NUMEXPR_NUM_THREADS"] = "4"

```text

### Memory Management

```python

# Clear memory between experiments
import gc
gc.collect()

# Use memory-mapped arrays for large datasets
np.save('temp_data.npy', data, allow_pickle=False)
data_mmap = np.load('temp_data.npy', mmap_mode='r')

```text

### Batch Processing

```python

# Process configurations in batches to avoid memory issues
batch_size = 10
for i in range(0, len(configs), batch_size):
    batch = configs[i:i+batch_size]
    process_batch(batch)
    gc.collect()

```text

## 📈 Scaling Estimates

| # CPUs | Parallel Frameworks | Est. Time (Validation) |
|--------|---------------------|------------------------|
| 2 | Sequential | ~3.5 hours |
| 4 | Parallel (3 frameworks) | ~1.5 hours |
| 8 | Parallel + multiprocess | ~1 hour |

## ⚠️ Known Limitations

### Qiskit Transpilation Overhead
- **Issue**: Qiskit transpiler adds 20-30% overhead vs. PennyLane
- **Mitigation**: Transpilation caching, reduced circuit depth
- **Expected**: Qiskit ~30% slower than PennyLane (acceptable)


### Cirq Density Matrix Memory
- **Issue**: Density matrices require 4× memory vs. state vectors
- **Mitigation**: Use state vectors where possible, gc.collect() after each run
- **Expected**: Cirq RAM usage ~2× PennyLane (acceptable with 8GB+ RAM)


### Shot Noise Variance
- **Issue**: Reduced shots (512 vs 1024) increases variance
- **Mitigation**: Multiple seeds (already in config), statistical corrections
- **Expected**: ±5% variance (within acceptable ±3-5% tolerance)


## ✅ Validation Checklist

Before running experiments:

- [ ] Dependencies installed: `pip list | grep -E "pennylane|qiskit|cirq"`
- [ ] Config file exists: `ls configs/experiment_cpu_optimized.yaml`
- [ ] Output directories created: `ls -d results/{pennylane,qiskit,cirq}`
- [ ] Disk space available: `df -h .` (check >2 GB free)
- [ ] RAM available: `free -h` (check >6 GB free)


## 🐛 Troubleshooting

### Out of Memory

```bash

# Reduce batch size and shots
sed -i 's/batch_size: 16/batch_size: 8/' configs/experiment_cpu_optimized.yaml
sed -i 's/shots: 512/shots: 256/' configs/experiment_cpu_optimized.yaml

```text

### Too Slow

```bash

# Use minimal grid
bash scripts/run_multiframework_minimal.sh

```text

### Import Errors

```bash

# Reinstall dependencies
pip install --upgrade -r requirements.txt

```

## 📝 Summary

**Key Finding**: All frameworks are CPU-compatible and ready for execution. No GPU required.


**Recommendation**: Start with **Validation Run** (432 configs, ~3h) to verify infrastructure before full production run.


**Next Step**: Execute `bash scripts/run_multiframework_validation.sh` when ready.

