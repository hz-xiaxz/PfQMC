# Replica Implementation Design for PfQMC

## Overview

This document describes the design for implementing replica support in the Pfaffian QMC codebase. The goal is to simulate 4 identical copies (replicas) of a t-V model, resulting in a block-diagonal `4N × 4N` Green's function matrix instead of a single `N × N` block.

## Architecture Analysis

### Current Class Hierarchy

```
SpinlessTvUtils (base configuration)
    ↓
SpinlessTvChainUtils (chain-specific config)
    ↓
Chain_tV (walker that builds operator sequence)
    ↓
PfQMC (algorithm driver)
```

**Key Observations:**
1. `PfQMC` is **dimension-agnostic** - it only cares about `nDim` from the walker
2. Operators (`DenseOperator`, `SpinlessVOperator`) work with matrices of any size
3. Matrix generation happens in Utils classes (`KineticGenerator()`, `InteractionBGenerator()`)
4. The `nDim` parameter propagates through: `Utils → Walker → PfQMC → Operators`

## Design Decision: Direct Parameter Approach

**Selected Approach:** Add replica parameter directly to existing classes (no separate wrapper)

### Pros
- ✅ Minimal code changes (~6 functions modified)
- ✅ Leverages existing dimension-agnostic architecture
- ✅ Easy to maintain and understand
- ✅ Can be enabled/disabled with single parameter (`nReplica=1` for single model)
- ✅ Works with all existing observables and analysis tools

### Cons
- ❌ Mixes single-replica and multi-replica logic in same class
- ❌ Need to update each model separately (but models already differ)

### Alternative Approaches Considered

**Option A: Wrapper Class**
```cpp
class ReplicaUtils {
    int nReplica;
    SpinlessTvUtils* singleModel;
    // Wrap all matrix generators
};
```
- ❌ Extra abstraction layer
- ❌ Need to forward all method calls
- ❌ Complicates code structure

**Option B: Inheritance**
```cpp
class ReplicaChainUtils : public SpinlessTvChainUtils {
    // Override matrix generators
};
```
- ❌ Need separate replica class for each model
- ❌ Code duplication
- ✅ Clean separation of concerns

## Implementation Plan

### 1. Modify Utils Classes

Add `nReplica` parameter and update dimension calculation:

**File:** `inc/kitaevChain.h:10` (SpinlessTvChainUtils)

```cpp
class SpinlessTvChainUtils : public SpinlessTvUtils {
public:
    int nReplica;  // NEW: number of replicas (default 1)

    SpinlessTvChainUtils(int Lx_, double dt_, double V_, int LTau_,
                         int boundaryType_, double delta_, double mu_,
                         int hsScheme_, int nReplica_ = 1)  // NEW PARAMETER
        : SpinlessTvUtils(Lx_, 1, dt_, V_, LTau_, 2*Lx_*nReplica_, false, hsScheme_),
          nReplica(nReplica_) {  // nDim now accounts for replicas
        // ... rest of initialization
    }
};
```

**Key Change:** `nDim = 2 * Lx * nReplica` instead of `nDim = 2 * Lx`

### 2. Modify Matrix Generators

Create block-diagonal structure for kinetic and interaction terms.

#### A. Kinetic Term (inc/kitaevChain.h:47)

**Current:** Generates `H` of size `nDim × nDim` (with `nDim = 2*Lx`)

**Modified:**
```cpp
inline void KineticGenerator(MatType &H) const {
    H.setZero();

    if (nReplica == 1) {
        // Original single-replica code
        generateKineticSingleReplica(H);
    } else {
        // Build single-replica kinetic matrix
        int nDim_single = 2 * Lx;
        MatType H_single = MatType::Zero(nDim_single, nDim_single);
        generateKineticSingleReplica(H_single);

        // Replicate to block diagonal
        for (int r = 0; r < nReplica; r++) {
            H.block(r*nDim_single, r*nDim_single, nDim_single, nDim_single) = H_single;
        }
    }
}

// Helper function (extract existing logic)
void generateKineticSingleReplica(MatType &H) const {
    // All existing hopping, pairing, chemical potential code goes here
    int idx1, idx2;
    DataType tmp = (1.0i);
    // ... lines 50-125 of current implementation
}
```

**Matrix Structure:**
```
H = [ H_single    0         0         0      ]
    [ 0           H_single  0         0      ]
    [ 0           0         H_single  0      ]
    [ 0           0         0         H_single]
```

#### B. Interaction Term (inc/spinless_tV.h:91)

**Current:** Generates `B` of size `nDim × nDim` based on auxiliary fields `s`

**Critical Design Decision:** Are auxiliary fields shared or independent?

##### Option 1: Shared Auxiliary Fields (Recommended)
All replicas use the **same** auxiliary field configuration `s`. This maintains replica symmetry and is physically motivated for studying replica tricks.

```cpp
inline void InteractionBGenerator(MatType &B, const iVecType &s,
                                  const int bondType, bool inv = false) const {
    B.setIdentity();

    if (nReplica == 1) {
        generateInteractionSingleReplica(B, s, bondType, inv);
    } else {
        int nDim_single = 2 * Lx;  // For spinless models
        MatType B_single = MatType::Identity(nDim_single, nDim_single);
        generateInteractionSingleReplica(B_single, s, bondType, inv);

        // Replicate to all blocks
        for (int r = 0; r < nReplica; r++) {
            B.block(r*nDim_single, r*nDim_single, nDim_single, nDim_single) = B_single;
        }
    }
}
```

**Auxiliary field size:** `s.size() = nBond` (unchanged - same fields for all replicas)

##### Option 2: Independent Auxiliary Fields (Alternative)
Each replica has **independent** auxiliary fields. Requires 4× auxiliary fields and modified indexing.

```cpp
// Not recommended - breaks replica symmetry
// Would need: s.size() = nBond * nReplica
// And aux2MajoranaIdx would need replica indexing
```

**Recommendation:** Use **Option 1 (shared fields)** for physical relevance and simplicity.

### 3. Update Auxiliary Field Indexing

**File:** `inc/kitaevChain.h:25` (aux2MajoranaIdx)

With **shared auxiliary fields**, the current implementation works unchanged because:
- `auxIdx` maps to bond within a single replica
- Same mapping applies to all replicas due to block-diagonal structure
- Each block gets identical interaction matrices

**No changes needed** if using shared fields.

If using independent fields (not recommended):
```cpp
void aux2MajoranaIdx(int auxIdx, int k, int bondType, int &idx1, int &idx2) const {
    // Extract replica and local bond index
    int nBond_single = (boundaryType == 0) ? Lx : (Lx - 1);
    int replica_idx = auxIdx / nBond_single;
    int local_auxIdx = auxIdx % nBond_single;

    // Get indices within single replica
    int local_idx1, local_idx2;
    getSingleReplicaIndices(local_auxIdx, k, bondType, local_idx1, local_idx2);

    // Offset by replica
    int nDim_single = 2 * Lx;
    idx1 = replica_idx * nDim_single + local_idx1;
    idx2 = replica_idx * nDim_single + local_idx2;
}
```

### 4. Update Observable Calculations

Observables need to handle block-diagonal Green's function.

#### Option A: Average Over Replicas (Recommended)

```cpp
DataType energyFromGreensFunc(MatType &g) const {
    if (nReplica == 1) {
        return energyFromGreensFuncSingleReplica(g);
    }

    DataType energy_total = 0.0;
    int nDim_single = 2 * Lx;

    for (int r = 0; r < nReplica; r++) {
        MatType g_replica = g.block(r*nDim_single, r*nDim_single,
                                     nDim_single, nDim_single);
        energy_total += energyFromGreensFuncSingleReplica(g_replica);
    }

    return energy_total / nReplica;
}
```

Apply to all observables:
- `energyFromGreensFunc()` (kitaevChain.h:309)
- `particleNumber()` (kitaevChain.h:228 or similar)
- `EdgeCorrelator()` (kitaevChain.h:129)
- `StructureFactorCDW()` (kitaevChain.h:136)
- `Z2FermionParity()` (kitaevChain.h:238)

#### Option B: Per-Replica Observables

For debugging or studying replica differences:
```cpp
std::vector<DataType> energyPerReplica(MatType &g) const {
    std::vector<DataType> energies(nReplica);
    int nDim_single = 2 * Lx;

    for (int r = 0; r < nReplica; r++) {
        MatType g_replica = g.block(r*nDim_single, r*nDim_single,
                                     nDim_single, nDim_single);
        energies[r] = energyFromGreensFuncSingleReplica(g_replica);
    }
    return energies;
}
```

### 5. Verify Sign Calculation

**File:** `src/pfqmc.cpp:125` (getSignRaw)

**Key Question:** Does sign calculation need modification?

**Analysis:**
- Sign is computed from Pfaffians of 4N×4N blocks (pfqmc.cpp:125-192)
- For block-diagonal structure: `Pf(block_diag(A, A, A, A)) = Pf(A)^4`
- Current implementation builds full 4N×4N matrices in `pfaffianForSignOfProduct()`

**Potential Optimization:**
```cpp
DataType PfQMC::getSignRaw() {
    if (replica_mode && nReplica == 4) {
        // Compute sign from single replica block
        DataType sign_single = computeSignSingleReplica();
        // Pfaffian of block-diagonal is product: Pf^4
        // For sign: sign^4 = ±1 (always +1 for even power)
        return 1.0;  // or more carefully handle sign_single^4
    }
    // ... existing implementation
}
```

**Recommendation:** Leave unchanged initially - current implementation should work correctly, just slower. Profile later if performance is an issue.

### 6. Update Main Function

**File:** `main.cpp` (chain runner around line 214)

Add `nReplica` parameter:

```cpp
void runChain(int Lx, int LTau, double dt, double V, double delta,
              int nthreads, int nseed, int evaluationLength,
              const std::string &filepath, double mu, int hsScheme,
              int boundary, int nReplica = 1) {  // NEW PARAMETER

    // Create configuration with replica support
    SpinlessTvChainUtils config(Lx, dt, V, LTau, boundary, delta, mu, hsScheme, nReplica);

    // Rest unchanged - PfQMC automatically handles larger dimension
    rdGenerator rd(nseed);
    Chain_tV walker(&config, &rd);
    PfQMC pfqmc(&walker, stabilizationTime);

    // ... thermalization and measurement loops
}
```

Command line:
```bash
# Single replica (original)
mpirun ./build/main --chain 12 400 0.1 0.5 0.0 4 42 1000 output.txt 0.0 0 0

# Four replicas (new)
mpirun ./build/main --chain 12 400 0.1 0.5 0.0 4 42 1000 output.txt 0.0 0 0 4
```

## Summary of Required Changes

| File | Function/Location | Line | Change |
|------|------------------|------|--------|
| `inc/kitaevChain.h` | Constructor | ~15 | Add `nReplica` parameter, update `nDim = 2*Lx*nReplica` |
| `inc/kitaevChain.h` | `KineticGenerator()` | 47 | Create block-diagonal kinetic matrix |
| `inc/spinless_tV.h` | `InteractionBGenerator()` | 91 | Create block-diagonal interaction matrix (shared aux fields) |
| `inc/kitaevChain.h` | `energyFromGreensFunc()` | 309 | Average energy over replicas |
| `inc/kitaevChain.h` | Other observables | Various | Average over replica blocks |
| `main.cpp` | Chain runner | ~214 | Add `nReplica` parameter |

**Total:** ~6 function modifications + 1 parameter addition

## Implementation Checklist

- [ ] Add `nReplica` parameter to `SpinlessTvChainUtils` constructor
- [ ] Update `nDim` calculation to account for replicas
- [ ] Refactor `KineticGenerator()` to support block-diagonal structure
- [ ] Refactor `InteractionBGenerator()` with shared auxiliary fields
- [ ] Update `energyFromGreensFunc()` to average over replicas
- [ ] Update other observable functions (correlators, structure factors, etc.)
- [ ] Add `nReplica` parameter to main function and CLI parsing
- [ ] Test with `nReplica=1` to verify backward compatibility
- [ ] Test with `nReplica=4` for correctness
- [ ] Verify sign calculation remains stable
- [ ] Profile performance and optimize if needed
- [ ] Update documentation and examples

## Testing Strategy

### 1. Backward Compatibility
```bash
# Should give identical results to original code
./build/main --chain 8 200 0.1 0.5 0.0 1 42 1000 test_single.txt 0.0 0 0 1
# vs
./build/main --chain 8 200 0.1 0.5 0.0 1 42 1000 test_original.txt 0.0 0 0
```

### 2. Replica Consistency
```bash
# Energy per replica should match single-replica simulation
./build/main --chain 8 200 0.1 0.5 0.0 1 42 1000 test_replica.txt 0.0 0 0 4
```
Verify: `E_avg(4 replicas) ≈ E(1 replica)` (within statistical error)

### 3. Sign Problem
Monitor `sign` vs `signRaw` - should remain ±1 for sign-problem-free cases.

### 4. Green's Function Structure
Add diagnostic to verify block-diagonal structure:
```cpp
// After thermalization, check off-diagonal blocks are zero
for (int r1 = 0; r1 < nReplica; r1++) {
    for (int r2 = r1+1; r2 < nReplica; r2++) {
        double offdiag = g.block(r1*nDim_single, r2*nDim_single,
                                  nDim_single, nDim_single).norm();
        assert(offdiag < 1e-10);  // Should be numerically zero
    }
}
```

## Physics Motivation

Why simulate 4 identical replicas?

1. **Replica trick in quantum systems**: Study entanglement, disorder averaging
2. **Enhanced statistics**: 4× measurements per sweep (if treating independently)
3. **Symmetry probes**: Test replica symmetry breaking phenomena
4. **Methodological**: Template for replica-based algorithms (e.g., parallel tempering between replicas)

With **shared auxiliary fields**, all replicas evolve identically (deterministically), so observables are perfectly correlated. This is useful for:
- Testing numerical stability of block-diagonal implementation
- Template for future replica-based methods with coupling between replicas
- Studying entanglement properties

## Future Extensions

1. **Inter-replica coupling**: Add terms coupling different replicas (breaks block-diagonal structure)
2. **Replica-specific parameters**: Different `μ`, `V`, `delta` per replica
3. **Replica exchange**: Swap configurations between replicas (parallel tempering)
4. **Independent auxiliary fields**: Break replica symmetry for disorder studies

## References

- Original PfQMC paper: [arXiv:2408.10311](https://arxiv.org/abs/2408.10311)
- Replica trick: Statistical mechanics and quantum field theory textbooks
- Block-diagonal matrices: Linear algebra references
