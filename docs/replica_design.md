# 4-Replica PFQMC Design Strategy

## Overview

This document describes the strategy for implementing 4 replicas in the PFQMC simulation. The core idea is to enlarge every Green's function matrix from `nDim × nDim` to `4*nDim × 4*nDim` (or equivalently a block-diagonal structure of 4 independent `nDim × nDim` blocks).

## Current Architecture Summary

### Key Data Structures

| Structure | Location | Current Size | Description |
|-----------|----------|--------------|-------------|
| `g` (Green's function) | `PfQMC` class (`inc/pfqmc.h:12`) | `nDim × nDim` | Main Green's function matrix |
| `UDT.U, UDT.T` | `qr_udt.h` | `nDim × nDim` | Stabilization matrices |
| `UDT.D` | `qr_udt.h` | `nDim` | Diagonal scaling vector |
| `B, B_inv` | `SpinlessVOperator` | `nDim × nDim` | Interaction matrices |
| `g0, g0_inv` | `DenseOperator` | `nDim × nDim` | Single-slice Green's functions |
| `mat, mat_inv` | `DenseOperator` | `nDim × nDim` | Dense operator matrices |

### Key Methods Affected

The Green's function flows through the following critical paths:

1. **Initialization**: `rightInit()`, `leftInit()` in `PfQMC`
2. **Sweeps**: `rightSweep()`, `leftSweep()` in `PfQMC`
3. **Updates**: `update()`, `singleFlip()` in `SpinlessVOperator`
4. **Propagation**: `left_propagate()`, `right_propagate()` in operators
5. **Stabilization**: `UDT` decomposition and `onePlusInv()`
6. **Sign calculation**: `getSignRaw()` in `PfQMC`

## Replica Design Strategy

### Option A: Block-Diagonal Matrices (Recommended)

Enlarge all matrices to `4*nDim × 4*nDim` with block-diagonal structure:

```
G_replica = | G_1   0    0    0  |
            |  0   G_2   0    0  |
            |  0    0   G_3   0  |
            |  0    0    0   G_4 |
```

**Advantages:**
- Minimal code changes - existing matrix operations mostly work
- Natural parallelization per replica block
- Memory locality for cache efficiency
- UDT decomposition naturally handles block structure

**Disadvantages:**
- Increased memory footprint (16× for matrices)
- Some operations need block-aware implementations

### Option B: Vector of Matrices

Store 4 separate `nDim × nDim` matrices:
```cpp
std::vector<MatType> g_replicas(4);  // or std::array<MatType, 4>
```

**Advantages:**
- No memory waste from storing zeros
- Clean separation of replicas

**Disadvantages:**
- Requires significant refactoring of all interfaces
- More complex code paths

### Recommendation: Hybrid Approach

Use **Option A (block-diagonal)** for the core `PfQMC` class internals where matrix algebra dominates, but provide **replica-indexed accessors** for measurements and updates that operate per-replica.

## Detailed Implementation Plan

### Phase 1: Core Data Structure Changes

#### 1.1 Type Definitions (`inc/types.h`)

Add replica count constant and helper macros:
```cpp
constexpr int N_REPLICAS = 4;

// Helper to convert replica-local index to global index
inline int replicaIdx(int replica, int localIdx, int nDim) {
    return replica * nDim + localIdx;
}
```

#### 1.2 PfQMC Class (`inc/pfqmc.h`, `src/pfqmc.cpp`)

**Changes needed:**

| Member | Current | New | Notes |
|--------|---------|-----|-------|
| `nDim` | int | int | Keep as single-replica dimension |
| `nDimTotal` | (new) | int | = `4 * nDim` |
| `g` | `MatType(nDim, nDim)` | `MatType(4*nDim, 4*nDim)` | Block-diagonal |
| `udtL`, `udtR` | `vector<UDT>` | `vector<UDT>` | Each UDT enlarged to `4*nDim` |
| `sign` | `DataType` | `std::array<DataType, 4>` | Per-replica signs |

**Methods to update:**

- `PfQMC::PfQMC()` - constructor, line 5-30 in `pfqmc.cpp`
- `PfQMC::rightInit()` - lines 24-49 in `pfqmc.h`
- `PfQMC::leftInit()` - lines 51-77 in `pfqmc.h`
- `PfQMC::rightSweep()` - lines 32-79 in `pfqmc.cpp`
- `PfQMC::leftSweep()` - lines 81-123 in `pfqmc.cpp`
- `PfQMC::getSignRaw()` - lines 125-192 in `pfqmc.cpp`

#### 1.3 UDT Class (`inc/qr_udt.h`)

The UDT decomposition should work transparently on larger matrices. However, verify:

- Constructor `UDT(MatType &A)` - lines 62-99
- `onePlusInv()` - lines 178-191
- `operator*` overloads - lines 194-213
- Global `onePlusInv(UDT&, UDT&)` - lines 216-253

**Potential optimization:** Block-aware UDT that decomposes each replica independently for better numerical stability.

### Phase 2: Operator Class Changes

#### 2.1 Base Operator Interface (`inc/operator.h`)

All virtual methods need to handle `4*nDim` matrices:

| Method | Lines | Change Required |
|--------|-------|-----------------|
| `left_multiply()` | 15 | Auto-works with larger matrices |
| `right_multiply()` | 18 | Auto-works with larger matrices |
| `left_propagate()` | 21 | Auto-works with larger matrices |
| `right_propagate()` | 22 | Auto-works with larger matrices |
| `update()` | 23 | **Needs block-aware implementation** |
| `getGreensMat()` | 27 | Must generate block-diagonal |
| `stabilizedLeftMultiply()` | 31 | Auto-works with larger matrices |

#### 2.2 DenseOperator (`inc/operator.h:38-130`)

**Changes:**
- `mat`, `mat_inv`: Enlarge to `4*nDim × 4*nDim` block-diagonal
- `g0`, `g0_inv`: Enlarge to `4*nDim × 4*nDim` block-diagonal
- Constructor needs to replicate single-replica matrix 4 times

#### 2.3 SpinlessVOperator (`inc/spinless_tV.h:141-405`)

**Critical changes:**

| Member/Method | Lines | Change |
|---------------|-------|--------|
| `B`, `B_inv` | 158-159 | Enlarge to `4*nDim × 4*nDim` |
| `singleFlip()` | 189-303 | Update for each replica independently |
| `singleFlipSingleMajorana()` | 305-344 | Update for each replica independently |
| `update()` | 346-364 | Loop over replicas |
| `getGreensMat()` | 389-392 | Generate block-diagonal |

**Key insight for `singleFlip()`**: The auxiliary field `s` is shared across replicas, but each replica's Green's function block is updated independently. The acceptance ratio `r` can be computed per-replica or as a product across replicas (depending on your MCMC scheme).

### Phase 3: Model-Specific Classes

#### 3.1 SpinlessTvUtils (`inc/spinless_tV.h:8-139`)

Methods that generate matrices need to produce block-diagonal structures:

- `InteractionTanhGenerator()` - lines 50-87
- `InteractionBGenerator()` - lines 90-138

Add helper method:
```cpp
// Replicate a single-replica matrix to block-diagonal form
void replicateToBlockDiagonal(const MatType& single, MatType& block, int nReplicas);
```

#### 3.2 Lattice-Specific Utilities

Update all lattice configuration classes:
- `SpinlessTvHoneycombUtils` (`inc/honeycomb.h`)
- `SpinlessTvSquareUtils` (`inc/square.h`)
- `SpinlessTvChainUtils` (`inc/kitaevChain.h`)
- `SpinlessTvChain1dUtils` (`inc/chain1d_tV.h`)

### Phase 4: Measurement Updates

#### 4.1 Energy and Observables

All measurement functions need to:
1. Extract individual replica blocks from the full Green's function
2. Compute observables per replica
3. Optionally average or correlate across replicas

Example pattern:
```cpp
DataType energyFromGreensFunc(const MatType& g_full) {
    DataType total = 0.0;
    for (int r = 0; r < N_REPLICAS; r++) {
        MatType g_r = g_full.block(r*nDim, r*nDim, nDim, nDim);
        total += computeEnergySingleReplica(g_r);
    }
    return total / N_REPLICAS;
}
```

#### 4.2 Inter-Replica Correlations (New)

Add new measurement functions to compute correlations between replicas:
```cpp
DataType interReplicaCorrelation(const MatType& g_full, int r1, int r2);
```

## Files to Update Summary

### Must Modify

| File | Priority | Estimated Complexity |
|------|----------|---------------------|
| `inc/types.h` | High | Low |
| `inc/pfqmc.h` | High | Medium |
| `src/pfqmc.cpp` | High | High |
| `inc/qr_udt.h` | High | Medium |
| `inc/operator.h` | High | Medium |
| `inc/spinless_tV.h` | High | High |

### May Need Updates

| File | Notes |
|------|-------|
| `inc/honeycomb.h` | Measurement functions |
| `inc/square.h` | Measurement functions |
| `inc/kitaevChain.h` | Measurement functions |
| `inc/chain1d_tV.h` | Measurement functions |
| `main.cpp` | Simulation loops and output |
| `inc/skewMatUtils.h` | Pfaffian computations (if signs are per-replica) |

## Implementation Order

1. **types.h** - Add constants and helpers
2. **qr_udt.h** - Verify UDT works with larger matrices (may need block-aware version)
3. **operator.h** - Update base class interface
4. **spinless_tV.h** - Update SpinlessVOperator (most complex)
5. **pfqmc.h/cpp** - Update PfQMC class
6. **Model files** - Update lattice configs and measurements
7. **main.cpp** - Update simulation drivers

## Testing Strategy

1. **Unit tests**: Verify block-diagonal matrix operations
2. **Regression**: Single replica should reproduce existing results
3. **Symmetry**: All 4 replicas should give statistically identical results for independent observables
4. **Performance**: Benchmark to ensure overhead is acceptable

## Open Questions

1. **Shared vs. independent auxiliary fields**: Should all replicas share the same HS field `s`, or should each have independent fields?
   - Shared: Simpler, correlates replicas
   - Independent: True replica exchange, but 4× more MC updates

2. **Sign computation**: Per-replica signs or combined product?

3. **Replica exchange moves**: Should we implement swap moves between replicas?

4. **Memory optimization**: Consider block-aware storage to avoid storing zeros explicitly.
