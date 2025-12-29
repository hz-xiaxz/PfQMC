# 4-Replica PFQMC Design Strategy

## Overview

This document describes the strategy for implementing 4 replicas in the PFQMC simulation. The key features are:

1. **4 independent replicas** with **distinct auxiliary fields** (not shared)
2. **Block 1-2 mixing** and **Block 3-4 mixing** operators applied **only at τ=0**
3. Green's function enlarged from `nDim × nDim` to `4*nDim × 4*nDim`

### Physical Picture

```
τ = β  ←─────────────────────────────────────────────────── τ = 0

       ┌─────────────────────────────────────────────────┐
Rep 1: │  B_V[s1] ─── B_K ─── B_V[s1] ─── ... ─── B_V[s1]│──┐
       └─────────────────────────────────────────────────┘  │ M_12 mixing
       ┌─────────────────────────────────────────────────┐  │
Rep 2: │  B_V[s2] ─── B_K ─── B_V[s2] ─── ... ─── B_V[s2]│──┘
       └─────────────────────────────────────────────────┘

       ┌─────────────────────────────────────────────────┐
Rep 3: │  B_V[s3] ─── B_K ─── B_V[s3] ─── ... ─── B_V[s3]│──┐
       └─────────────────────────────────────────────────┘  │ M_34 mixing
       ┌─────────────────────────────────────────────────┐  │
Rep 4: │  B_V[s4] ─── B_K ─── B_V[s4] ─── ... ─── B_V[s4]│──┘
       └─────────────────────────────────────────────────┘
```

The mixing operators M_12 and M_34 couple the off-diagonal blocks of the Green's function at τ=0, enabling inter-replica correlations.

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

### Green's Function Block Structure

The 4-replica Green's function has the structure:

```
G = | G_11  G_12   0     0   |
    | G_21  G_22   0     0   |
    |  0     0    G_33  G_34 |
    |  0     0    G_43  G_44 |
```

Where:
- **Diagonal blocks** `G_ii`: Intra-replica Green's functions (always present)
- **Off-diagonal blocks** `G_12, G_21`: Inter-replica correlations between replicas 1-2 (induced by M_12)
- **Off-diagonal blocks** `G_34, G_43`: Inter-replica correlations between replicas 3-4 (induced by M_34)
- **Zero blocks**: No mixing between (1,2) pair and (3,4) pair

This is a **2×2 block structure** at the pair level:
```
G = | G_pair12    0       |
    |    0     G_pair34   |
```

Each pair block is `2*nDim × 2*nDim`.

### Auxiliary Field Structure

Each replica has its **own independent auxiliary field**:

```cpp
// Current: single aux field per operator
iVecType* s;  // size = nBonds

// New: 4 independent aux fields
std::array<iVecType*, 4> s_replica;  // s_replica[r] has size = nBonds
```

During MC updates, each replica's auxiliary field is updated independently, but the acceptance criterion may involve the full determinant/pfaffian ratio.

---

## Mixing Operators at τ=0

### Operator Array Structure

The current `op_array` is a linear sequence of operators applied from τ=0 to τ=β:

```
Current: op_array = [B_0, B_1, B_2, ..., B_{L-1}]
         τ:          0    Δτ   2Δτ       (L-1)Δτ
```

With mixing operators, we need a **special operator at τ=0** that couples replica pairs:

```
New: op_array = [M_mix, B_0, B_1, B_2, ..., B_{L-1}]
                  ↑
            Mixing operator at τ=0
            (applied once per sweep, not at every time slice)
```

### Mixing Operator Definition

The mixing operator `M_mix` has a `4*nDim × 4*nDim` block structure:

```
M_mix = | M_12    0   |
        |   0   M_34  |
```

Where `M_12` and `M_34` are `2*nDim × 2*nDim` matrices that mix replica pairs:

```
M_12 = | A_12  B_12 |      M_34 = | A_34  B_34 |
       | C_12  D_12 |             | C_34  D_34 |
```

Each sub-block is `nDim × nDim`. The specific form of A, B, C, D depends on the physical mixing you want (e.g., swap operators, entanglement operators, etc.).

### New Operator Class: `MixingOperator`

```cpp
// inc/mixing_operator.h (NEW FILE)

class MixingOperator : public Operator {
public:
    int nDim;           // single-replica dimension
    int nDimTotal;      // = 4 * nDim

    // M_12 block (2*nDim × 2*nDim) - mixes replicas 1 and 2
    MatType M_12;
    MatType M_12_inv;

    // M_34 block (2*nDim × 2*nDim) - mixes replicas 3 and 4
    MatType M_34;
    MatType M_34_inv;

    // Full mixing matrix (4*nDim × 4*nDim)
    MatType mat;
    MatType mat_inv;

    MixingOperator(int _nDim, const MatType& _M12, const MatType& _M34);

    // Operator interface
    void left_multiply(const MatType& A, MatType& B) override;
    void right_multiply(const MatType& A, MatType& B) override;
    void left_propagate(MatType& g, MatType& tmp) override;
    void right_propagate(MatType& g, MatType& tmp) override;

    // No auxiliary field updates for mixing operator
    DataType update(MatType& g) override { return 1.0; }

    void getGreensMat(MatType& g0) override;
    void stabilizedLeftMultiply(UDT& F) override;
};
```

### Implementation Details

```cpp
MixingOperator::MixingOperator(int _nDim, const MatType& _M12, const MatType& _M34)
    : nDim(_nDim), nDimTotal(4 * _nDim)
{
    M_12 = _M12;
    M_34 = _M34;
    M_12_inv = _M12.inverse();
    M_34_inv = _M34.inverse();

    // Build full block matrix
    mat = MatType::Zero(nDimTotal, nDimTotal);
    mat.block(0, 0, 2*nDim, 2*nDim) = M_12;
    mat.block(2*nDim, 2*nDim, 2*nDim, 2*nDim) = M_34;

    mat_inv = MatType::Zero(nDimTotal, nDimTotal);
    mat_inv.block(0, 0, 2*nDim, 2*nDim) = M_12_inv;
    mat_inv.block(2*nDim, 2*nDim, 2*nDim, 2*nDim) = M_34_inv;
}

void MixingOperator::left_multiply(const MatType& A, MatType& B) {
    B = mat * A;
}

void MixingOperator::right_multiply(const MatType& A, MatType& B) {
    B = A * mat;
}

void MixingOperator::left_propagate(MatType& g, MatType& tmp) {
    tmp = mat * g;
    g = tmp * mat_inv;
}

void MixingOperator::right_propagate(MatType& g, MatType& tmp) {
    tmp = mat_inv * g;
    g = tmp * mat;
}

void MixingOperator::getGreensMat(MatType& g0) {
    // Green's function for mixing operator: g0 = 2(1+M)^{-1} - 1
    MatType identity = MatType::Identity(nDimTotal, nDimTotal);
    g0 = 2.0 * (identity + mat).inverse() - identity;
}

void MixingOperator::stabilizedLeftMultiply(UDT& F) {
    F = mat * F;
}
```

---

## Modified Operator Array for 4 Replicas

### SpinlessVOperator Changes

The interaction operator now handles 4 replicas with distinct auxiliary fields:

```cpp
class SpinlessVOperator4Replica : public Operator {
protected:
    const SpinlessTvUtils* config;

public:
    const int nDim;         // single-replica dimension
    const int nDimTotal;    // = 4 * nDim
    const int bondType;

    // 4 independent auxiliary field arrays
    std::array<iVecType*, 4> s;  // s[r] for replica r

    // 4 independent B matrices (each nDim × nDim)
    std::array<MatType, 4> B_replica;
    std::array<MatType, 4> B_inv_replica;

    // Full B matrix (4*nDim × 4*nDim, block-diagonal)
    MatType B;
    MatType B_inv;

    rdGenerator* rd;

    SpinlessVOperator4Replica(const SpinlessTvUtils* _config,
                               std::array<iVecType*, 4> _s,
                               int _bondType,
                               rdGenerator* _rd);

    void rebuildFullB();  // Rebuild B from B_replica blocks

    // Update each replica independently
    DataType update(MatType& g) override;

    // Single flip for specific replica
    void singleFlip(MatType& g, int replica, int idxAux, double rand,
                    bool& flag, DataType& signCur);
};
```

### Update Method for 4 Replicas

```cpp
DataType SpinlessVOperator4Replica::update(MatType& g) {
    DataType signCur = 1.0;
    bool flag;
    double rand;

    // Update each replica independently
    for (int r = 0; r < 4; r++) {
        for (int i = 0; i < s[r]->size(); i++) {
            rand = rd->rdUniform01();
            singleFlip(g, r, i, rand, flag, signCur);
        }
    }
    return signCur;
}

void SpinlessVOperator4Replica::singleFlip(MatType& g, int replica,
                                            int idxAux, double rand,
                                            bool& flag, DataType& signCur) {
    // Extract the relevant 2×2 block for this replica pair
    // Replica 0,1 share pair block [0:2*nDim, 0:2*nDim]
    // Replica 2,3 share pair block [2*nDim:4*nDim, 2*nDim:4*nDim]

    int pairOffset = (replica < 2) ? 0 : 2*nDim;
    int localReplica = replica % 2;
    int replicaOffset = localReplica * nDim;
    int globalOffset = pairOffset + replicaOffset;

    // Get indices in the global matrix
    int idx1_global, idx2_global;
    // ... compute based on aux field index and replica ...

    // The acceptance ratio now involves the 2×2 pair block
    // because replicas within a pair are coupled by mixing operator

    // ... rest of update logic using pair-block of g ...
}
```

---

## Imaginary Time Flow with Mixing

### Sweep Structure

The sweep now accounts for the mixing operator at τ=0:

```
rightSweep():
  1. Start at τ=0 with G at τ=0
  2. Apply mixing operator M_mix (propagate G through mixing)
  3. For l = 0 to L-1:
     a. Update auxiliary fields for all 4 replicas at time slice l
     b. Propagate G: G → B_l * G * B_l^{-1}
     c. Stabilization at checkpoints

leftSweep():
  1. Start at τ=β with G at τ=β
  2. For l = L-1 down to 0:
     a. Propagate G: G → B_l^{-1} * G * B_l
     b. Update auxiliary fields for all 4 replicas
     c. Stabilization at checkpoints
  3. Apply inverse mixing operator M_mix^{-1}
```

### Modified PfQMC Class

```cpp
class PfQMC4Replica {
public:
    int stb;
    int nDim;           // single-replica dimension
    int nDimTotal;      // = 4 * nDim

    MatType g;          // 4*nDim × 4*nDim Green's function

    std::vector<Operator*> op_array;  // Regular operators (B_V, B_K, etc.)
    MixingOperator* mixingOp;          // Special mixing operator at τ=0

    int op_length;
    std::vector<bool> need_stabilization;
    int checkpoints;
    std::vector<UDT> udtL;
    std::vector<UDT> udtR;

    DataType sign;  // Combined sign (or array of per-replica signs)

    PfQMC4Replica(Spinless_tV_4Replica* walker,
                  const MatType& M12, const MatType& M34,
                  int _stb = 10);

    void rightInit();
    void leftInit();
    void rightSweep();
    void leftSweep();
    DataType getSignRaw();
};
```

### rightSweep with Mixing

```cpp
void PfQMC4Replica::rightSweep() {
    MatType tmp = MatType::Identity(nDimTotal, nDimTotal);
    MatType Aseg = MatType::Identity(nDimTotal, nDimTotal);
    int curSeg = 0;
    DataType signCur;

    // Apply mixing operator at τ=0 first
    mixingOp->left_propagate(g, tmp);
    mixingOp->left_multiply(Aseg, tmp);
    std::swap(Aseg, tmp);

    // Then sweep through regular time slices
    for (int l = 0; l < op_length; l++) {
        signCur = op_array[l]->update(g);
        this->sign *= signCur;

        op_array[l]->left_multiply(Aseg, tmp);
        std::swap(Aseg, tmp);

        if (need_stabilization[(l + 1) % op_length]) {
            // ... stabilization code (same structure as before) ...
        } else {
            op_array[l]->left_propagate(g, tmp);
        }
    }
}
```

---

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

### New Files to Create

| File | Description |
|------|-------------|
| `inc/mixing_operator.h` | MixingOperator class for τ=0 replica pair mixing |
| `inc/spinless_tV_4replica.h` | 4-replica version of SpinlessVOperator |
| `inc/pfqmc_4replica.h` | 4-replica version of PfQMC class |

### Must Modify

| File | Priority | Changes |
|------|----------|---------|
| `inc/types.h` | High | Add `N_REPLICAS=4`, helper macros |
| `inc/qr_udt.h` | High | Verify works with `4*nDim` matrices |
| `inc/operator.h` | Medium | Base class unchanged, but verify interface |

### May Need Updates (for measurements)

| File | Notes |
|------|-------|
| `inc/honeycomb.h` | Add 4-replica measurement functions |
| `inc/square.h` | Add 4-replica measurement functions |
| `inc/kitaevChain.h` | Add 4-replica measurement functions |
| `inc/chain1d_tV.h` | Add 4-replica measurement functions |
| `main.cpp` | Add 4-replica simulation entry points |
| `inc/skewMatUtils.h` | May need pair-block pfaffian computations |

---

## Detailed Function Changes

### Functions in `src/pfqmc.cpp` to Update/Create

| Function | Line | Change Type | Description |
|----------|------|-------------|-------------|
| `PfQMC::PfQMC()` | 5-30 | **New version** | Create `PfQMC4Replica` constructor handling mixing operator |
| `PfQMC::rightInit()` | 24-49 (in .h) | **New version** | Account for mixing operator in UDT chain |
| `PfQMC::leftInit()` | 51-77 (in .h) | **New version** | Account for mixing operator in UDT chain |
| `PfQMC::rightSweep()` | 32-79 | **New version** | Apply mixing at τ=0, update 4 replicas |
| `PfQMC::leftSweep()` | 81-123 | **New version** | Apply inverse mixing, update 4 replicas |
| `PfQMC::getSignRaw()` | 125-192 | **New version** | Compute sign with pair-block structure |

### Functions in `inc/spinless_tV.h` to Update/Create

| Function | Line | Change Type | Description |
|----------|------|-------------|-------------|
| `SpinlessVOperator` constructor | 165-177 | **New class** | `SpinlessVOperator4Replica` with 4 aux fields |
| `singleFlip()` | 189-303 | **New version** | Per-replica flip with pair-block acceptance |
| `update()` | 346-364 | **New version** | Loop over 4 replicas |
| `rebuildFullB()` | (new) | **New function** | Build 4*nDim block-diagonal B from per-replica B |
| `left_multiply()` | 370-372 | **Inherit** | Works with larger matrix |
| `right_multiply()` | 366-368 | **Inherit** | Works with larger matrix |
| `left_propagate()` | 374-378 | **Modify** | Use rebuilt full B and B_inv |
| `right_propagate()` | 379-383 | **Modify** | Use rebuilt full B and B_inv |
| `getGreensMat()` | 389-392 | **New version** | Return 4*nDim block-diagonal g0 |

### Functions in `inc/qr_udt.h` to Verify

| Function | Line | Status | Notes |
|----------|------|--------|-------|
| `UDT(MatType& A)` | 62-99 | **Verify** | Should work with 4*nDim matrices |
| `onePlusInv()` | 178-191 | **Verify** | May need block-aware version for stability |
| `operator*(UDT, UDT)` | 194-203 | **Verify** | Should work transparently |
| `operator*(MatType, UDT)` | 206-213 | **Verify** | Should work transparently |
| `onePlusInv(UDT&, UDT&)` | 216-253 | **Verify** | May need block-aware version |

---

## Key Implementation Insights

### 1. Pair-Block Structure is Essential

Since mixing only occurs within pairs (1-2) and (3-4), many operations can be optimized:

```cpp
// Instead of full 4*nDim × 4*nDim operations:
MatType g_pair12 = g.block(0, 0, 2*nDim, 2*nDim);
MatType g_pair34 = g.block(2*nDim, 2*nDim, 2*nDim, 2*nDim);

// Operations on pairs can be done independently (parallelizable)
```

### 2. Acceptance Ratio for Coupled Replicas

When replicas 1 and 2 are coupled by M_12, flipping an auxiliary field in replica 1 affects the acceptance ratio through the off-diagonal blocks G_12 and G_21:

```cpp
// Simplified acceptance ratio for coupled pair
DataType computePairAcceptanceRatio(const MatType& g_pair, int replica,
                                     int idx1, int idx2, DataType thlV, int auxCur) {
    // g_pair is 2*nDim × 2*nDim
    int offset = replica * nDim;

    // Diagonal contribution (same as single replica)
    DataType r_diag = 1.0 - 1.0i * thlV * auxCur * g_pair(offset+idx1, offset+idx2);

    // Off-diagonal contribution from inter-replica correlations
    int other_offset = (1 - replica) * nDim;
    DataType r_offdiag = /* contribution from G_12/G_21 blocks */;

    return r_diag * r_offdiag * etaM;
}
```

### 3. Sign Computation

The sign computation in `getSignRaw()` needs to account for:
- Product of pfaffians/determinants from each replica
- Contribution from mixing operator
- Proper handling of pair-block structure

---

## Testing Strategy

### Unit Tests
1. **MixingOperator**: Verify `M * M_inv = I`, propagation correctness
2. **Block structure**: Verify zero blocks remain zero after sweeps
3. **Single replica limit**: When M_12 = M_34 = I (identity), should reduce to 4 independent replicas

### Integration Tests
1. **Sign consistency**: `getSignRaw()` should match accumulated sign
2. **Energy**: Compare against ED for small systems
3. **Replica symmetry**: Statistical equivalence of replica pairs

### Regression Tests
1. **No mixing case**: Verify results match 4× single-replica runs

