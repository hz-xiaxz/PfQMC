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

## Implementation Approach: Direct Modification (Recommended)

### Why Direct Modification?

Instead of creating new classes (`SpinlessVOperator4Replica`, `PfQMC4Replica`), we **modify existing classes** with a `nReplicas` parameter that defaults to 1. This approach:

1. **Minimal code duplication** - Single code path to maintain
2. **Easy regression testing** - `nReplicas=1` must match original behavior exactly
3. **Incremental development** - Can test each change against known results
4. **Simpler op_array** - No need to choose between class variants

### Design Pattern

Add `nReplicas` parameter (default=1) to existing classes:

```cpp
class SpinlessVOperator : public Operator {
protected:
    const SpinlessTvUtils* config;
    const double etaM;

public:
    const int nReplicas;      // NEW: default = 1 for backward compatibility
    const int nDimSingle;     // NEW: single-replica dimension (was nDim)
    const int nDim;           // CHANGED: = nReplicas * nDimSingle
    const int bondType;

    // CHANGED: from single pointer to vector
    std::vector<iVecType*> s;           // s[r] for replica r, size = nReplicas

    // NEW: per-replica B matrices
    std::vector<MatType> B_replica;     // size = nReplicas, each nDimSingle × nDimSingle

    // EXISTING: full block-diagonal B (nDim × nDim)
    MatType B;
    MatType B_inv;

    rdGenerator* rd;

    // CHANGED: constructor signature
    SpinlessVOperator(const SpinlessTvUtils* _config,
                      std::vector<iVecType*> _s,    // vector instead of single pointer
                      int _bondType,
                      rdGenerator* _rd,
                      int _nReplicas = 1);          // NEW: default = 1
    // ...
};
```

---

## Files to Update Summary

### New Files to Create

| File | Description |
|------|-------------|
| `inc/mixing_operator.h` | MixingOperator class for τ=0 replica pair mixing (only new class needed) |

### Must Modify (Direct Changes)

| File | Priority | Changes |
|------|----------|---------|
| `inc/types.h` | High | Add helper macros for replica indexing |
| `inc/spinless_tV.h` | High | Add `nReplicas` param, vector of aux fields, per-replica B |
| `inc/pfqmc.h` | High | Add `nReplicas` param, handle mixing operator |
| `src/pfqmc.cpp` | High | Update sweeps for replicas and mixing |
| `inc/operator.h` | Medium | Add `nReplicas` to DenseOperator |
| `inc/qr_udt.h` | Low | Verify works with larger matrices |

### May Need Updates (for measurements)

| File | Notes |
|------|-------|
| `inc/honeycomb.h` | Measurement functions (extract per-replica blocks) |
| `inc/square.h` | Measurement functions |
| `inc/kitaevChain.h` | Measurement functions |
| `inc/chain1d_tV.h` | Measurement functions |
| `main.cpp` | Add 4-replica simulation entry points |

---

## Detailed Function Changes

### Functions in `inc/spinless_tV.h` to Modify

| Function | Line | Change Type | Description |
|----------|------|-------------|-------------|
| `SpinlessVOperator()` | 165-177 | **Modify** | Add `nReplicas` param, init vector of aux fields |
| `singleFlip()` | 189-303 | **Modify** | Add `replica` param, compute global indices |
| `singleFlipSingleMajorana()` | 305-344 | **Modify** | Add `replica` param |
| `update()` | 346-364 | **Modify** | Loop over `nReplicas` |
| `rebuildFullB()` | (new) | **Add** | Build block-diagonal B from B_replica[] |
| `left_propagate()` | 374-378 | **Modify** | Call `rebuildFullB()` first if nReplicas > 1 |
| `right_propagate()` | 379-383 | **Modify** | Call `rebuildFullB()` first if nReplicas > 1 |
| `getGreensMat()` | 389-392 | **Modify** | Build block-diagonal g0 |
| Destructor | 179 | **Modify** | Delete all s[r] pointers |

### Functions in `src/pfqmc.cpp` to Modify

| Function | Line | Change Type | Description |
|----------|------|-------------|-------------|
| `PfQMC::PfQMC()` | 5-30 | **Modify** | Add `nReplicas`, `mixingOp` params |
| `PfQMC::rightInit()` | 24-49 (in .h) | **Modify** | Include mixing operator in UDT chain |
| `PfQMC::leftInit()` | 51-77 (in .h) | **Modify** | Include mixing operator in UDT chain |
| `PfQMC::rightSweep()` | 32-79 | **Modify** | Apply mixing at τ=0 |
| `PfQMC::leftSweep()` | 81-123 | **Modify** | Apply inverse mixing at τ=0 |
| `PfQMC::getSignRaw()` | 125-192 | **Modify** | Handle pair-block structure |

### Functions in `inc/qr_udt.h` to Verify

| Function | Line | Status | Notes |
|----------|------|--------|-------|
| `UDT(MatType& A)` | 62-99 | **Verify** | Should work with larger matrices |
| `onePlusInv()` | 178-191 | **Verify** | Test numerical stability |
| `operator*(UDT, UDT)` | 194-203 | **Verify** | Should work transparently |
| `operator*(MatType, UDT)` | 206-213 | **Verify** | Should work transparently |
| `onePlusInv(UDT&, UDT&)` | 216-253 | **Verify** | Test numerical stability |

---

## Code Sketches for Direct Modification

### Modified SpinlessVOperator Constructor

```cpp
// inc/spinless_tV.h - Modified constructor

SpinlessVOperator(const SpinlessTvUtils* _config,
                  std::vector<iVecType*> _s,
                  int _bondType,
                  rdGenerator* _rd,
                  int _nReplicas = 1)
    : config(_config),
      etaM(_config->etaM),
      nReplicas(_nReplicas),
      nDimSingle(_config->nDim),
      nDim(_nReplicas * _config->nDim),  // total dimension
      bondType(_bondType),
      singleMaj(_config->singleMaj),
      hsScheme(_config->hsScheme)
{
    assert(_s.size() == _nReplicas);
    s = _s;
    rd = _rd;

    // Initialize per-replica B matrices
    B_replica.resize(nReplicas);
    for (int r = 0; r < nReplicas; r++) {
        B_replica[r] = MatType::Identity(nDimSingle, nDimSingle);
        config->InteractionBGenerator(B_replica[r], *s[r], bondType, false);
    }

    // Build full block-diagonal B
    rebuildFullB();
}

// Backward-compatible constructor (single replica)
SpinlessVOperator(const SpinlessTvUtils* _config,
                  iVecType* _s,
                  int _bondType,
                  rdGenerator* _rd)
    : SpinlessVOperator(_config, std::vector<iVecType*>{_s}, _bondType, _rd, 1)
{}
```

### rebuildFullB() Implementation

```cpp
void SpinlessVOperator::rebuildFullB() {
    if (nReplicas == 1) {
        // Single replica: B_replica[0] is the full B
        B = B_replica[0];
        return;
    }

    // Multiple replicas: build block-diagonal
    B = MatType::Zero(nDim, nDim);
    for (int r = 0; r < nReplicas; r++) {
        B.block(r * nDimSingle, r * nDimSingle, nDimSingle, nDimSingle) = B_replica[r];
    }
}

void SpinlessVOperator::rebuildFullBInv() {
    if (nReplicas == 1) {
        B_inv = MatType::Identity(nDimSingle, nDimSingle);
        config->InteractionBGenerator(B_inv, *s[0], bondType, true);
        return;
    }

    B_inv = MatType::Zero(nDim, nDim);
    for (int r = 0; r < nReplicas; r++) {
        MatType B_inv_r = MatType::Identity(nDimSingle, nDimSingle);
        config->InteractionBGenerator(B_inv_r, *s[r], bondType, true);
        B_inv.block(r * nDimSingle, r * nDimSingle, nDimSingle, nDimSingle) = B_inv_r;
    }
}
```

### Modified singleFlip() with Replica Parameter

```cpp
void SpinlessVOperator::singleFlip(MatType& g, int replica, int idxAux,
                                    double rand, bool& flag, DataType& signCur) {
    // Compute offset for this replica in the global matrix
    int offset = replica * nDimSingle;

    DataType r;
    int auxCur = (*s[replica])(idxAux);
    int idx1, idx2, idx3, idx4;
    DataType tmp[2];
    const int inc = 1;
    DataType alpha;

    // Get local indices (within single replica)
    config->aux2MajoranaIdx(idxAux, 0, bondType, idx1, idx2);
    config->aux2MajoranaIdx(idxAux, 1, bondType, idx3, idx4);

    // Convert to global indices
    int idx1_g = offset + idx1;
    int idx2_g = offset + idx2;
    int idx3_g = offset + idx3;
    int idx4_g = offset + idx4;

    if (hsScheme == 0) {
        tmp[0] = (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx1_g, idx2_g)));
        tmp[1] = (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx3_g, idx4_g)));
        r = tmp[0] * tmp[1];
        r += (config->thlV * config->thlV) *
             ((g(idx1_g, idx3_g) * g(idx2_g, idx4_g)) - (g(idx2_g, idx3_g) * g(idx1_g, idx4_g)));
        r *= etaM;
    } else if (hsScheme == 1) {
        // ... similar with global indices ...
    }

    flag = rand < std::abs(r);

    if (flag) {
        signCur *= (r / std::abs(r));

        if (hsScheme == 0) {
            for (int imaj = 0; imaj < 2; imaj++) {
                config->aux2MajoranaIdx(idxAux, imaj, bondType, idx1, idx2);
                int idx1_g = offset + idx1;
                int idx2_g = offset + idx2;

                if (imaj == 1) {
                    tmp[1] = (1.0 - ((1.0i) * (config->thlV) * double(auxCur) * g(idx1_g, idx2_g)));
                }

                // Update aux field and B_replica matrix
                (*s[replica])(idxAux) = -auxCur;
                B_replica[replica](idx1, idx2) = -B_replica[replica](idx1, idx2);
                B_replica[replica](idx2, idx1) = -B_replica[replica](idx2, idx1);

                // Update Green's function (using global indices)
                cVecType x1 = -g.col(idx1_g);
                cVecType x2 = -g.col(idx2_g);
                x1(idx1_g) += 2;
                x2(idx2_g) += 2;
                alpha = (+1.0i) * double(auxCur) * (config->thlV) / tmp[imaj];
                zgeru(&nDim, &nDim, &alpha, x1.data(), &inc, x2.data(), &inc, g.data(), &nDim);
                alpha = -alpha;
                zgeru(&nDim, &nDim, &alpha, x2.data(), &inc, x1.data(), &inc, g.data(), &nDim);
            }
        }
        // ... hsScheme == 1 case ...
    }
}
```

### Modified update() Loop

```cpp
DataType SpinlessVOperator::update(MatType& g) {
    double rand;
    bool flag;
    DataType signCur = 1.0;

    // Loop over all replicas
    for (int r = 0; r < nReplicas; r++) {
        if (singleMaj) {
            for (int i = 0; i < s[r]->size(); i++) {
                rand = rd->rdUniform01();
                singleFlipSingleMajorana(g, r, i, rand, flag, signCur);
            }
        } else {
            for (int i = 0; i < s[r]->size(); i++) {
                rand = rd->rdUniform01();
                singleFlip(g, r, i, rand, flag, signCur);
            }
        }
    }
    return signCur;
}
```

### Modified PfQMC Class

```cpp
// inc/pfqmc.h - Modified class

class PfQMC {
public:
    int stb;
    int nReplicas;        // NEW: number of replicas (default 1)
    int nDimSingle;       // NEW: single-replica dimension
    int nDim;             // CHANGED: = nReplicas * nDimSingle

    MatType g;            // nDim × nDim Green's function

    std::vector<Operator*> op_array;
    Operator* mixingOp;   // NEW: mixing operator at τ=0 (nullptr if nReplicas==1)

    int op_length;
    std::vector<bool> need_stabilization;
    int checkpoints;
    std::vector<UDT> udtL;
    std::vector<UDT> udtR;

    DataType sign;

    // CHANGED: constructor with optional nReplicas and mixing
    PfQMC(Spinless_tV* walker, int _stb = 10,
          int _nReplicas = 1, Operator* _mixingOp = nullptr);

    void rightInit();
    void leftInit();
    void rightSweep();
    void leftSweep();
    DataType getSignRaw();
};
```

### Modified rightSweep with Mixing

```cpp
void PfQMC::rightSweep() {
    MatType tmp = MatType::Identity(nDim, nDim);
    MatType Aseg = MatType::Identity(nDim, nDim);
    int curSeg = 0;
    DataType signCur;

    // NEW: Apply mixing operator at τ=0 (if present)
    if (mixingOp != nullptr) {
        mixingOp->left_propagate(g, tmp);
        mixingOp->left_multiply(Aseg, tmp);
        std::swap(Aseg, tmp);
    }

    // Rest is same as before, but operates on larger matrices
    for (int l = 0; l < op_length; l++) {
        signCur = op_array[l]->update(g);
        this->sign *= signCur;

        op_array[l]->left_multiply(Aseg, tmp);
        std::swap(Aseg, tmp);

        if (need_stabilization[(l + 1) % op_length]) {
            if (curSeg == 0) {
                udtR[curSeg] = UDT(Aseg);
            } else {
                udtR[curSeg] = Aseg * udtR[curSeg - 1];
            }
            Aseg = MatType::Identity(nDim, nDim);

            if (curSeg == (checkpoints - 1)) {
                udtR[curSeg].onePlusInv(g);
            } else {
                g = onePlusInv(udtL[curSeg + 1], udtR[curSeg]);
            }
            curSeg++;
        } else {
            op_array[l]->left_propagate(g, tmp);
        }
    }
}
```

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

---

## Testing Strategy for Direct Modification

### Phase 1: Backward Compatibility (Critical)

Before any replica features, verify `nReplicas=1` exactly matches original code:

```cpp
// Test 1: Constructor equivalence
SpinlessVOperator op_old(config, s_ptr, bondType, rd);           // original
SpinlessVOperator op_new(config, {s_ptr}, bondType, rd, 1);     // new with nReplicas=1

// Test 2: B matrix equivalence
assert((op_old.B - op_new.B).norm() < 1e-14);

// Test 3: Full simulation equivalence
// Run identical simulation with old and new code, compare:
// - Energy at each step
// - Sign at each step
// - Final Green's function
```

### Phase 2: Independent Replicas (No Mixing)

Test `nReplicas=4` with `mixingOp=nullptr` (identity mixing):

```cpp
// Should be equivalent to running 4 independent simulations
// Each replica evolves independently, results should be statistically equivalent

// Test: Run 4 separate single-replica simulations
// Compare average observables against 4-replica run
// They should match within statistical error
```

### Phase 3: With Mixing Operators

Test with actual M_12 and M_34 mixing:

```cpp
// Test 1: Verify off-diagonal blocks G_12, G_21 are non-zero after mixing
// Test 2: Verify zero blocks (between pair 12 and pair 34) remain zero
// Test 3: Compare against ED for small systems
```

### Quick Validation Checklist

| Test | Condition | Expected Result |
|------|-----------|-----------------|
| `nReplicas=1` energy | Compare to original | Exact match |
| `nReplicas=1` sign | Compare to original | Exact match |
| `nReplicas=4, no mixing` | Independent replicas | Statistical equivalence |
| `nReplicas=4, with mixing` | Off-diagonal blocks | Non-zero G_12, G_21 |
| `nReplicas=4, with mixing` | Zero blocks | G_13, G_14, etc. = 0 |

---

## Detailed Design: Ambiguous Parts

### 1. Acceptance Ratio for Coupled Replicas

When replicas within a pair are coupled by the mixing operator, the acceptance ratio for flipping an auxiliary field must account for the full pair-block structure.

**Key Insight**: The Green's function update formula `G' = G + (G - 2I) * Δ` where Δ is rank-2, still applies but now operates on the full `4*nDim × 4*nDim` matrix. The acceptance ratio `r` is computed from elements within a single replica's diagonal block.

**Why this works**: Even though G has off-diagonal blocks (G_12, G_21), flipping auxiliary field in replica 1 only directly modifies the diagonal block B_11. The Sherman-Morrison-Woodbury update propagates through the full matrix correctly.

```cpp
// Acceptance ratio computation (hsScheme == 0)
// For replica r with offset = r * nDimSingle

DataType computeAcceptanceRatio(const MatType& g, int replica,
                                 int idx1, int idx2, int idx3, int idx4,
                                 DataType thlV, int auxCur, double etaM) {
    int offset = replica * nDimSingle;
    int idx1_g = offset + idx1;
    int idx2_g = offset + idx2;
    int idx3_g = offset + idx3;
    int idx4_g = offset + idx4;

    // These G elements are from the diagonal block G_rr
    DataType tmp0 = 1.0 - 1.0i * thlV * double(auxCur) * g(idx1_g, idx2_g);
    DataType tmp1 = 1.0 - 1.0i * thlV * double(auxCur) * g(idx3_g, idx4_g);

    DataType r = tmp0 * tmp1;
    r += thlV * thlV * (g(idx1_g, idx3_g) * g(idx2_g, idx4_g)
                      - g(idx2_g, idx3_g) * g(idx1_g, idx4_g));
    r *= etaM;

    return r;
}
```

**Important**: The rank-2 update (`zgeru`) operates on the full matrix, which correctly updates both diagonal blocks (G_rr) and off-diagonal blocks (G_rs where s is the paired replica).

### 2. Sign Computation with Replicas

The sign computation in `getSignRaw()` must handle the enlarged matrix structure.

**Approach A: Product of Pair Signs**

Since the matrix is block-diagonal at the pair level, the total pfaffian factorizes:

```
Pf(G_full) = Pf(G_pair12) × Pf(G_pair34)
```

```cpp
DataType PfQMC::getSignRaw() {
    if (nReplicas == 1) {
        // Original single-replica code
        return getSignRawSingleReplica();
    }

    // For 4 replicas with pair structure
    DataType sign_pair12 = getSignRawForPair(0);  // replicas 0,1
    DataType sign_pair34 = getSignRawForPair(1);  // replicas 2,3

    return sign_pair12 * sign_pair34;
}

DataType PfQMC::getSignRawForPair(int pairIdx) {
    int pairOffset = pairIdx * 2 * nDimSingle;
    int pairDim = 2 * nDimSingle;

    // Extract pair block of the product of B matrices
    // Compute pfaffian for this pair
    // ... similar to original getSignRaw but on 2*nDimSingle dimension ...
}
```

**Approach B: Direct Computation**

Alternatively, compute on the full `4*nDim` matrix if the pfaffian routine handles block structure efficiently.

**Recommendation**: Use Approach A (factorized) for:
- Better numerical stability
- Easier debugging (can check each pair independently)
- Parallelization opportunity

### 3. Mixing Operator Specifics

The mixing operator M couples replicas at τ=0. Several physical choices exist:

#### Option A: SWAP Operator

Swaps Majorana fermions between replicas:

```cpp
// SWAP between replica 0 and 1 for site i:
// γ_i^{(0)} ↔ γ_i^{(1)}

MatType buildSwapMixing(int nDimSingle, const std::vector<int>& swapSites) {
    MatType M = MatType::Identity(2 * nDimSingle, 2 * nDimSingle);

    for (int site : swapSites) {
        // Swap rows/cols for this site between replica blocks
        // M[site, site] = 0, M[site, nDimSingle+site] = 1
        // M[nDimSingle+site, nDimSingle+site] = 0, M[nDimSingle+site, site] = 1
        M(site, site) = 0;
        M(site, nDimSingle + site) = 1;
        M(nDimSingle + site, nDimSingle + site) = 0;
        M(nDimSingle + site, site) = 1;
    }
    return M;
}
```

#### Option B: Rotation Operator

Continuous rotation between replicas:

```cpp
// Rotation by angle θ between replica 0 and 1 for site i:
// γ_i^{(0)} → cos(θ) γ_i^{(0)} + sin(θ) γ_i^{(1)}
// γ_i^{(1)} → -sin(θ) γ_i^{(0)} + cos(θ) γ_i^{(1)}

MatType buildRotationMixing(int nDimSingle, double theta,
                             const std::vector<int>& rotateSites) {
    MatType M = MatType::Identity(2 * nDimSingle, 2 * nDimSingle);
    double c = cos(theta);
    double s = sin(theta);

    for (int site : rotateSites) {
        M(site, site) = c;
        M(site, nDimSingle + site) = s;
        M(nDimSingle + site, site) = -s;
        M(nDimSingle + site, nDimSingle + site) = c;
    }
    return M;
}
```

#### Option C: Partial SWAP (for Entanglement Entropy)

For computing Rényi entropy, use partial transpose / partial SWAP:

```cpp
// SWAP only sites in region A between replicas
MatType buildPartialSwapMixing(int nDimSingle,
                                const std::vector<int>& regionA_sites) {
    return buildSwapMixing(nDimSingle, regionA_sites);
}
```

### 4. DenseOperator Replica Support

The kinetic operator (DenseOperator) must also support replicas:

```cpp
class DenseOperator : public Operator {
public:
    int nReplicas;
    int nDimSingle;
    int nDim;  // = nReplicas * nDimSingle

    MatType mat_single;      // single-replica matrix (nDimSingle × nDimSingle)
    MatType mat;             // full block-diagonal (nDim × nDim)
    MatType mat_inv;
    MatType g0, g0_inv;

    DenseOperator(const MatType& mat_single_, DataType _s, int _nReplicas = 1)
        : nReplicas(_nReplicas),
          nDimSingle(mat_single_.rows()),
          nDim(_nReplicas * mat_single_.rows())
    {
        mat_single = mat_single_;

        if (nReplicas == 1) {
            mat = mat_single;
        } else {
            // Build block-diagonal
            mat = MatType::Zero(nDim, nDim);
            for (int r = 0; r < nReplicas; r++) {
                mat.block(r*nDimSingle, r*nDimSingle, nDimSingle, nDimSingle) = mat_single;
            }
        }
        mat_inv = mat.inverse();

        // Build g0 block-diagonal
        MatType g0_single = (MatType::Identity(nDimSingle, nDimSingle) + mat_single).inverse()
                            * 2.0 - MatType::Identity(nDimSingle, nDimSingle);
        if (nReplicas == 1) {
            g0 = g0_single;
        } else {
            g0 = MatType::Zero(nDim, nDim);
            for (int r = 0; r < nReplicas; r++) {
                g0.block(r*nDimSingle, r*nDimSingle, nDimSingle, nDimSingle) = g0_single;
            }
        }
        g0_inv = g0.inverse();
        signOfWeight = _s;
    }
};
```

### 5. Walker Class (Spinless_tV) Replica Support

The walker class that builds the operator array needs replica support:

```cpp
class Spinless_tV {
public:
    std::vector<Operator*> op_array;
    int nDim;           // total dimension = nReplicas * nDimSingle
    int nDimSingle;     // single-replica dimension
    int nReplicas;

    // For 4 replicas, each operator in op_array handles all 4 replicas
    // Auxiliary fields are stored per-replica in each SpinlessVOperator
};
```

Example walker construction for 4 replicas:

```cpp
// In lattice-specific walker (e.g., Chain_tV)
Chain_tV(const SpinlessTvChainUtils* _config, rdGenerator* _rd, int _nReplicas = 1) {
    nReplicas = _nReplicas;
    nDimSingle = _config->nDim;
    nDim = nReplicas * nDimSingle;

    for (int l = 0; l < _config->l; l++) {
        // Create kinetic operator (same for all replicas)
        MatType B_K = /* kinetic matrix */;
        op_array.push_back(new DenseOperator(B_K, 1.0, nReplicas));

        // Create interaction operators with per-replica aux fields
        for (int bondType = 0; bondType < nBondTypes; bondType++) {
            std::vector<iVecType*> s_replicas(nReplicas);
            for (int r = 0; r < nReplicas; r++) {
                s_replicas[r] = new iVecType(nBonds);
                // Initialize aux fields randomly
                for (int b = 0; b < nBonds; b++) {
                    (*s_replicas[r])(b) = _rd->rdZ2();
                }
            }
            op_array.push_back(new SpinlessVOperator(_config, s_replicas,
                                                      bondType, _rd, nReplicas));
        }
    }
}
```

---

## Implementation Order (Revised)

1. **Phase 1: Backward Compatibility**
   - Modify `SpinlessVOperator` with `nReplicas` param (default=1)
   - Modify `DenseOperator` with `nReplicas` param (default=1)
   - Verify `nReplicas=1` matches original exactly

2. **Phase 2: Multi-Replica Without Mixing**
   - Test `nReplicas=4` with `mixingOp=nullptr`
   - Verify 4 independent replicas behave correctly
   - Check block-diagonal structure preserved

3. **Phase 3: Add Mixing Operator**
   - Implement `MixingOperator` class
   - Integrate into `PfQMC` sweeps
   - Test with simple SWAP mixing

4. **Phase 4: Sign Computation**
   - Implement pair-factorized `getSignRaw()`
   - Verify sign consistency

5. **Phase 5: Measurements**
   - Add per-replica measurement extraction
   - Add inter-replica correlation measurements

