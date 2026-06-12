# Replica Implementation Quick Reference

Branch: `replica-tV`

## What Changes

To implement 4 identical replicas of a t-V model (block-diagonal 4N×4N Green's function):

### Answer: **No separate class needed** - add `nReplica` parameter to existing classes

## Design Choice

**Direct parameter approach** (not wrapper/inheritance):
- Add `nReplica` parameter to Utils constructor
- Modify 6 functions to handle block-diagonal matrices
- All other code works automatically (dimension-agnostic)

## Functions to Modify

### 1. Utils Constructor (inc/kitaevChain.h:~15)
```cpp
// Before: nDim = 2 * Lx
// After:  nDim = 2 * Lx * nReplica
```

### 2. KineticGenerator() (inc/kitaevChain.h:47)
Build single-replica matrix, then replicate to block diagonal:
```cpp
MatType H_single = buildSingleReplicaKinetic();
for (int r = 0; r < nReplica; r++) {
    H.block(r*nDim_single, r*nDim_single, nDim_single, nDim_single) = H_single;
}
```

### 3. InteractionBGenerator() (inc/spinless_tV.h:91)
Same block-diagonal approach with **shared auxiliary fields**

### 4. Observable Functions (inc/kitaevChain.h:309+)
Average over replica blocks:
```cpp
for (int r = 0; r < nReplica; r++) {
    MatType g_replica = g.block(r*nDim_single, r*nDim_single, ...);
    result += computeObservable(g_replica);
}
return result / nReplica;
```

### 5. Main Function (main.cpp:~214)
Add `nReplica` parameter to CLI

### 6. aux2MajoranaIdx() (inc/kitaevChain.h:25)
**No changes needed** if using shared auxiliary fields

## Key Design Decisions

| Question | Decision | Rationale |
|----------|----------|-----------|
| Separate class? | **No** | Existing architecture is dimension-agnostic |
| Auxiliary fields? | **Shared across replicas** | Maintains replica symmetry, simpler implementation |
| Observable handling? | **Average over replicas** | Physical expectation value |
| Sign calculation? | **No changes initially** | Current implementation works, can optimize later |

## Matrix Structure

```
Green's function g (4N × 4N):

g = [ g_single    0         0         0      ]
    [ 0           g_single  0         0      ]
    [ 0           0         g_single  0      ]
    [ 0           0         0         g_single]
```

All operators (kinetic, interaction) have same block-diagonal structure.

## Backward Compatibility

Setting `nReplica = 1` (default) recovers original single-replica behavior.

## See Also

- Full design document: `REPLICA_IMPLEMENTATION_DESIGN.md`
- Implementation checklist in design doc
- Testing strategy in design doc
