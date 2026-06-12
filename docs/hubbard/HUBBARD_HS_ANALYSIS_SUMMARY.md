# Hubbard Model HS Transformation Analysis Summary

## Key Findings

### 1. The aux2MajoranaIdx() Mapping is CORRECT

**Initial suspicion**: The mapping couples same-spin Majoranas incorrectly.

**Reality**: For on-site Hubbard interaction `U n_{i,↑} n_{i,↓}`, the HS transformation of the quartic term:
```
exp(dt U/4 · γ^1_{i,↑} γ^2_{i,↑} γ^1_{i,↓} γ^2_{i,↓})
```

decouples as:
```
(1/2) Σ_{s=±1} exp(λ s · i γ^1_{i,↑} γ^2_{i,↑}) · exp(λ s · i γ^1_{i,↓} γ^2_{i,↓})
```

This **requires** coupling `γ^1_{i,↑}` with `γ^2_{i,↑}` AND `γ^1_{i,↓}` with `γ^2_{i,↓}` using the **same** auxiliary field s_i.

The current implementation (`inc/hubbardChain.h:76-77`) correctly returns:
- For imaj=0: (4i+0, 4i+1) = (γ^1_{i,↑}, γ^2_{i,↑}) ✓
- For imaj=1: (4i+2, 4i+3) = (γ^1_{i,↓}, γ^2_{i,↓}) ✓

**Status**: ✅ CORRECT

### 2. Suspect Issue: HS Parameter λ

The critical question is the value of λ in:
```cpp
lambdaV = acosh(exp(factor * U * dt));
```

**Current implementation** (`inc/hubbardChain.h:40`):
```cpp
lambdaV = acosh(exp(0.5 * U * dt));  // factor = 1/2
```

**Theoretical derivation**:
The Hubbard interaction in Majorana basis is:
```
U n_{i,↑} n_{i,↓} = (U/4)[1 - i γ^1_{i,↑} γ^2_{i,↑} - i γ^1_{i,↓} γ^2_{i,↓} - γ^1_{i,↑} γ^2_{i,↑} γ^1_{i,↓} γ^2_{i,↓}]
```

The imaginary-time evolution is:
```
exp(-dt U n_{i,↑} n_{i,↓})
= exp(-dt U/4) · exp(i dt U/4 · γ^1_↑ γ^2_↑) · exp(i dt U/4 · γ^1_↓ γ^2_↓) · exp(dt U/4 · γ^1_↑ γ^2_↑ γ^1_↓ γ^2_↓)
```

Only the **quartic term** needs HS decoupling:
```
exp(dt U/4 · γ^1_↑ γ^2_↑ γ^1_↓ γ^2_↓) ≈ (1/2) Σ_{s=±1} exp(λ s · i γ^1_↑ γ^2_↑) · exp(λ s · i γ^1_↓ γ^2_↓)
```

For the discrete HS transformation:
```
exp(α AB) = (1/2) Σ_{s=±1} exp(λ s A) exp(λ s B)
```

where `cosh²(λ) = exp(α)`, giving:
```
λ = acosh(exp(α/2))
```

With `α = dt U/4`:
```
λ = acosh(exp(dt U/8))  // factor = 1/8
```

**Discrepancy**: Current uses factor=1/2, theory suggests factor=1/8.

### 3. Complications: Two-Point and Constant Terms

The full exponential includes:
- **Constant**: `exp(-dt U/4)`
- **Two-point up**: `exp(i dt U/4 · γ^1_↑ γ^2_↑)`
- **Two-point down**: `exp(i dt U/4 · γ^1_↓ γ^2_↓)`
- **Four-point**: `exp(dt U/4 · γ^1_↑ γ^2_↑ γ^1_↓ γ^2_↓)`

Different QMC implementations handle these terms differently:

**Option A**: Absorb two-point terms into kinetic part or handle separately, apply HS only to four-point.
- Would use `λ = acosh(exp(dt U/8))`

**Option B**: Package terms together in a specific way that changes the effective α.
- Could lead to a different factor

**Option C**: Use a modified HS scheme that handles multiple terms simultaneously.
- Might explain the factor=1/2

### 4. Comparison with Kitaev Chain

The Kitaev chain has NN interaction:
```
V n_i n_{i+1} = (V/4)[1 - i γ^1_i γ^2_i - i γ^1_{i+1} γ^2_{i+1} - γ^1_i γ^2_i γ^1_{i+1} γ^2_{i+1}]
```

**Exact same mathematical structure!**

Looking at `inc/spinless_tV.h:32`:
```cpp
lambdaV = acosh(exp(0.5 * V * dt));  // Also factor = 1/2
```

So both models use the same factor=1/2. This suggests:
- Either there's a consistent convention/transformation being used
- Or there's a consistent bug in both implementations

**Test**: If the Kitaev chain results are **verified correct** against ED or analytical solutions, then factor=1/2 is correct and my theoretical derivation missed something.

**Test**: If the Kitaev chain results are also **wrong**, then factor=1/2 is indeed the bug.

### 5. Energy Formula: Already Fixed

Previously, the energy formula was using incorrect signs. After fixing to match Kitaev pattern:
```cpp
energyFromGreensFunc(const MatType &g) {
    // E = E_kinetic + E_interaction
    // E_kinetic = (+,+,-) pattern for (hopping_↑, hopping_↓, chemical_potential)
    // E_interaction = -(U/4) <four-point interaction>
}
```

With current code (factor=1/2):
- Energy ≈ +1.395 (should be ≈ -1.95)

**Still wrong**, suggesting the λ parameter is indeed the issue.

## Recommended Action Plan

### Immediate Test: Try Different λ Factors

Modify `inc/hubbardChain.h:40` and test:

**Test 1**: Factor = 1/4
```cpp
lambdaV = acosh(exp(0.25 * U * dt));
```

**Test 2**: Factor = 1/8
```cpp
lambdaV = acosh(exp(0.125 * U * dt));
```

**Test 3**: Factor = 1/16
```cpp
lambdaV = acosh(exp(0.0625 * U * dt));
```

Run with same parameters (L=4, U=4.0) and compare energy against ED result (-1.947).

### Verify Kitaev Chain

Check if the Kitaev chain with V interaction has been validated:
- Look for test cases or papers comparing QMC vs exact/ED
- If Kitaev works with factor=1/2, understand why
- Apply same reasoning to Hubbard model

### Check for Additional Energy Contributions

Verify that the energy calculation includes ALL terms:
- Kinetic: hopping and chemical potential ✓
- Interaction constant: -U/4 per site
- Interaction two-point: ±(U/4) <i γ^1 γ^2> terms
- Interaction four-point: -(U/4) <γ^1_↑ γ^2_↑ γ^1_↓ γ^2_↓> ✓

The current energy formula may be **missing** the constant and two-point contributions!

### Check Sign Structure

Verify the sign during sampling remains close to 1:
```
sign ≈ signRaw ≈ 1.0
```

If sign becomes negative or complex, it indicates sign problem or fundamental issue with the HS transformation.

## Current Status

✅ **aux2MajoranaIdx mapping**: CORRECT
✅ **Energy formula pattern**: FIXED (matches Kitaev)
❌ **HS parameter λ**: SUSPECTED BUG (factor 1/2 vs 1/8 or other)
❓ **Energy constant/two-point terms**: NEED TO CHECK
❓ **Verification against working model**: PENDING

## Next Step

**Empirical test**: Run simulations with different λ factors to identify the correct value.
