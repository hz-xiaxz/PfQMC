# Correct Hubbard-Stratonovich Transformation for 1D Hubbard Model in Majorana Basis

## Problem Statement

The on-site Hubbard interaction is:
```
H_int = U Σ_i n_{i,↑} n_{i,↓}
```

where `n_{i,σ} = c†_{i,σ} c_{i,σ}` is the number operator for spin σ at site i.

## Step 1: Majorana Representation

Define Majorana fermions:
```
c_{i,σ} = (γ^1_{i,σ} + i γ^2_{i,σ}) / 2
c†_{i,σ} = (γ^1_{i,σ} - i γ^2_{i,σ}) / 2
```

where γ are Hermitian operators satisfying:
- `{γ^a_i, γ^b_j} = 2 δ_ij δ_ab`
- `(γ^a_i)† = γ^a_i`

## Step 2: Number Operator in Majorana Basis

```
n_{i,σ} = c†_{i,σ} c_{i,σ}
        = [(γ^1_{i,σ} - i γ^2_{i,σ})/2] · [(γ^1_{i,σ} + i γ^2_{i,σ})/2]
        = (1/4)[(γ^1_{i,σ})² + i γ^1_{i,σ} γ^2_{i,σ} - i γ^2_{i,σ} γ^1_{i,σ} - i²(γ^2_{i,σ})²]
```

Using `(γ^a)² = 1` and `γ^2 γ^1 = -γ^1 γ^2`:
```
n_{i,σ} = (1/4)[1 + i γ^1_{i,σ} γ^2_{i,σ} - i(-γ^1_{i,σ} γ^2_{i,σ}) + 1]
        = (1/4)[2 + 2i γ^1_{i,σ} γ^2_{i,σ}]
        = (1/2)[1 + i γ^1_{i,σ} γ^2_{i,σ}]
```

**But wait!** The standard convention in the literature uses:
```
n_{i,σ} = (1/2)[1 - i γ^1_{i,σ} γ^2_{i,σ}]
```

This discrepancy arises from different Majorana conventions. Let me verify which convention the codebase uses...

Looking at `inc/hubbardChain.h:21-22`, the code uses:
```cpp
// n_{i,σ} = (1 - i γ^1 γ^2) / 2
```

So the codebase convention is:
```
n_{i,σ} = (1 - i γ^1_{i,σ} γ^2_{i,σ}) / 2
```

This means the Majorana representation must be:
```
c_{i,σ} = (γ^2_{i,σ} + i γ^1_{i,σ}) / 2  [Note: swapped from above!]
c†_{i,σ} = (γ^2_{i,σ} - i γ^1_{i,σ}) / 2
```

Let me verify:
```
n = c† c = [(γ^2 - i γ^1)/2] · [(γ^2 + i γ^1)/2]
         = (1/4)[(γ^2)² + i γ^2 γ^1 - i γ^1 γ^2 + (γ^1)²]
         = (1/4)[1 + i γ^2 γ^1 + i γ^2 γ^1 + 1]  [using γ^1 γ^2 = -γ^2 γ^1]
         = (1/4)[2 + 2i γ^2 γ^1]
         = (1/2)[1 + i γ^2 γ^1]
         = (1/2)[1 - i γ^1 γ^2]  ✓
```

Good! The convention is settled.

## Step 3: Interaction Term Expansion

```
U n_{i,↑} n_{i,↓} = U · [(1 - i γ^1_{i,↑} γ^2_{i,↑})/2] · [(1 - i γ^1_{i,↓} γ^2_{i,↓})/2]
                  = (U/4) · [1 - i γ^1_{i,↑} γ^2_{i,↑}][1 - i γ^1_{i,↓} γ^2_{i,↓}]
```

Expanding the product:
```
= (U/4) · [1
           - i γ^1_{i,↑} γ^2_{i,↑}
           - i γ^1_{i,↓} γ^2_{i,↓}
           + i² γ^1_{i,↑} γ^2_{i,↑} γ^1_{i,↓} γ^2_{i,↓}]

= (U/4) · [1
           - i γ^1_{i,↑} γ^2_{i,↑}
           - i γ^1_{i,↓} γ^2_{i,↓}
           - γ^1_{i,↑} γ^2_{i,↑} γ^1_{i,↓} γ^2_{i,↓}]
```

## Step 4: Hubbard-Stratonovich Transformation

The imaginary-time evolution operator is:
```
exp(-dt · U n_{i,↑} n_{i,↓}) = exp(-dt U/4 · [1 - i γ^1_{i,↑} γ^2_{i,↑} - i γ^1_{i,↓} γ^2_{i,↓} - γ^1_{i,↑} γ^2_{i,↑} γ^1_{i,↓} γ^2_{i,↓}])
```

This can be decomposed:
```
= exp(-dt U/4) · exp(i dt U/4 · γ^1_{i,↑} γ^2_{i,↑}) · exp(i dt U/4 · γ^1_{i,↓} γ^2_{i,↓}) · exp(dt U/4 · γ^1_{i,↑} γ^2_{i,↑} γ^1_{i,↓} γ^2_{i,↓})
```

The first three factors are **constant** or **quadratic** (already in single-body form). The problematic term is the **quartic** four-Majorana interaction:
```
exp(dt U/4 · γ^1_{i,↑} γ^2_{i,↑} γ^1_{i,↓} γ^2_{i,↓})
```

### HS Decoupling of Quartic Term

For a product of two operators A and B, the discrete HS transformation is:
```
exp(α A B) = (1/2) Σ_{s=±1} exp(λ s A) exp(λ s B)
```

where λ satisfies:
```
cosh²(λ) = exp(α)
=> λ = acosh(√[exp(α)])
=> λ = acosh(exp(α/2))
```

In our case:
- `α = dt U/4`
- `A = γ^1_{i,↑} γ^2_{i,↑}`
- `B = γ^1_{i,↓} γ^2_{i,↓}`

Therefore:
```
λ = acosh(exp(dt U/8))
```

**Wait!** Let me double-check this formula. The standard discrete HS transformation is:
```
exp(α A B) = (1/2) Σ_{s=±1} exp(λ s [A + B])
```

where:
```
cosh²(λ) = exp(α)
```

So:
```
λ = acosh(exp(α/2))
```

With `α = dt U/4`:
```
λ = acosh(exp(dt U/8))
```

### Alternative HS Scheme

There's also a "density-density" style HS transformation that treats the product differently:
```
exp(-dt U n_↑ n_↓)
= exp(-dt U/4) · exp(dt U/4 (n_↑ + n_↓ - 2n_↑ n_↓))
```

Using the identity:
```
exp(α(A + B - 2AB)) = (1/2) Σ_{s=±1} exp(λ s [A - B])
```

This would give a **different** coupling pattern.

## Step 5: Correct Majorana Coupling

### Scheme 0 (Standard Decoupling):

Decouple the four-point term as:
```
exp(dt U/4 · γ^1_{i,↑} γ^2_{i,↑} γ^1_{i,↓} γ^2_{i,↓})
= (1/2) Σ_{s=±1} exp(λ s [i γ^1_{i,↑} γ^2_{i,↑} + i γ^1_{i,↓} γ^2_{i,↓}])
```

where:
```
λ = acosh(exp(dt U/8))
```

**Key observation:** This couples `γ^1_{i,↑}` with `γ^2_{i,↑}` AND `γ^1_{i,↓}` with `γ^2_{i,↓}` with the **same auxiliary field** s.

In the codebase indexing (site i, spin σ, Majorana index α):
- Majorana index for (i, ↑, 1): `4*i + 0`
- Majorana index for (i, ↑, 2): `4*i + 1`
- Majorana index for (i, ↓, 1): `4*i + 2`
- Majorana index for (i, ↓, 2): `4*i + 3`

So `aux2MajoranaIdx(i, imaj, ...)` should return:
- **First pair**: `(4*i + 0, 4*i + 1)` - couples γ^1_{i,↑} with γ^2_{i,↑}
- **Second pair**: `(4*i + 2, 4*i + 3)` - couples γ^1_{i,↓} with γ^2_{i,↓}

But these are coupled by the **same** auxiliary field s_i at site i!

### Current Code Issues

Looking at `inc/hubbardChain.h:76-77`:
```cpp
idx1 = i * 4 + imaj * 2 + 0;
idx2 = i * 4 + imaj * 2 + 1;
```

This returns:
- For `imaj=0`: `(4*i+0, 4*i+1)` - γ^1_{i,↑}, γ^2_{i,↑} ✓
- For `imaj=1`: `(4*i+2, 4*i+3)` - γ^1_{i,↓}, γ^2_{i,↓} ✓

This is **correct** for Scheme 0!

So why is the energy wrong?

## Step 6: Check the HS Parameter

Looking at `inc/hubbardChain.h:40`:
```cpp
lambdaV = acosh(exp(0.5 * V * dt));
```

Here `V` is the Hubbard `U`. So:
```
λ = acosh(exp(dt U/2))
```

But from our derivation, we should have:
```
λ = acosh(exp(dt U/8))
```

**This is the bug!** The factor is off by 4×.

### Why the Factor of 4?

The interaction term is:
```
U n_↑ n_↓ = (U/4)[1 - i γ^1_↑ γ^2_↑ - i γ^1_↓ γ^2_↓ - γ^1_↑ γ^2_↑ γ^1_↓ γ^2_↓]
```

The quartic term contributes:
```
-(U/4) γ^1_↑ γ^2_↑ γ^1_↓ γ^2_↓
```

In imaginary time:
```
exp(dt U/4 · γ^1_↑ γ^2_↑ γ^1_↓ γ^2_↓)
```

So `α = dt U/4`, giving:
```
λ = acosh(exp(dt U/8))
```

### Comparison with Kitaev Chain

For the Kitaev chain with NN interaction `V n_i n_{i+1}`:
```
V n_i n_{i+1} = V · [(1 - i γ^1_i γ^2_i)/2] · [(1 - i γ^1_{i+1} γ^2_{i+1})/2]
              = (V/4)[1 - i γ^1_i γ^2_i - i γ^1_{i+1} γ^2_{i+1} - γ^1_i γ^2_i γ^1_{i+1} γ^2_{i+1}]
```

Same structure! So the factor should also be `λ = acosh(exp(dt V/8))` for Kitaev chain.

But looking at `inc/spinless_tV.h`, let me check what the actual parameter is...

## Step 7: Full Exponential Including All Terms

Actually, we need to be more careful. The full imaginary-time operator is:
```
exp(-dt U n_↑ n_↓)
= exp(-dt U/4 · [1 - i γ^1_↑ γ^2_↑ - i γ^1_↓ γ^2_↓ - γ^1_↑ γ^2_↑ γ^1_↓ γ^2_↓])
```

Using Trotter decomposition (to first order in dt):
```
≈ exp(-dt U/4) · exp(i dt U/4 · γ^1_↑ γ^2_↑) · exp(i dt U/4 · γ^1_↓ γ^2_↓) · exp(dt U/4 · γ^1_↑ γ^2_↑ γ^1_↓ γ^2_↓)
```

The constant factor `exp(-dt U/4)` cancels in ratios.

The two-point terms `exp(i dt U/4 · γ^1_σ γ^2_σ)` are absorbed into the kinetic operator or can be handled separately.

**In practice**, many QMC codes absorb these terms differently. Let me check how the code constructs operators...

Looking at `inc/hubbardChain.h:88-106` (the walker constructor), I see:
- Kinetic operators: `exp(-dt/2 · H_K)` at boundaries, `exp(-dt · H_K)` in bulk
- Interaction operators: One per site, applied in sequence

The interaction operator `SpinlessVOperator` generates a matrix from the auxiliary fields. Let me check the matrix generator in `inc/spinless_tV.h`...

## Conclusion and Recommendation

**However**, there's complexity in how the constant and two-point terms are handled. The correct approach depends on:

1. Whether constant factors are tracked separately
2. Whether two-point terms are absorbed into kinetic or interaction operators
3. The exact Trotter decomposition scheme used

**Next Steps:**

1. Check `inc/spinless_tV.h` to see how `InteractionBGenerator()` constructs the interaction matrix
2. Compare with the Kitaev chain parameter to see if there's a consistent convention
3. Test different values of λ empirically to find the correct factor
4. Consider that the "correct" factor might involve additional considerations from the Trotter splitting

Let me investigate the interaction matrix generator next...
