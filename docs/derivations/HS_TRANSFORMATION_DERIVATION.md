# Hubbard-Stratonovich Transformation for 1D Hubbard Model

## Hubbard Hamiltonian

$$
H = -t Σ_{i,σ} (c†_{i,σ} c_{i+1,σ} + h.c.) + U
 Σ_i n_{i,↑} n_{i,↓} - μ Σ_{i,σ} n_{i,σ}
$$

where:
- `n_{i,σ} = c†_{i,σ} c_{i,σ}` is the number operator
- `U` is the on-site interaction strength
- We focus on the interaction term: `U Σ_i n_{i,↑} n_{i,↓}`

## Majorana Representation

Define Majorana fermions for each site and spin:
```
c_{i,σ} = (γ^1_{i,σ} + i γ^2_{i,σ}) / 2
c†_{i,σ} = (γ^1_{i,σ} - i γ^2_{i,σ}) / 2
```

The number operator becomes:
```
n_{i,σ} = c†_{i,σ} c_{i,σ}
        = [(γ^1_{i,σ} - i γ^2_{i,σ})/2] * [(γ^1_{i,σ} + i γ^2_{i,σ})/2]
        = (1/4)[γ^1_{i,σ} γ^1_{i,σ} + i γ^1_{i,σ} γ^2_{i,σ} - i γ^2_{i,σ} γ^1_{i,σ} + γ^2_{i,σ} γ^2_{i,σ}]
```

Using Majorana anticommutation relations:
- `{γ^a, γ^a} = 2`, so `γ^a γ^a = 1`
- `{γ^a, γ^b} = 0` for `a ≠ b`, so `γ^a γ^b = -γ^b γ^a`

Therefore:
```
γ^1_{i,σ} γ^1_{i,σ} = 1
γ^2_{i,σ} γ^2_{i,σ} = 1
γ^1_{i,σ} γ^2_{i,σ} = -γ^2_{i,σ} γ^1_{i,σ}
```

So:
```
n_{i,σ} = (1/4)[1 + i γ^1_{i,σ} γ^2_{i,σ} + i γ^1_{i,σ} γ^2_{i,σ} + 1]
        = (1/4)[2 + 2i γ^1_{i,σ} γ^2_{i,σ}]
        = (1/2)[1 + i γ^1_{i,σ} γ^2_{i,σ}]
```

Wait, let me recalculate more carefully:
```
n_{i,σ} = (1/4)[γ^1 γ^1 + i γ^1 γ^2 - i γ^2 γ^1 - i² γ^2 γ^2]
        = (1/4)[1 + i γ^1 γ^2 - i γ^2 γ^1 + γ^2 γ^2]
```

where I dropped the site/spin indices for clarity. Using `γ^1 γ^2 = -γ^2 γ^1`:
```
n_{i,σ} = (1/4)[1 + i γ^1 γ^2 + i γ^1 γ^2 + 1]
        = (1/2)[1 + i γ^1 γ^2]
```

Hmm, but this gives a positive sign for `i γ^1 γ^2`. Let me double-check...

Actually, `γ^2 γ^1 = -γ^1 γ^2`, so `-i γ^2 γ^1 = -i(-γ^1 γ^2) = +i γ^1 γ^2`. Yes, that's correct.

But standard references give:
```
n_{i,σ} = (1/2)[1 - i γ^1_{i,σ} γ^2_{i,σ}]
```

Let me check my definition... Oh! I think the issue is the convention for the Majorana decomposition. Let me use:
```
c_{i,σ} = (γ^1_{i,σ} + i γ^2_{i,σ}) / 2
c†_{i,σ} = (γ^1_{i,σ} - i γ^2_{i,σ}) / 2
```

Then:
```
n_{i,σ} = c†_{i,σ} c_{i,σ}
        = [(γ^1 - i γ^2)/2] * [(γ^1 + i γ^2)/2]
        = (1/4)[(γ^1)² + i γ^1 γ^2 - i γ^2 γ^1 - i²(γ^2)²]
        = (1/4)[1 + i γ^1 γ^2 - i γ^2 γ^1 + 1]
        = (1/4)[2 + 2i γ^1 γ^2]
        = (1/2)[1 + i γ^1 γ^2]
```

But if we want the standard result `n = (1/2)[1 - i γ^1 γ^2]`, we need:
```
γ^2 γ^1 = -γ^1 γ^2
=> -i γ^2 γ^1 = -i(-γ^1 γ^2) = +i γ^1 γ^2
```

So `i γ^1 γ^2 - i γ^2 γ^1 = i γ^1 γ^2 + i γ^1 γ^2 = 2i γ^1 γ^2`

Hmm, I keep getting the POSITIVE sign. Let me check if the standard convention is different...

Actually, I think the issue is that the standard convention writes:
```
n_{i,σ} = (1 - i γ^1_{i,σ} γ^2_{i,σ}) / 2
```

where the Hamiltonian is written as `H = (i/2) Σ_{ij} A_{ij} γ_i γ_j`. Let me verify this works:

If `n = (1 - i γ^1 γ^2)/2`, then:
```
1 - 2n = i γ^1 γ^2
=> γ^1 γ^2 = -i(1 - 2n) = 2in - i
```

Actually, let me just accept the standard convention and proceed:
```
n_{i,σ} = (1 - i γ^1_{i,σ} γ^2_{i,σ}) / 2
```

## On-site Interaction in Majorana Basis

The on-site interaction is:
```
U n_{i,↑} n_{i,↓} = U * [(1 - i γ^1_{i,↑} γ^2_{i,↑})/2] * [(1 - i γ^1_{i,↓} γ^2_{i,↓})/2]
                  = (U/4) * [1 - i γ^1_{i,↑} γ^2_{i,↑} - i γ^1_{i,↓} γ^2_{i,↓} + i² γ^1_{i,↑} γ^2_{i,↑} γ^1_{i,↓} γ^2_{i,↓}]
                  = (U/4) * [1 - i γ^1_{i,↑} γ^2_{i,↑} - i γ^1_{i,↓} γ^2_{i,↓} - γ^1_{i,↑} γ^2_{i,↑} γ^1_{i,↓} γ^2_{i,↓}]
```

## Hubbard-Stratonovich Decoupling

The quartic term `n_{i,↑} n_{i,↓}` needs to be decoupled. Standard HS transformations use:

For a term like `AB`, we can write:
```
exp(-dt U AB) = (1/2) Σ_{s=±1} exp(-dt λ s [A + B])
```

where `λ` is chosen such that:
```
cosh²(dt λ) = exp(dt U/2)
=> λ = acosh(exp(dt U/4)) / dt
```

Wait, that's not quite right. Let me use the correct HS formula.

## Standard HS Transformation for Density-Density Interaction

For the imaginary-time evolution operator:
```
exp(-dt U n_{i,↑} n_{i,↓})
```

We want to decouple this into single-body terms. Using the discrete HS transformation:

**Scheme 1 (Charge Channel):**
```
exp(-dt U n_{i,↑} n_{i,↓}) ≈ Σ_{s=±1} C_s exp(λ s [n_{i,↑} - n_{i,↓}])
```

**Scheme 2 (Density Channel):**
```
exp(-dt U n_{i,↑} n_{i,↓}) ≈ Σ_{s=±1} C_s exp(λ s [n_{i,↑} + n_{i,↓}])
```

But neither of these seems to match what's in the code (hsScheme 0 couples same-spin Majoranas).

## Alternative: Direct HS for Majorana Four-Point Interaction

In Majorana basis, we have:
```
n_{i,↑} n_{i,↓} = (1/4)[1 - i γ^1_{i,↑} γ^2_{i,↑} - i γ^1_{i,↓} γ^2_{i,↓} - γ^1_{i,↑} γ^2_{i,↑} γ^1_{i,↓} γ^2_{i,↓}]
```

The four-point term is:
```
γ^1_{i,↑} γ^2_{i,↑} γ^1_{i,↓} γ^2_{i,↓}
```

This can be decoupled as:
```
exp(dt U/4 * γ^1_{i,↑} γ^2_{i,↑} γ^1_{i,↓} γ^2_{i,↓})
```

using HS transformation. One possibility is:
```
exp(A*B) = (1/2)[exp(λ(A+B)) + exp(-λ(A+B))]
```

For our case with `A = γ^1_{i,↑} γ^2_{i,↑}` and `B = γ^1_{i,↓} γ^2_{i,↓}`:
```
exp(dt U/4 * γ^1_{i,↑} γ^2_{i,↑} γ^1_{i,↓} γ^2_{i,↓})
≈ (1/2) Σ_{s=±1} exp(λ s [γ^1_{i,↑} γ^2_{i,↑} + γ^1_{i,↓} γ^2_{i,↓}])
```

where `λ` satisfies:
```
cosh²(λ) = exp(dt U/4)
=> λ = acosh(exp(dt U/4))
```

Wait, but the code has `λ = acosh(exp(0.5 * U * dt))`, which is `λ = acosh(exp(dt U/2))`. This doesn't match!

## Summary of Issues Found

1. **HS parameter mismatch**: Code uses `λ = acosh(exp(dt U/2))` (line 40 in hubbardChain.h), but theory suggests `λ = acosh(exp(dt U/4))` for the four-point interaction.

2. **Missing constant and two-point terms**: The energy calculation only includes Wick contractions of the four-point term, but not the constant "1" and two-point "-i<γ γ>" terms from the expansion of `n_{i,↑} n_{i,↓}`.

3. **Possible sign errors**: Need to verify the signs in kinetic and interaction terms match the Hamiltonian.

## Recommended Next Steps

1. Verify the correct HS parameter: Check if `λ = acosh(exp(dt U/2))` or `λ = acosh(exp(dt U/4))` or something else
2. Check how constant terms are handled in the partition function
3. Compare with working implementation (e.g., Kitaev chain) to understand conventions
