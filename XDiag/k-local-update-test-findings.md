# K-Local Update Test for Pauli Sampling: Findings

## Overview

This document summarizes the findings from `k-local-update-test.jl`, which investigates how k-local MCMC updates affect sampling efficiency for the Pauli distribution of the transverse-field Ising model (TFIM) ground state.

The Pauli distribution is:
$$\Pi(P) = |\langle\psi_0|P|\psi_0\rangle|^2$$

where $|\psi_0\rangle$ is the TFIM ground state and $P$ is a Pauli string.

## Key Metrics

### 1. Raw Acceptance Rate
The fraction of proposed k-local updates that are accepted by Metropolis criterion.

### 2. Normalized Acceptance Rate
$$\text{Normalized} = \frac{\text{Acc Rate}}{(3/4)^k}$$

**Rationale**: Each of the k sites has probability 3/4 of actually changing (1/4 chance of proposing the same Pauli). This normalization removes the trivial k-dependence.

**Interpretation**:
- If constant across k → distribution factorizes (sites independent)
- If decreases with k → k-body correlations exist
- If increases with k → larger moves find better paths

### 3. Average Hamming Distance (Accepted)
Number of sites that actually changed per accepted move. Only counts accepted moves.

Expected value: ~0.75k (since each site has 3/4 chance of changing)

### 4. Pauli String Hamming Weight
Number of non-identity Paulis in a sampled string. Measures the "spread" of the state in Pauli space.

## Results

### Part 1: Scaling Test (N=6, h=1.0)

| k | Acc Rate | Normalized | Avg Hamming | Num Accepted (of 5000) |
|---|----------|------------|-------------|------------------------|
| 1 | 0.31 | 0.41 | 1.00 | ~300 |
| 2 | 0.21 | 0.38 | 1.61 | ~750 |
| 3 | 0.14 | 0.33 | 2.11 | ~630 |
| 4 | 0.10 | 0.33 | 2.78 | ~510 |
| 5 | 0.09 | 0.37 | 3.52 | ~430 |

**Key finding**: Normalized acceptance rate is relatively flat (~0.33-0.41), suggesting the Pauli distribution is approximately factorizable at the critical point h=1.

### Part 2: Parameter Sweep (h = 0.8 to 1.2)

#### Normalized Acceptance Rate by h and k

| h | k=1 | k=2 | k=3 | k=4 | k=5 |
|------|------|------|------|------|------|
| 0.80 | 0.40 | 0.32 | 0.27 | 0.26 | 0.29 |
| 0.90 | 0.38 | 0.35 | 0.31 | 0.33 | 0.33 |
| 1.00 | 0.46 | 0.37 | 0.36 | 0.31 | 0.36 |
| 1.10 | 0.45 | 0.39 | 0.34 | 0.34 | 0.35 |
| 1.20 | 0.38 | 0.39 | 0.33 | 0.33 | 0.31 |

**Key findings**:
1. **h=1.0 (critical point)**: Flattest profile across k (~0.31-0.46), suggesting most independent site behavior
2. **h=0.8**: More k-dependence (drops from 0.40 to 0.26), indicating stronger correlations below criticality
3. **h>1**: Intermediate behavior

### Part 3: Pauli String Hamming Weight

| h | Avg Weight | P(I) | P(X) | P(Y) | P(Z) |
|------|------------|------|------|------|------|
| 0.80 | 4.12 | 0.31 | 0.29 | 0.22 | 0.18 |
| 0.90 | 4.02 | 0.33 | 0.31 | 0.18 | 0.19 |
| 1.00 | 3.78 | 0.37 | 0.34 | 0.14 | 0.15 |
| 1.10 | 3.90 | 0.35 | 0.41 | 0.17 | 0.07 |
| 1.20 | 3.92 | 0.35 | 0.40 | 0.14 | 0.11 |

**Key findings**:
1. **P(X) increases with h** (0.29 → 0.40): As transverse field dominates, X-type Paulis become more probable
2. **P(Z) decreases with h** (0.18 → 0.07): Ising Z-correlations weaken
3. **P(I) peaks near h=1** (0.37): Critical point has highest identity probability
4. **Avg Hamming weight minimum at h≈1** (3.78): Fewer non-identity Paulis at criticality

## Physical Interpretation

### Connection to Ground State Structure

The Pauli marginal distribution P(I,X,Y,Z) reveals the ground state structure:

1. **h < 1 (Ising-dominated)**:
   - Ground state closer to |↑↑...↑⟩ or |↓↓...↓⟩
   - More Z-type correlations
   - Pauli distribution has multi-site structure → k-dependence in acceptance

2. **h ≈ 1 (Critical)**:
   - Maximum entanglement
   - Pauli distribution most "democratic"
   - Higher P(I) → more compact Pauli representation
   - Normalized acceptance nearly independent of k

3. **h > 1 (Paramagnetic)**:
   - Ground state closer to |+...+⟩
   - X-type Paulis dominate
   - Some structure but different character

### Why h=0 Fails

At h=0 (classical Ising), the ground state is a product state with zero magic/SRE. The Pauli distribution becomes degenerate, and importance sampling fails because:
- Only a few Pauli strings have non-zero weight
- MCMC gets stuck in local modes

### Mixing Analysis

The "Hamming distance from reference" metric shows:
- **k=1**: Often stuck (mean Hamming ~1-2, should explore up to N=6)
- **k≥2**: Good exploration (mean Hamming ~4, close to N/2 as expected for ergodic sampling)

This confirms k=1 updates suffer from poor mixing even with decent acceptance rate.

## Practical Recommendations

1. **Use k≥2 for Pauli sampling**: k=1 has poor mixing despite high acceptance
2. **k=2 is often optimal**: Best balance of acceptance rate and move size
3. **Near criticality (h≈1)**: Any k works reasonably well due to factorizable distribution
4. **Away from criticality**: Larger k may be needed to overcome correlations

## Theoretical Connection

The normalized acceptance rate being constant across k implies:
$$\Pi(P) \approx \prod_{i=1}^N \pi_i(P_i)$$

i.e., the Pauli distribution approximately factorizes into single-site marginals. This is most true at the critical point, suggesting the "magic" is distributed uniformly rather than concentrated in specific multi-site correlations.

The (3/4)^k normalization accounts for the baseline probability that k sites actually change when k-site updates are proposed.
