# Single Flip Update in PFQMC

This document explains the `singleFlip` function in `SpinlessVOperator`, which performs the local update of auxiliary fields.

## Mathematical Context

The interaction term is decomposed using a discrete Hubbard-Stratonovich (HS) transformation. For a single bond, the exponential of the interaction is written as:

$$ e^{-\Delta\tau V n_i n_j} \propto \sum_{s=\pm 1} e^{\gamma s (n_i - n_j)^2} \dots $$

In the Majorana representation, this becomes an operator $B(s)$ that depends on the auxiliary field $s$.

When we propose a flip $s \to -s$, the operator changes $B(s) \to B(-s)$. The acceptance ratio $R$ for this move is given by the ratio of the Pfaffians (or determinants) of the new and old Green's functions.

## The `singleFlip` Algorithm

The function `singleFlip` calculates this ratio and updates the Green's function if accepted.

### 1. Calculate Acceptance Ratio $R$

The ratio of weights for flipping field $s_k$ is:

$$ R = \frac{\text{Pf}(1 + B(-s_k) \dots)}{\text{Pf}(1 + B(s_k) \dots)} = \det(1 + (B(-s_k) B(s_k)^{-1} - 1) G)_{kk} $$

In the code, this is computed efficiently using the Sherman-Morrison formula logic. The change $B(-s)B(s)^{-1}$ is a low-rank update (rank-2 for standard scheme).

The code computes:
```cpp
tmp[0] = 1.0 - i * tanh(lambda/2) * s * G(idx1, idx2)
tmp[1] = 1.0 - i * tanh(lambda/2) * s * G(idx3, idx4)
r = tmp[0] * tmp[1] + ... (cross terms)
```

### 2. Metropolis Step

We generate a random number `rand` $\in [0, 1)$.
- **Accept** if `rand < |R|`.
- **Reject** otherwise.

### 3. Update Green's Function (Sherman-Morrison)

If accepted, we must update the Green's function $G$ to $G_{new}$. Since the change in the matrix is low-rank (rank-2), we can update $G$ in $O(N^2)$ time instead of recomputing it ($O(N^3)$).

The update formula is roughly:
$$ G_{new} = G - \frac{G u v^T G}{1 + v^T G u} $$

In the code, this is implemented via `zgeru` (BLAS rank-1 update) calls:
```cpp
// Update columns x1, x2
x1 = -G.col(idx1) + ...
x2 = -G.col(idx2) + ...

// Rank-1 updates
G += alpha * x1 * x2^T
G += -alpha * x2 * x1^T
```

### 4. Coupled Replica Updates

In the 4-replica scheme, replicas are coupled in pairs (0-1 and 2-3). The update is proposed for a pair of replicas simultaneously.

1.  **Generate Random Number**: A single random number `rand` is generated for the pair.
2.  **Calculate Ratios**: Calculate $R_A$ for replica A and $R_B$ for replica B using `getRatio`.
3.  **Combined Ratio**: $R_{total} = R_A \times R_B$.
4.  **Acceptance**: Accept if `rand < |R_{total}|`.
5.  **Update**: If accepted, update both replicas.

### 5. Block Diagonal Optimization

When `assumeBlockDiagonal` is true, we assume the Green's function has a block-diagonal structure (2x2 blocks in replica space).

-   **`updateBlockPair`**: This helper function performs the update for a pair of replicas efficiently by only updating the relevant 2N x 2N block of the Green's function.
    > [!IMPORTANT]
    > The 2N * 2N block update strategy needs to be reconsidered. The current assumption might not be fully optimal or correct for all cases.

-   **`updateReplica`**: This helper function performs the standard dense update for a single replica.

### 6. Update State

Finally, we flip the auxiliary field and the corresponding operator matrix elements for both replicas in the pair if accepted.

