# Development Guide: Mixing Operator Implementation

This guide documents the implementation of the Mixing Operator (Swap Operator) for the Replica method in PfQMC.

## 1. Completed Implementation Tasks

### 1.1 Mixing Operator Logic Correction
*   **Issue**: The original implementation of `right_multiply` and `inv_right_multiply` had swapped logic for $S$ and $S^{-1}$.
*   **Fix**: Corrected the matrix multiplication logic to ensure consistency with $B = A \times S$ and $B = A \times S^{-1}$.

### 1.2 Direct Swap Matrix Implementation
*   **Context**: The operator is updated to be a direct swap in imaginary time, represented by $B = i \sigma_m \tau_y$.
*   **Matrix Representation**:
    $$ B = \sigma_m \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} $$
    This operator swaps the replicas with a sign change.
*   **Skew-Symmetry**: Unlike the previous "Sqrt-Swap" ($\pi/4$ rotation), this operator is a **monomial matrix** (signed permutation). It preserves the skew-symmetry of the Green's function ($G^T = -G$) during propagation ($G' = B G B$).
*   **Affected Functions**: `left_multiply`, `right_multiply`. The superposition factors ($1/\sqrt{2}$) were removed.

### 1.3 Independent Auxiliary Fields
*   **Change**: Replaced the global `is_swapped` flag (which applied to an entire replica pair) with fine-grained auxiliary fields `sigma[p][m]`.
*   **Structure**: `sigma` is now a `vector<vector<int>>` where `sigma[p][m]` controls the swap state for mode `m` of replica pair `p`.
*   **Initialization**: Fields are initialized randomly to $\pm 1$.

### 1.4 Local Metropolis Updates
*   **Logic**: Refactored `update(g)` to perform a sweep over all auxiliary fields.
*   **Update Step**: Implemented `localUpdate` using the Sherman-Morrison-Woodbury formula (rank-2 update) for the specific `Delta` matrix associated with flipping $\sigma_m$.

### 1.5 Block Update Implementation
*   **Optimization**: Implemented `Blockupdate` structure to optimize the update process by exploiting the block-diagonal structure of the Green's function where applicable.
*   **Implementation**:
    *   Removed `assumeBlockDiagonal` flag; the block update is now the enforced default.
    *   Refactored `localUpdate` to use `zgeru` (BLAS rank-1 update) for efficient matrix updates on the specific replica pair blocks.
    *   This mimics the efficient update structure found in `inc/spinless_tV.h`.

## 2. Technical Details

### 2.1 Index Definitions for Local Updates
In the context of the `MixingOperator`, we define the indices as follows:
*   **$m$**: The spatial mode index (site index). This is a "local" index relative to a single replica.
*   **$j$**: The global matrix index corresponding to mode $m$ in **Replica 1** (or $\alpha$).
*   **$k$**: The global matrix index corresponding to mode $m$ in **Replica 2** (or $\beta$).

The operator $S_m^{12} = e^{\frac{\pi}{4} \sigma_m \gamma_m^1 \gamma_m^2}$ couples these two specific indices.

### 2.2 Ratio Calculation
The Metropolis acceptance ratio for flipping the auxiliary field $\sigma_m$ associated with this operator is given by:
$$ R = \frac{w'}{w} = 1 + \sigma_m G_{j,k} $$

Where:
*   $\sigma_m$ is the current value of the auxiliary field (before the flip).
*   $G_{j,k}$ is the Green's function element connecting the two replicas at site $m$.
    *   $G_{j,k} = \langle T \gamma_j \gamma_k \rangle$
    *   In the code: `G12 = g(idxA, idxB)` where `idxA` is the index for Replica 1, site $m$, and `idxB` is the index for Replica 2, site $m$.

This confirms that the update depends on the **inter-replica correlation** at the specific site $m$.

**Code Reference (`inc/mixing_operator.h`):**
```cpp
// getRatio function
int idxA = offsetA + m;  // Global index for Replica 1, site m (j)
int idxB = offsetB + m;  // Global index for Replica 2, site m (k)
DataType G12 = g(idxA, idxB);
return 1.0 + (double)s * G12;
```

## 3. Pending / Undone Tasks

### 3.1 Sign/Phase Calculation
*   **Current State**: The code calculates `signRatio *= (r / prob)`.
*   **Issue**: The ratio $r$ is calculated as `1 + s * G12`. This is a real number (assuming G is real). However, the general Pfaffian ratio can be complex or negative. The current logic might not fully capture the topological phase or sign changes associated with the swap operator, especially if $G_{12}$ becomes large or complex.
*   **Action**: Verify if `1 + s * G12` is sufficient for the phase, or if a more rigorous Pfaffian phase tracking is needed.

### 3.2 Stabilization Logic
*   **Current State**: `stabilizedLeftMultiply` falls back to full matrix multiplication because the block-diagonal assumption of `UDT` might be violated by the swap operator (which mixes rows between blocks).
*   **Action**: Optimize `stabilizedLeftMultiply` to handle the specific structure of the swap operator without full matrix multiplication if performance becomes a bottleneck.

### 3.3 TrivialSwapOperator
*   **Current State**: The `TrivialSwapOperator` class was updated to do nothing (Identity) to serve as a baseline.
*   **Action**: Confirm if this class is still needed or if it should be removed/updated to match the new `MixingOperator` interface.

### 3.4 Validation
*   **Action**: Run extensive tests to verify that the Sqrt-Swap implementation produces correct physical results (e.g., Second Renyi Entropy).

### 3.5 Mixing Operator Sweep in `src/pfqmc.cpp`
*   **Issue**: In `src/pfqmc.cpp`, the mixing operator is not fully integrated into the sweep logic.
    *   **Right Sweep**: There is a `!FIXME: update Mixing for the sign` comment. The mixing operator is propagated, but the Metropolis update (sweep) and the corresponding Green's function update are missing or incomplete.
    *   **Green's Function Update**: Because the mixing operator is not being "swept" (updated via Metropolis) in the right sweep, the Green's function `g` is not being updated to reflect changes in the mixing operator's auxiliary fields.
*   **Action**: Implement the Metropolis sweep for the mixing operator in `src/pfqmc.cpp` (both left and right sweeps) and ensure the Green's function is correctly updated.

## 4. Green's Function Definition and Update Formula

### 4.1 Shifted Green's Function
The `PfQMC` simulation uses a **shifted** definition for the Green's function matrix `g` stored in the walker:
161950 g = 2(I + B)^{-1} = I + G_{std} 161950
where {std} = 2(I+B)^{-1} - I$ is the standard physical Green's function (antisymmetric for Majoranas).

This shift affects how update vectors are constructed in the Sherman-Morrison formula.

### 4.2 Update Vector Construction
When performing a rank-2 update ' = G + U V^T$, the vectors $ and $ are typically derived from columns of $ or  \pm G$.

*   **Standard Formula**: If we need columns of {std}$, we must extract them from `g` as:
    161950 G_{std}.col(k) = g.col(k) - e_k 161950
    In code: `g.col(k)` followed by `g(k,k) -= 1.0`.

*   **Complement Formula**: If we need columns of  - G_{std}$ (common in Sherman-Morrison for +B$), we calculate:
    161950 I - G_{std} = I - (g - I) = 2I - g 161950
    In code: `2.0 * e_k - g.col(k)`.
    This explains the factor of **2.0** seen in `spinless_tV.h` (e.g., `x1(idx1) += 2.0`).

*   **Swap Operator Formula**: The derived kernel $ for the Swap operator uses vectors  = I - G_{std}$.
    Therefore, in `MixingOperator::localUpdate`, we construct:
    ```cpp
    cVecType c0 = -g.col(idxA);
    c0(idxA) += 2.0; // Represents (2I - g).col(idxA)
    ```
    This ensures consistency with the underlying physical Green's function.
