# Development Memory: Mixing Operator Update

## 1. What's Done
* Fixed skew-symmetry violation in `MixingOperator::localUpdate`.
* Correctly applied the Sherman-Morrison-Woodbury (SMW) formula for the rank-2 update $\Delta B$.
* Verified the update locally using `reproduce_mixing_update.cpp`. Tests pass for both ground truth match ($\sim 10^{-14}$ error) and round-trip identity.

## 2. What Fix Applied
* **Reverted `$Q^T O Q$` form:** The symmetric form failed because the shifted Green's function $g = I + G$ breaks the symmetry between rows and columns ($col_j(g) \neq row_j(g)^T$).
* **Implemented Asymmetric SMW:** Used the explicit formula `g' = g - gU * M^{-1} * VTg_half`, where left side uses columns (`gU`) and right side uses rows (`VTg_half`). 
* Used explicit matrix generation (`Eigen::Matrix2cd`) which maps easily to BLAS `zgeru` (as demonstrated in `spinless_tV.h`) while remaining readable.

## 3. What's the Remaining Issue
* **Integration into Sweep:** The mixing operator is not updated during the right sweep in `src/pfqmc.cpp` (see `FIXME: update Mixing for the sign`).
* **Ratio Normalization/Sign:** The determinant ratio $1 + \sigma G_{jk}$ doesn't analytically match $\det(I+B')/\det(I+B)$ for the isolated 2x2 pure mixing block. Need to prove if it is globally valid inside the full trace.
* **Stabilization Performance:** `stabilizedLeftMultiply` falls back to dense matrix multiplication and needs optimization for sparse structures.

## 4. TODO
1. Integrate `MixingOperator` Metropolis sweep into `src/pfqmc.cpp`.
2. Optimize `MixingOperator::stabilizedLeftMultiply`.
3. Resolve theoretical ratio formulation for the operator $B = i\sigma\tau_y$.
4. Test full interacting model to confirm sign statistics.
