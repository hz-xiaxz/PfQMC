# Swap Operator Design & Implementation Plan

This document outlines the derivation and implementation details for the Swap Operator in the PfQMC algorithm, specifically for Renyi entropy calculations.

## 1. Update Matrix of the Swap Operator

**Goal**: Derive the update formula for switching between Identity ($I$) and Swap ($S$) at the boundary ($\tau=0$) without introducing an obvious sign problem.

### Derivation
Let the Green's function be $G = (1 + U)^{-1}$, where $U = T e^{-\int_0^\beta H d\tau}$.
When we insert a Swap operator $S$ at $\tau=0$, the propagator becomes $U' = U S$.
The new Green's function is $G' = (1 + U S)^{-1}$.

We can relate $G'$ to $G$:
$$
\begin{aligned}
G (G')^{-1} &= (1 + U)^{-1} (1 + U S) \\
&= (1 + U)^{-1} (1 + U + U(S - I)) \\
&= I + (1 + U)^{-1} U (S - I) \\
&= I + G (G^{-1} - I) (S - I) \\
&= I + (I - G) (S - I)
\end{aligned}
$$
Thus,
$$ G' = [I + (I - G)(S - I)]^{-1} G $$

Let $\Delta = S - I$. The update matrix is $M = I + (I - G)\Delta$.
Since $S$ only acts on the replicas being swapped (say replica $A$ and $B$), $\Delta$ is sparse. If we have $N_{rep}$ replicas and total dimension $D_{tot}$, $\Delta$ is non-zero only in the $2N \times 2N$ block corresponding to $A$ and $B$ (where $N$ is the single-replica dimension).

### Ratio Calculation
The acceptance ratio for inserting $S$ is:
$$ r = \frac{\det(1 + U S)}{\det(1 + U)} = \det(G (1 + U S)) = \det(I + (I - G)(S - I)) $$
This reduces to the determinant of the $2N \times 2N$ block:
$$ r = \det( I_{2N} + [(I - G)(S - I)]_{2N \times 2N} ) $$

### Sign Problem
The Swap operator $S$ for fermions can have a sign. In the Majorana representation, $S$ is a permutation matrix.
For a single site with Majoranas $\gamma_1, \gamma_2$, swapping two copies involves:
$\gamma_{A,1} \to \gamma_{B,1}, \gamma_{A,2} \to \gamma_{B,2}$ and vice versa.
The determinant of the permutation matrix $S$ depends on the dimension.
However, the ratio $r$ calculated above includes the physics of the transition. If the sign of $r$ is complex or negative, it contributes to the phase.
The user notes "no obvious sign problem", which implies we should ensure the formulation doesn't introduce *artificial* sign problems (e.g. from bad basis choices). Using the determinant formula above is standard and robust.

## 2. Simultaneous 2N * 2N Block Update

**Goal**: Efficiently update $G$ using the Woodbury matrix identity for the $2N \times 2N$ block.

### Formula
We need to compute $G' = (I + X)^{-1} G$ where $X = (I - G)\Delta$.
Using the Woodbury identity $(I + UV)^{-1} = I - U(I + VU)^{-1}V$.
Here, $\Delta$ is sparse, but $(I-G)\Delta$ is not necessarily sparse in the rows.
However, $\Delta$ can be factorized or we can just restrict to the subspace.
Let $P$ be the projection onto the $2N$ indices of replicas $A, B$.
$\Delta = P^T \delta P$, where $\delta$ is the $2N \times 2N$ block of $S-I$.
Then $X = (I - G) P^T \delta P$.
Let $U = (I - G) P^T$ and $V = \delta P$.
$$ G' = G - U (I_{2N} + V U)^{-1} V G $$
$$ G' = G - [(I - G) P^T] (I_{2N} + \delta P (I - G) P^T)^{-1} \delta P G $$

Let $K = I_{2N} + \delta P (I - G) P^T$. This is a $2N \times 2N$ matrix.
Note that $P (I - G) P^T = I_{2N} - G_{sub}$, where $G_{sub}$ is the $2N \times 2N$ submatrix of $G$.
So $K = I_{2N} + \delta (I_{2N} - G_{sub})$.
Since $\delta = S_{sub} - I_{2N}$,
$K = I + (S_{sub} - I)(I - G_{sub}) = I + S_{sub} - S_{sub} G_{sub} - I + G_{sub} = S_{sub} + G_{sub} - S_{sub} G_{sub} = S_{sub} (I - G_{sub}) + G_{sub}$.
Wait, simpler:
$K = I + (S-I)(I-G) = I + S - SI - I + G = S - SG + G$ (subscript implied).
Actually, let's stick to $K = I + \delta (I - G_{sub})$.

**Algorithm**:
1. Extract $G_{sub}$ (the $2N \times 2N$ block of $G$ for replicas $A, B$).
2. Compute $\delta = S_{sub} - I_{2N}$.
3. Compute kernel $K = I_{2N} + \delta (I_{2N} - G_{sub})$.
4. Compute inverse $K^{-1}$.
5. Construct update matrices:
   $L = (I - G) P^T$ (The columns of $I-G$ corresponding to $A, B$).
   $R = \delta P G$ (The rows of $G$ multiplied by $\delta$).
6. $G' = G - L K^{-1} R$.

This updates all elements of $G$ simultaneously.

## 3. Boundary Application and Swap Operator Design

**Goal**: Implement `MixingOperator` (or `SwapOperator`) class.

### Design
The `MixingOperator` should be a subclass of `Operator`.
It sits at index 0 (or end) of the operator list.
It manages the state of swaps between replica pairs.

*   **State**: A vector `is_swapped` of size `nReplicas / 2`.
*   **Matrix**: It effectively represents a matrix $M_{mix}$.
    *   If `is_swapped[p]` is false, block is $I$.
    *   If `is_swapped[p]` is true, block is $S$.
*   **Update**:
    *   Propose flipping `is_swapped[p]`.
    *   Calculate ratio $r$.
    *   If accepted, perform the $2N \times 2N$ update derived above.

### Boundary Condition
In `PfQMC`, the `mixingOp` is applied at the boundary.
The code in `pfqmc.cpp` already has placeholders:
```cpp
// Apply mixing operator at τ=0 first if it exists
if (mixingOp != nullptr) {
    mixingOp->left_propagate(g, tmp);
    // ...
}
```
We need to ensure `MixingOperator` implements `update(g)` which performs the Metropolis step for the swaps.

## 4. Better Sign Calculation with Local Update

**Goal**: Efficiently update the global sign.

### Sign Update
The global sign is updated by multiplying by the sign of the ratio $r$.
$Sign_{new} = Sign_{old} \times \text{sgn}(r)$.
$r = \det(K)$.
So we just need the sign of the determinant of the $2N \times 2N$ kernel matrix $K$.
Since $2N$ is small (relative to total size), this is fast.

### Local vs Global
The `update` method in `Operator` returns `DataType` which is the ratio (or sign of ratio?).
In `PfQMC::rightSweep`:
```cpp
signCur = op_array[l]->update(g);
this->sign *= signCur;
```
If `MixingOperator::update` performs the Metropolis step and updates $G$, it should return the sign change $\frac{W_{new}}{W_{old}} / |\frac{W_{new}}{W_{old}}|$.
Actually, usually `update` returns the ratio sign.

We need to make sure the sign of the determinant is computed correctly.
For complex Hermitian matrices, determinant is real.
For general matrices, it's complex.
$r$ is complex in general.
The "sign" usually refers to $r / |r|$.

### Implementation Detail
In `MixingOperator::update(g)`:
1. Iterate over pairs.
2. For each pair, propose flip.
3. Compute $K$ and $\det(K)$.
4. Metropolis check.
5. If accept:
   Update $G$.
   Accumulate sign change.
6. Return total sign change.
