# Derivation of Swap Operator Matrix Form

## Problem Statement
We are given the transformation rule for the swap operator $\mathcal{S}$ acting on a pair of Majorana fermions $\psi^{(\alpha)}$ and $\psi^{(\beta)}$:

$$ \mathcal{S}^\dagger \begin{pmatrix} \psi^{(\alpha)} \\ \psi^{(\beta)} \end{pmatrix} \mathcal{S} = \sigma_m \begin{pmatrix} \psi^{(\beta)} \\ -\psi^{(\alpha)} \end{pmatrix} $$

We want to find the matrix representation $S$ such that this transformation holds.

## Derivation

Let the vector of operators be $\Psi = \begin{pmatrix} \psi^{(\alpha)} \\ \psi^{(\beta)} \end{pmatrix}$.
The transformation is given by $\mathcal{S}^\dagger \Psi \mathcal{S} = S_{target} \Psi$, where $S_{target} = \sigma_m \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$.

However, in the single-particle basis (which is what the matrices in the code represent), the operator $\mathcal{S}$ acts as a rotation matrix $S$ on the vector space spanned by the fermions.
The relation is defined as:
$$ \mathcal{S}^\dagger \psi_i \mathcal{S} = \sum_j S_{ij} \psi_j $$

Substituting our specific case:
1.  **First Component:**
    $$ \mathcal{S}^\dagger \psi^{(\alpha)} \mathcal{S} = \sigma_m \psi^{(\beta)} $$
    This means the first row of $S$ must pick out the second component ($\beta$).
    $$ S_{11} = 0, \quad S_{12} = \sigma_m $$

2.  **Second Component:**
    $$ \mathcal{S}^\dagger \psi^{(\beta)} \mathcal{S} = -\sigma_m \psi^{(\alpha)} $$
    This means the second row of $S$ must pick out the first component ($\alpha$) with a minus sign.
    $$ S_{21} = -\sigma_m, \quad S_{22} = 0 $$

## Resulting Matrix

Combining these, the matrix $S$ is:

$$ S = \sigma_m \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} $$

### Case $\sigma_m = 1$ (Swapped)
$$ S = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} $$

### Case $\sigma_m = -1$ (Inverse Swapped)
$$ S = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} $$

## Verification

Let's verify $S \Psi$:
$$ \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} \begin{pmatrix} \psi^{(\alpha)} \\ \psi^{(\beta)} \end{pmatrix} = \begin{pmatrix} \psi^{(\beta)} \\ -\psi^{(\alpha)} \end{pmatrix} $$
This matches the target vector on the RHS of the original equation (for $\sigma_m=1$).

## Implementation in Code

### Left Multiply ($B = S \times A$)
Acting on rows of $A$:
*   Row 1 of $B$ = Row 2 of $A$
*   Row 2 of $B$ = -Row 1 of $A$

### Right Multiply ($B = A \times S$)
Acting on columns of $A$:
$$ [C_1, C_2] \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} = [-C_2, C_1] $$
*   Col 1 of $B$ = -Col 2 of $A$
*   Col 2 of $B$ = Col 1 of $A$

### new Verification
$$
S = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & \sigma_m \\ -\sigma_m & 1 \end{pmatrix}
$$
