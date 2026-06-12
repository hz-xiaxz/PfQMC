# Why the Mixing Operator Update Changes the Sign

## 1. General AFQMC Update Logic
In Auxiliary Field QMC, we sample auxiliary fields (or discrete operator choices) $\mathbf{x}$ with a probability distribution proportional to the absolute value of the weight:
$$ P(\mathbf{x}) \propto |W(\mathbf{x})| $$
where $W(\mathbf{x}) = \text{Tr}[\prod B(\mathbf{x}_i)]$ (or $\langle \psi_L | \prod B | \psi_R \rangle$).

When we propose a change $\mathbf{x} \to \mathbf{x}'$, we calculate the ratio of weights:
$$ r = \frac{W(\mathbf{x}')}{W(\mathbf{x})} $$

The Metropolis acceptance probability is:
$$ P_{acc} = \min\left(1, \left| \frac{W(\mathbf{x}')}{W(\mathbf{x})} \right| \right) = \min(1, |r|) $$

Since we sample based on $|W|$, but the physical weight is $W$, we must track the **phase** (or sign) of the weight separately.
$$ W_{new} = W_{old} \times r $$
$$ \text{Phase}_{new} = \text{Phase}_{old} \times \frac{r}{|r|} $$

This phase is used in the final measurement:
$$ \langle \mathcal{O} \rangle = \frac{\sum_i \mathcal{O}(\mathbf{x}_i) \sigma_i}{\sum_i \sigma_i} $$
where $\sigma_i$ is the accumulated phase of sample $i$.

## 2. Mixing Operator Specifics
The Mixing Operator introduces a discrete choice at time $\tau=0$:
*   **State 0**: Identity $I$ (no swap)
*   **State 1**: Swap $S$ (swap replicas)

When we propose a flip (e.g., $I \to S$), we are changing one operator in the chain.
The ratio $r$ is calculated using the Green's function $G$:
$$ r = \frac{\det(I + U_{new} \dots)}{\det(I + U_{old} \dots)} = \det(I + (1-G)\Delta) $$
(Note: For Majoranas, the ratio is the **Pfaffian**, which is the square root of this determinant).

### Why is there a sign change?
The ratio $r$ is not guaranteed to be positive real.
1.  **Complex Weights**: If the Hamiltonian or the swap operator involves complex numbers (or if the path integral leads to negative weights), $r$ can be negative or complex.
2.  **Swap Operator**: The swap operator $S = \frac{1}{\sqrt{2}}(I + \sigma J)$ is unitary but can introduce phases.
3.  **Fermionic Sign Problem**: Even with real matrices, the determinant/Pfaffian can be negative.

Therefore, whenever we accept a move, we must update the global phase by multiplying by $r/|r|$.

## 3. The Ambiguity in `sqrt(det)`
In the code, the ratio is calculated as:
```cpp
DataType det_K = K.determinant(); // Ratio squared
DataType r = std::sqrt(det_K);    // Ratio
```
The determinant gives $r^2$. Taking the square root introduces a sign ambiguity ($\pm r$).
*   If the physics is continuous, we might track the branch.
*   Here, since it's a discrete flip ($I \leftrightarrow S$), the sign of the square root is non-trivial to determine purely from `sqrt(det)`.
*   However, the magnitude $|r|$ is correct, so the sampling probability is correct.
*   The **error** in the current code is that `std::sqrt` arbitrarily chooses the principal root (positive real part). If the true physical ratio should be negative (e.g., crossing a node), this code might miss the sign change.

## Summary
We compute the sign because we are performing importance sampling with respect to the **modulus** of the wavefunction/partition function. The phase of the ratio $r$ must be accumulated to correct for this reweighting.
