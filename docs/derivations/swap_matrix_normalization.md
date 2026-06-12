# Normalization of the Swap Operator Matrix

## Question
Does the factor $\frac{1}{\sqrt{2}}$ in the operator expansion $S = \frac{1}{\sqrt{2}}(1 + \sigma_m \gamma_\alpha \gamma_\beta)$ imply a factor of $\frac{1}{\sqrt{2}}$ in the matrix representation used in `left_multiply`?

## Derivation

1.  **Operator Expansion**:
    $$ S = \frac{1}{\sqrt{2}} (1 + X) $$
    where $X = \sigma_m \gamma_\alpha \gamma_\beta$ and $X^\dagger = -X$, $X^2 = -1$.
    Note that $S$ is unitary:
    $$ S^\dagger S = \frac{1}{2} (1 - X)(1 + X) = \frac{1}{2} (1 + X - X - X^2) = \frac{1}{2} (1 - (-1)) = 1 $$

2.  **Adjoint Action (Matrix Representation)**:
    The matrix $M$ used in `left_multiply` is defined by the transformation of the field operators:
    $$ \Psi' = S^\dagger \Psi S $$
    
    Let's compute this explicitly for $\gamma_\alpha$:
    $$ S^\dagger \gamma_\alpha S = \frac{1}{2} (1 - X) \gamma_\alpha (1 + X) $$
    $$ = \frac{1}{2} (1 - X) (\gamma_\alpha + \gamma_\alpha X) $$
    
    Using commutation relations:
    *   $X \gamma_\alpha = \sigma_m \gamma_\alpha \gamma_\beta \gamma_\alpha = -\sigma_m \gamma_\beta$
    *   $\gamma_\alpha X = \sigma_m \gamma_\alpha \gamma_\alpha \gamma_\beta = \sigma_m \gamma_\beta$
    
    Substitute back:
    $$ = \frac{1}{2} (1 - X) (\gamma_\alpha + \sigma_m \gamma_\beta) $$
    $$ = \frac{1}{2} (\gamma_\alpha + \sigma_m \gamma_\beta - X \gamma_\alpha - X \sigma_m \gamma_\beta) $$
    $$ = \frac{1}{2} (\gamma_\alpha + \sigma_m \gamma_\beta - (-\sigma_m \gamma_\beta) - \sigma_m (\sigma_m \gamma_\alpha)) $$
    $$ = \frac{1}{2} (\gamma_\alpha + \sigma_m \gamma_\beta + \sigma_m \gamma_\beta - \gamma_\alpha) $$
    $$ = \sigma_m \gamma_\beta $$

    Similarly for $\gamma_\beta$:
    $$ S^\dagger \gamma_\beta S = -\sigma_m \gamma_\alpha $$

3.  **Conclusion**:
    The transformation is:
    $$ \begin{pmatrix} \gamma_\alpha' \\ \gamma_\beta' \end{pmatrix} = \begin{pmatrix} 0 & \sigma_m \\ -\sigma_m & 0 \end{pmatrix} \begin{pmatrix} \gamma_\alpha \\ \gamma_\beta \end{pmatrix} $$
    
    The factor of $\frac{1}{\sqrt{2}}$ from $S$ and $\frac{1}{\sqrt{2}}$ from $S^\dagger$ combine to form $\frac{1}{2}$, which exactly cancels the doubling of terms in the expansion.
    **Therefore, the matrix representation has elements $0, \pm 1$ and NO factor of $\frac{1}{\sqrt{2}}$.**
