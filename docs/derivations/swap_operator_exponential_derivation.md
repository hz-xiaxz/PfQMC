# Derivation of Swap Operator from Exponential Form

## Definition
The swap operator is given by:
$$ \mathcal{S} = \exp\left( \frac{\pi}{4} \sigma_m \gamma_\alpha \gamma_\beta \right) $$

We want to find the matrix $S$ that represents this operator in the single-particle Majorana basis, defined by the transformation:
$$ \mathcal{S}^\dagger \Psi \mathcal{S} = S \Psi $$
where $\Psi = (\gamma_\alpha, \gamma_\beta)^T$.

## Derivation via Matrix Exponential

The operator is given by:
$$ \mathcal{S} = \exp\left( \frac{\pi}{4} \sigma_m \gamma_\alpha \gamma_\beta \right) $$

This operator generates a rotation in the subspace spanned by $\gamma_\alpha$ and $\gamma_\beta$.
The term $\frac{1}{2} \gamma_\alpha \gamma_\beta$ is the generator of rotations in the $\alpha\beta$ plane.
Specifically, for an operator $U(\phi) = \exp(\frac{\phi}{2} \gamma_\alpha \gamma_\beta)$, the vector transformation is a rotation by angle $\phi$:
$$ \begin{pmatrix} \gamma_\alpha' \\ \gamma_\beta' \end{pmatrix} = \begin{pmatrix} \cos\phi & \sin\phi \\ -\sin\phi & \cos\phi \end{pmatrix} \begin{pmatrix} \gamma_\alpha \\ \gamma_\beta \end{pmatrix} $$

Comparing our operator $\mathcal{S}$ with $U(\phi)$:
$$ \frac{\phi}{2} = \frac{\pi}{4} \sigma_m \implies \phi = \frac{\pi}{2} \sigma_m $$

Now we compute the matrix exponential for the rotation matrix $S$ with angle $\phi = \frac{\pi}{2} \sigma_m$:
$$ S = \exp\left( \phi \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} \right) = \begin{pmatrix} \cos(\frac{\pi}{2}\sigma_m) & \sin(\frac{\pi}{2}\sigma_m) \\ -\sin(\frac{\pi}{2}\sigma_m) & \cos(\frac{\pi}{2}\sigma_m) \end{pmatrix} $$

Using $\cos(\pm \pi/2) = 0$ and $\sin(\pm \pi/2) = \pm 1$:
$$ S = \begin{pmatrix} 0 & \sigma_m \\ -\sigma_m & 0 \end{pmatrix} = \sigma_m \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} $$

This confirms the matrix representation used in the code.
