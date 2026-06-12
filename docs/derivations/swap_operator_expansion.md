# Expansion of Swap Operator $S$

## Objective
Expand the operator $S = \exp\left( \frac{\pi}{4} \sigma_m \gamma_\alpha \gamma_\beta \right)$ to eliminate the exponential function.

## 1. Properties of the Exponent
Let $X = \sigma_m \gamma_\alpha \gamma_\beta$.
First, let's compute the square of $X$:
$$ X^2 = (\sigma_m \gamma_\alpha \gamma_\beta) (\sigma_m \gamma_\alpha \gamma_\beta) $$
Since $\sigma_m^2 = 1$:
$$ X^2 = \gamma_\alpha \gamma_\beta \gamma_\alpha \gamma_\beta $$
Using the anticommutation relation $\{\gamma_\alpha, \gamma_\beta\} = 0 \implies \gamma_\beta \gamma_\alpha = -\gamma_\alpha \gamma_\beta$:
$$ X^2 = \gamma_\alpha (-\gamma_\alpha \gamma_\beta) \gamma_\beta = -(\gamma_\alpha^2) (\gamma_\beta^2) $$
Since Majorana fermions satisfy $\gamma^2 = 1$:
$$ X^2 = -(1)(1) = -1 $$

Thus, the operator $X$ behaves like the imaginary unit $i$ in the Taylor series expansion.

## 2. Taylor Series Expansion
The exponential of an operator $A$ is defined as:
$$ e^{\theta A} = \sum_{k=0}^\infty \frac{(\theta A)^k}{k!} $$
Here, let $\theta = \frac{\pi}{4}$ and $A = X$.
$$ S = e^{\theta X} = \sum_{k=0}^\infty \frac{\theta^k X^k}{k!} $$

Separating into even and odd powers:
$$ S = \sum_{n=0}^\infty \frac{\theta^{2n} X^{2n}}{(2n)!} + \sum_{n=0}^\infty \frac{\theta^{2n+1} X^{2n+1}}{(2n+1)!} $$

Using $X^2 = -1$, we have:
*   $X^{2n} = (X^2)^n = (-1)^n$
*   $X^{2n+1} = X (X^2)^n = X (-1)^n$

Substituting these back:
$$ S = \left( \sum_{n=0}^\infty \frac{(-1)^n \theta^{2n}}{(2n)!} \right) I + \left( \sum_{n=0}^\infty \frac{(-1)^n \theta^{2n+1}}{(2n+1)!} \right) X $$

## 3. Identifying Trigonometric Functions
The series in the parentheses are the Taylor expansions for cosine and sine:
*   $\cos(\theta) = \sum_{n=0}^\infty \frac{(-1)^n \theta^{2n}}{(2n)!}$
*   $\sin(\theta) = \sum_{n=0}^\infty \frac{(-1)^n \theta^{2n+1}}{(2n+1)!}$

Thus:
$$ S = \cos(\theta) I + \sin(\theta) X $$

## 4. Final Result
Substitute $\theta = \frac{\pi}{4}$ and $X = \sigma_m \gamma_\alpha \gamma_\beta$:
*   $\cos(\frac{\pi}{4}) = \frac{1}{\sqrt{2}}$
*   $\sin(\frac{\pi}{4}) = \frac{1}{\sqrt{2}}$

$$ S = \frac{1}{\sqrt{2}} I + \frac{1}{\sqrt{2}} \sigma_m \gamma_\alpha \gamma_\beta $$

$$ S = \frac{1}{\sqrt{2}} (1 + \sigma_m \gamma_\alpha \gamma_\beta) $$

This is the expanded form of the swap operator with the exponential eliminated.
