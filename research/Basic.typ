#import "@preview/physica:0.9.7": *
#set math.equation(numbering: "(1)")
Based on assaad

For simplicity, we focus on Hubbard model,
$
  H   &= H_t + H_U\
  H_U &= U sum_i (n_(i,↑) - 1/2)(n_(i,↓) - 1/2) \
$

Trotterization:
#let Δτ = $Δ τ$
$
  Z = tr (e^(-β H)) = tr (e^(-Δτ H_t) e^(-Δτ H_U))^M + O(Δτ^2) \
$

== HS decomposition:

Grounded as
#let phia = $phi.alt$
$
  integral_( -∞ )^oo dd(phia) e^(-(phia+A)^2/2) = sqrt(2 pi)
$
rewritten as
$
  e^(A^2/2) = 1/sqrt(2 pi) integral_( -∞ )^oo dd(phia) e^(-phia^2/2 -phia A)
$

rewrite $H_U$ as
#let upa = $arrow.t$
#let downa = $arrow.b$
$
  H_U &= U(n_upa - 1/2)(n_downa - 1/2) \
      &= U n_upa n_downa - U/2 (n_upa + n_downa) + U/4 \
      &= -U/2 (n_upa - n_downa)^2 + U/4 \
$

We want to compute $e^(-Delta tau H_U)$

Let's introduce an decomposition ansatz, by
$
  e^(-Delta tau H_U) = gamma sum_(s = ±1) e^(alpha s (n_upa - n_downa))
$
where $alpha$ and $gamma$ are parameters to be determined.

If we decompose the operator on $ket(0), ket(upa), ket(downa), ket(upa downa)$, we have
$
  H_U ket(0)         &= U/4 ket(0) \
  H_U ket(upa)       &= -U/4 ket(upa) \
  H_U ket(downa)     &= -U/4 ket(downa) \
  H_U ket(upa downa) &= U/4 ket(upa downa) \
$

then the $e^(-Delta tau H_U)$
$
  e^(-Delta tau H_U) ket(0)         &= e^(-Delta tau U/4) ket(0) = 2gamma ket(0)\
  e^(-Delta tau H_U) ket(upa)       &= e^(Delta tau U/4) ket(upa) = 2 gamma cosh(alpha) ket(upa)\
  e^(-Delta tau H_U) ket(downa)     &= e^(Delta tau U/4) ket(downa) = 2gamma cosh(alpha) ket(downa)\
  e^(-Delta tau H_U) ket(upa downa) &= e^(-Delta tau U/4) ket(upa downa) = 2 gamma ket(upa downa)\
$

Thus,
$
  gamma = 1/2 e^(-Delta tau U/4) \
  cosh(alpha) = e^(Delta tau U/2)
$

Another form of HS transformation is, without explicitly introducing $n_upa -n_downa$ which breaks $S U(2)$ symmetry
$
  e^(-Delta tau H_U) = tilde(gamma) sum_(s = ±1) e^(i tilde(alpha) s (n_upa + n_downa -1))
$
where $cos(tilde(alpha)) = e^(-Delta tau U/2) $ and $tilde(gamma) = 1/2 e^(Delta tau U/4)$

However, this approach introduces complex numbers in the simulation

In practice we use
$
  e^(Delta tau lambda A^2) = sum_(l = plus.minus 1, plus.minus 2) gamma(l) e^(sqrt(Delta tau lambda) eta(l) A) + cal(O)(Delta tau^4)
$
where
$
  gamma(plus.minus 1) = 1+sqrt(6)/3, gamma(plus.minus 2) = 1 - sqrt(6)/3 \
  eta(plus.minus 1) = plus.minus sqrt(2(3 - sqrt(6))), eta(plus.minus 2) = plus.minus sqrt(2(3 + sqrt(6)))
$

== AFQMC
Ok fine now, we can move back to the AFQMC.

#let trotU = $e^(-Δ τ H_U)$
We can now write $trotU$ as
$
  trotU = C sum_({s_i}) exp(alpha sum_i s_i (n_(i, upa) - n_(i, downa)))
$

denote
$
  H_t                                         &= c^dagger T c\
  alpha sum_i s_i (n_(i, upa) - n_(i, downa)) &= c^dagger V(s) c\
$
#let cd = $c^dagger$
$
  U_s(tau_2, tau_1) &= product_(n=n_1+1)^n_2 e^(cd V(s_n) c) e^(-Delta_tau cd T c)\
  B_s(tau_2, tau_1) &= product_(n=n_1+1)^n_2 e^(V(s_n)) e^(-Delta_tau T)\
$
where $n_1 Delta_tau = tau_1$ and $n_2 Delta_tau = tau_2$

== Slater Determinant and their properties

$
  tr(product_i e^(cd A_i c)) = det(1 + product_i e^(A_i)) \
$

Thus,
$
  Z = C^m sum_({s_i}) det(1 + B_s (beta, 0)) \
$
