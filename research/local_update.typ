#import "@preview/physica:0.9.7":*
#import "@preview/theorion:0.4.1":*
#show:show-theorion
#import cosmos.fancy:*
= "Local Update" issue for Pauli String sampling
We want to compute $S_("SRE") = - log Omega_N = - log ( sum_P 1/2^N mel(psi, P, psi)^4)$ where $rho= ketbra(psi)$
== MCMC scheme
Let $ omega_P = mel(psi, P, psi)^2 $

We can prove by the propetry of Pauli Group:
$
  sum_P omega_P = 2^N
$

Thus we have
$
  Omega_N = sum_P (omega_P/2^N) omega_P = 1/M sum_(i=1)^M omega_(i)
$
where $P_(i)$ are sampled from the distribution $pi(P) = omega_P/2^N$.

The update scheme is
$
  "acc" = min(1, omega_(P"new")/omega_(P"old"))
$

== "Local Update" issue
How to sample through Pauli strings e.g.
$
  P = X_1 Y_2 I_3 Z_4 ...
$

We Propose
$
  P' = Z_1 Y_2 I_3 Z_4 ... = -i Y_1 P
$
which locally update one pauli character

thus the update ratio is
$
  omega_(P"new")/omega_(P"old") = mel(psi, P', psi)^2/(mel(psi, P, psi)^2 )
$
// Then the problem is, what if
// $ P ket(psi) "is a Stabilizer State of" i Y "Operator?" $

== Hamming Distance
#definition[
  The minimal number of positions at which the corresponding symbols are different between two strings of equal length.
]

To avoid the "Local Update" issue, we can define a hamming distance $d_H(P, P')$ between two Pauli strings $P$ and $P'$.

Then we can propose updates with $d_H(P, P') >= 2$

So local update fails when the importance landscapes acts like:
$
  omega_P = cases(0 "some states", 1 quad P in A)
$

where $forall p, p' in A, d_H (p,p') >> 1$

== Magic state as linear combination of stabilizer states
An $n-$qubit magic state can be expressed as a linear combination of stabilizer states:
$
  ket(psi) = 1/sqrt(k) sum_(i=1)^k c_i ket(phi_i)
$
These $ket(phi_i)$ should not share a common stabilizer group and they are orthonormal. $abs(c_i)^2=1$.

Let's ask what components contribute to $mel(psi, P, psi)$:
$
  mel(psi, P, psi) 
  &= sum_(i,j) c_i^* c_j mel(phi_i, P, phi_j) \
  &= sum_i abs(c_i)^2 mel(phi_i, P, phi_i) + sum_(i != j) c_i^* c_j mel(phi_i, P, phi_j)
$
#let pm = $plus.minus$
1. when $i=j$ , $P$ is the stabilize operator of $ket(phi_i)$ , $mel(phi_i, P, phi_i) = 1$. 

