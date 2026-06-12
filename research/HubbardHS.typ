$
  H = -t sum_(i,σ) (c^†_(i,σ) c_(i+1,σ) + h.c.) + U sum_i n_(i,↑) n_(i,↓) - μ sum_(i,σ) n_(i,σ)
$

#let upa = $arrow.t$
#let downa = $arrow.b$
Majorana fermion representation:
$
  c_(σ)        &= (γ^1_(σ) + i γ^2_(σ)) / 2\
  c^dagger_(σ) &= (γ^1_(σ) - i γ^2_(σ)) / 2\
$
#let gamone = $gamma^1$
#let gamtwo = $gamma^2$

$
  c^dagger_(i, sigma) c_(i+1, sigma) &= 1/4 (gamone_(i, sigma) - i gamtwo_(i, sigma)) (gamone_(i+1, sigma) + i gamtwo_(i+1, sigma)) \
                                     &= 1/4 [gamone_(i, sigma) gamone_(i+1, sigma) + i gamone_(i, sigma) gamtwo_(i+1, sigma) - i gamtwo_(i, sigma) gamone_(i+1, sigma) + gamtwo_(i, sigma) gamtwo_(i+1, sigma)]\
$
It's Hermitian conjugate is:
$
  c^dagger_(i+1, sigma) c_(i, sigma) &= 1/4 [gamone_(i+1, sigma) gamone_(i, sigma) + i gamone_(i+1, sigma) gamtwo_(i, sigma) - i gamtwo_(i+1, sigma) gamone_(i, sigma) + gamtwo_(i+1, sigma) gamtwo_(i, sigma)]\
$

summation is
$
  c^dagger_(i, sigma) c_(i+1, sigma) + c^dagger_(i+1, sigma) c_(i, sigma) &= 1/2 [gamone_(i, sigma) gamone_(i+1, sigma) + gamtwo_(i, sigma) gamtwo_(i+1, sigma)]\
$

The kinetic term
$
  T = -t/2 i (gamma^1_(i, sigma) gamma^2_(i+1, sigma) + gamma^1_(i+1, sigma) gamma^2_(i, sigma))
$

HS transformation for the interaction term:
$
  U (n_(i,↑)-1/2) (n_(i,↓)-1/2) = -U/4 (gamone_(i, upa) gamtwo_(i, upa) gamone_(i, downa) gamtwo_(i, downa))

$

I will keep the Hamiltonian half-filled now on.
