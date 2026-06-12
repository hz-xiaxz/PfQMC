#import "@preview/physica:0.9.7": *

#set text(font: "New Computer Modern", size: 11pt)
#set par(justify: true)

// --- Title ---
#align(center)[
  #text(size: 17pt, weight: "bold")[Mathematical Deduction: Magic Density Estimation]
  \
  #text(size: 12pt, style: "italic")[Derivation using Monte Carlo Sampling]
]

#v(1em)

= 1. Problem Definition

We aim to estimate the **Magic Density** ($M_2$) for an $N$-qubit state $ket(psi)$. The theoretical definition is:

$ M_2 = - 1/N log(Omega) $

Where $Omega$ is defined as the average over Pauli strings $P in cal(P)_N$:

$ Omega = 1/2^N sum_(P in cal(P)_N) abs(mel(psi, P, psi))^4 $

For notational convenience, let us define the "weight" $w(P)$ as the squared expectation value:

$ w(P) = abs(mel(psi, P, psi))^2 $

Thus, our target quantity simplifies to:
$ Omega = 1/2^N sum_(P) w(P)^2 $

= 2. Sampling Distribution

We utilize a Monte Carlo method by sampling Pauli strings $P$ from a probability distribution $pi(P)$ proportional to their weights.

$ pi(P) = w(P) / Z_"sampling" $

where $Z_"sampling" = sum_P w(P)$ is the normalization constant.

== Proof of Normalization ($Z_"sampling" = 2^N$)
To find $Z_"sampling"$, we assume the density matrix is pure, $rho = ket(psi)bra(psi)$.

1. Expand the weight definition:
  $ w(P) = mel(psi, P, psi) mel(psi, P^dagger, psi) $

  Substituting $rho$ into the expression:
  $ w(P) = mel(psi, P rho P^dagger, psi) $

2. Sum over all Pauli strings ($P = P^dagger$ for Paulis):
  $ Z_"sampling" = sum_P w(P) = angle.l psi | (sum_P P rho P) | psi angle.r $

3. Apply the **Pauli Twirl Identity**:
  Averaging a state over the full Pauli group results in the maximally mixed state scaled by the dimension $d^2 = 4^N$:
  $ sum_P P rho P = 2^N bb(I) $

4. Calculate the final expectation:
  $ Z_"sampling" = mel(psi, 2^N bb(I), psi) = 2^N braket(psi, psi) = 2^N $

Thus, the specific sampling probability is:
$ pi(P) = w(P) / 2^N $

= 3. Derivation of the Estimator

We want to estimate $sum_P w(P)^2$. We rewrite this as an expectation value over $pi(P)$:

$ sum_P w(P)^2 &= sum_P [ w(P)/2^N dot 2^N w(P) ] \
             &= sum_P [ pi(P) dot (2^N w(P)) ] \
             &= bb(E)_pi [ 2^N w(P) ] $

Using Monte Carlo approximation with $M$ samples ${P_1, ..., P_M}$ drawn from $pi(P)$:

$ sum_P w(P)^2 approx 1/M sum_(i=1)^M 2^N w(P_i) $

= 4. Final Calculation of $Omega$

Substitute the approximation back into the definition of $Omega$:

$ Omega &= 1/2^N sum_P w(P)^2 \
      &approx 1/2^N ( 1/M sum_(i=1)^M 2^N w(P_i) ) $

The factor $2^N$ cancels out:

$ Omega approx 1/M sum_(i=1)^M w(P_i) $

= 5. Conclusion

The quantity inside the logarithm is the arithmetic mean of the weights of the sampled states.

#align(center)[
  #rect(inset: 12pt, radius: 4pt, stroke: 1pt + luma(100), fill: luma(245))[
    $ M_2 approx - 1/N log ( 1/M sum_(i=1)^M w(P_i) ) $
  ]
]

Where:
- $w(P_i) = abs(mel(psi, P_i, psi))^2$
- $M$ is the number of Monte Carlo samples.