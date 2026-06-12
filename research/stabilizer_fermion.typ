#import "@preview/physica:0.9.7":*
#import "@preview/theorion:0.4.1":*
#import cosmos.clouds: *
#set-inherited-levels(1)
#set-theorion-numbering("1.1")
#show: show-theorion

= Stabilizer state
We'd like to see under Pauli basis, what states have Stabilizer Renyi Entropy zero
#let pm = $plus.minus$
#let ox = $times.circle$
#definition(title: [Pauli Group])[
  $
    cal(P) = {pm 1, pm i} ox {I, X, Y, Z}
  $
]

#remark[
  Without $pm i$, this would never form a group, since $X Y = i Z $
]

#definition(title: [Stabilizer Group])[
  A Stabilizer Group $cal(S)$ is an *abelian* subgroup of the Pauli group $cal(P)_n$ that does not contain $-I$.

]

#definition(title: [Stabilizer Code])[
  A Stabilizer Code is defined by a Stabilizer Group $cal(S)$, such that
  $
    Q(cal(S)) ={ket(psi) | P ket(psi) = ket(psi) forall P in cal(S)}
  $
  Then $Q(cal(S))$ is the Stabilizer code and $cal(S)$ is its Stabilizer.
]

#definition(title: [Stabilizer State])[
  If a Stabilizer code of $N-$ qubit has $N$ independent generators, then the code space is one-dimensional, and the unique state in the code space is called a Stabilizer State.
]

#problem[How to exhaust all Stabilizer States of N qubits?]
