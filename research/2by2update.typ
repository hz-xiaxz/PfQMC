#import "@preview/physica:0.9.7":*
#import "@preview/theorion:0.4.1":*
#import cosmos.clouds:*
#show: show-theorion
= $2N$ by $2N$ block update

We need to update inside the replica group (12) and (34).

According to Appendix F in Han's paper;

== Hamiltonian part

We consider $t-V$ model,

which gives flipping a rank-4 skew symmetric matrix,
$
  m_((1j),(1k)) &= + 2i lambda sigma\
  m_((2j),(2k)) &= + 2i lambda sigma\
$

Consider
$
  B_m^(1) &= exp(-m^(1))\
  B_m^(2) &= exp(-m^(2))\
$

The total $B$ matrix is block diagonal
$
  B_m = diag(B_m^((1)), B_m^((2)))
$

$
  G'_l &= 2 [I + A_l B_m]^(-1)-I\
       &= 2 [(I+A_l + (I + A_l -I) (B_m -I))]^(-1)-I\
       &= [I + (I-G_l)/2 (B_m -I)]^(-1) (I+G_l) -I\
       &= [I-U(C^(-1)+V U)^(-1) V] (I+G_l) -I\
       &= G_l - Q^T O Q
$
here $C = B_m - I = diag(B_m^((1))-I, B_m^((2))-I)$

$
  U = mat(
    ..., -G^(00)_(1 j), ..., -G^(00)_(1 k), ..., -G^(01)_(1,j+N ), ..., -G^(01)_(1,k+N), ...;..., -G^(00)_(2 j), ..., -G^(00)_(2 k), ..., -G^(01)_(2,j+N ), ..., -G^(01)_(2,k+N), ...;dots.v, dots.v, dots.v, dots.v, dots.v, dots.v, dots.v, dots.v, dots.v;
  ) "/" 2
$

let's see if Q could be block diagonalized

Do this manually
$
  C^(-1) = diag((B_m^((1))-I)^(-1), (B_m^((2))-I)^(-1)) \
$

For tv model
$
  B_m - I = mat(cosh(2 lambda -1), -i sigma sinh(2 lambda);i sigma sinh(2 lambda), cosh(2 lambda) -1) \
$


$
  F= (B_m-I)^(-1) = mat(-1/2, -i sigma/2 coth(lambda);i sigma/2 coth(lambda), -1/2) \
$

$
  C^(-1) = diag(F, F) \
$

$
  V U = mat(
    1, - G^(00)_(j k), -G^(01)_(j, j+N), -G^(01)_(j, k+N);-G^(00)_(k j), 1, -G^(01)_(k, j+N), -G^(01)_(k, k+N);-G^(10)_(j+N, j), -G^(10)_(j+N, k), 1, -G^(11)_(j+N, k+N);-G^(10)_(k+N, j), -G^(10)_(k+N, k), -G^(11)_(k+N, j+N), 1
  ) "/" 2 \
$
since $G_(i i)= 0$

== Sequential update
Instead, we can update replica one first, that saving time, and no need to calculate the $4times 4$ matrix inverse again.
let
$
  B_m = diag(B_m^((1)), 0)
$

Thus
$
  U = mat(
    ..., -G^(00)_(1 j), ..., -G^(00)_(1 k), ;..., -G^(00)_(2 j), ..., -G^(00)_(2 k);dots.v, dots.v, dots.v, dots.v;..., 1-G^00_(j j), dots.v, -G^(00)_(j k);dots.v, dots.v, dots.v, dots.v;..., -G^(10)_(N+1, j), ..., ...
  ) "/" 2
$

I stil has concern about the analytical version of $O$, so derive now.

$
  V U = mat(1, -G^(00)_(j k);-G^(00)_(k j), 1)
$
Since the second replica is irrelevant in this first replica update, $G^(01)$ is not involved in the kernel

Let $g$ be the relevant block in $G_l$,
$
  g= mat(0, G^(00)_(j k);-G^(00)_(j k), 0)
$
Note that $G$ is skew-symmetric, so $G_(k j) = - G_(j k)$.

The update is then
$
  G' = G - Q^T O Q
$
where $Q = (1+G)_((j,k), (:))$.
$Q$ is a $2 times 2N$ matrix consisting of the $j$-th and $k$-th rows of $1+G$.

From the Woodbury identity on
$
  G' = (I + (I-G)/2 (B_m-I))^(-1) (I+G) - I
$
we can identify
$
  U = Q^T / 2
$
and
$
  V Q^T = I - g
$
where $g = mat(0, y;-y, 0)$ and $y = G_(j k)$.

The update matrix $O$ is given by
$
  O &= (2 C^(-1) + I - g)^(-1) \
    &= (mat(-1, -2z;2z, -1) + mat(1, -y;y, 1))^(-1) \
    &= mat(0, -2z-y;2z+y, 0)^(-1) \
    &= 1/(2z+y) mat(0, 1;-1, 0)
$
where $z = i sigma/2 coth(lambda)$.

This is a very compact result.

== Update of off-diagonal block $G^(01)$

Although $G^(01)$ does not enter the kernel $O$, it still gets updated through the $Q$ matrix.

The matrix $Q$ spans all replica columns. For indices $j, k$ in replica 0:
$
  Q = (1+G)_((j,k), (:)) = mat((1+G)^(00)_(j, :), G^(01)_(j, :);(1+G)^(00)_(k, :), G^(01)_(k, :)) = [Q_0 | Q_1]
$
where:
- $Q_0$ is $2 times N$: rows $j, k$ of $(1+G)^(00)$ (includes $+1$ on diagonal)
- $Q_1$ is $2 times N$: rows $j, k$ of $G^(01)$ (pure off-diagonal block)

The rank-2 update $Q^T O Q$ has block structure:
$
  (Q^T O Q)^((r s)) = Q_r^T O Q_s
$

Therefore, the off-diagonal block updates as:
$
  G'^(01) = G^(01) - Q_0^T O Q_1
$

Expanding with $O = 1/(2z+y) mat(0, 1;-1, 0)$:
$
  G'^(01)_(a b) = G^(01)_(a b) - 1/(2z+y) [(1+G)^(00)_(j a) G^(01)_(k b) - (1+G)^(00)_(k a) G^(01)_(j b)]
$

This can be written in matrix form. Let $q_j, q_k$ be rows $j, k$ of $(1+G)^(00)$ (as column vectors after transpose), and $r_j, r_k$ be rows $j, k$ of $G^(01)$:
$
  G'^(01) = G^(01) - 1/(2z+y) (q_j r_k^T - q_k r_j^T)
$

=== Key observation

The kernel $O$ only depends on $G^(00)_(j k)$ (the intra-replica Green's function element), but the off-diagonal block $G^(01)$ is still updated via the outer product structure of $Q^T O Q$.

This means a block-restricted implementation (updating only the diagonal block) would miss the coupling to off-diagonal blocks.

== update ratio
#let Pf = $"Pf"$
$
  R = (-1)^N eta_m "Pf" mat(mat(G_m, 0;0, 0), -I;I, mat(G^00, G^01;-G^01, G^11)))
$
So the acceptance ratio is used similarly.
