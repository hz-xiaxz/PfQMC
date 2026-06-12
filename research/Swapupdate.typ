#import "@preview/physica:0.9.7":*
#import "@preview/theorion:0.4.1":*
#import cosmos.clouds:*
#show: show-theorion
= Swap Operator update
The total swap operator is
#let pm = $plus.minus$
$
  sum_( sigma_m = pm 1 ) S_m^(12)[sigma_m] * S_m^(34)[sigma_m] = 1+gamma_m^1 gamma_m^2 gamma_m^3 gamma_m^4
$

$
  S_m^(12)[sigma_m]
    &= e^(pi/4 sigma_m gamma_m^1 gamma_m^2) \
    &= exp(pi/8 sigma_m gamma_m^T mat(0, 1;-1, 0) gamma_m)\
    &= exp(-1/4 gamma_m^T h_(12) gamma_m) \
$

where
$
  h_(12) = - i pi/2 sigma_m tau_y
$

While same approach should apply for $S_m^(34)$, however, with a sign flip we will get a trivially sign problem free update.
$
  S_m^(34)[sigma_m]
    &= e^(-pi/4 sigma_m gamma_m^4 gamma_m^3) \
$

under the basis of
$ gamma_m = vec(gamma_m^4, gamma_m^3) $

we have $ h_(43) = - h_(12) $

So with total majorana vector
$
  Gamma_m = vec(gamma_m^1, gamma_m^2, gamma_m^4, gamma_m^3)
$

#remark[
  This tweak can also be viewed as changing the $1+gamma^1 gamma^2 gamma^3 gamma^4$ Majorana Strings to $1-gamma^1 gamma^2 gamma^3 gamma^4 = 1 + gamma^1 gamma^2 gamma^4 gamma^3$ , without physical consequences.
]

we have

$
  S_m^(12)[sigma_m] * S_m^(34)[sigma_m] = exp(-1/4 Gamma_m^T H_m Gamma_m) \
$

where
$
  H = -i pi/2 sigma_m tau_y times.o sigma_z
$

Let's check if this is sign problem free

Take
$
  T^+ &= I times.o sigma_z K\
  T^- &= tau_x times.o i sigma_y K
$

thus this falls into the majorana class of sign-free problem.

== update
$
  B = e^(-H) = exp(i pi/2 sigma_m tau_y times.o sigma_z) \
$
as
$
  (i tau_y times.o sigma_z)^2 = -1
$

we can expand the exponential as
$
  B &= sin(pi/2 sigma_m) (i tau_y times.o sigma_z) + cos(pi/2 sigma_m) I = i sigma_m tau_y times.o sigma_z\
$

$
  G_m &= 2(1+B)^(-1) -1 \
      &= 2 (I + i sigma_m tau_y times.o sigma_z)^(-1) -1 \
      & =2 mat(1, sigma_m, 0, 0;-sigma_m, 1, 0, 0;0, 0, 1, -sigma_m;0, 0, sigma_m, 1)^(-1) -1 \
      & = mat(1, -sigma_m, 0, 0;sigma_m, 1, 0, 0;0, 0, 1, sigma_m;0, 0, -sigma_m, 1) -1 \
      &= -i sigma_m tau_y times.o sigma_z \
$

$
  sinh(H/4) &= sinh(-i pi/8 sigma_m tau_y times.o sigma_z) \
            &= sigma_m sin(pi/8) (-i tau_y times.o sigma_z) \
$

#let Pf = $"Pf"$
#let IfN = $I_(4N times 4N)$

#definition[
  $
    eta_m = (-1)^N Pf mat(sqrt(2) sinh(H/4), -IfN;IfN, sqrt(2) sinh(H/4)) \
  $
]

let
$
  A = mat(0, -IfN;IfN, sqrt(2) sinh(H/4))
$

$
  Pf A = sqrt(det(A)) = 1
$
hint: this is only valid with dimension of even number $4N$

$
  A^(-1) = mat(sqrt(2) sinh(H/4), IfN;-IfN, 0)
$

$
  C      &= sqrt(2) sinh(H/4) \
  C^(-1) &= 1/(sigma_m sin(pi/8)) (i tau_y times.o sigma_z) \
$

To calculate $Pf C^(-1)$

We use
#theorem[
  $
    A = mat(0, B;-B^T, 0)
  $

  gives
  $
    Pf A = (-1)^(n(n-1)/2) det(B)
  $
]

and
$
  Pf (lambda M) = lambda^n Pf M
$
where $M$ is $2n times 2n$ matrix
$
  Pf C_m^(-1) = 1/(sqrt(2)sigma_m sin(pi/8))^2 Pf mat(0, sigma_z;-sigma_z, 0) = 1/(2 sin^2(pi/8)) (-1)(-1) = 1/(2 sin^2(pi/8)) \
$
thus
$
  Pf C^(-1) = (1/(2 sin^2(pi/8)))^N \
$

$
  eta &= (-1)^N (Pf (A+B C B^T)) /( Pf A) = (Pf ( C^(-1) + B^T A^(-1) B )) /(Pf(C^(-1)))\
      &= (-1)^N product_i^m eta_m
$
$
  eta_m = 1/(2 sin^2(pi/8))
$

== Swap Operator
$
  S_m^(12)[sigma_m]
    &= 1/sqrt(2) (1+ sigma_m gamma_m^1 gamma_m^2) \
$

We know this
$
  S_m^(12)[sigma_m]^dagger
$

$
  S^(12)_m [sigma_m] vec(gamma^((alpha))_m, gamma_m^((beta))) = 1/(sqrt(2)) mat(1, sigma_m;-sigma_m, 1) vec(gamma^((alpha))_m, gamma_m^((beta))) = 1/sqrt(2) (vec(gamma^((alpha))_m, gamma^((beta))_m)+ sigma_m vec(gamma^((beta))_m, -gamma^((alpha))_m)) \
$

How does $S_m^12$ acts on sector of the total $G$ matrix?

$
  S_m^12 [sigma_m] mat(0, G^(alpha beta)_(m , m+N);G^(beta alpha)_(m+N, m), 0) =1/(sqrt(2)) mat(-sigma_m G, G;-G, sigma_m G) \
$

Looks like the diagonal terms are no-longer zero! Does that make any sense?

== Local update
Copying from the `2by2update.typ` file, we have the kernel

$ 
  C^(-1) = (B_("mix")-I)^(-1) = mat(-1, -sigma_m; sigma_m ,-1)^(-1) = 1/2 mat(-1, sigma_m;-sigma_m, -1) \
 $

$
  O &= (2 C^(-1) + I - g)^(-1) \
    &= (2 * 1/2mat(-1, sigma_m;-sigma_m, -1) + mat(1, -G^12;-G^21, 1))^(-1)\
    &=mat(0, sigma_m -G^12;- sigma_m - G^21, 0)^(-1) \
    &= 1/(sigma_m - G^12 ) mat(0, -1;1, 0) \ 
$

Yes but is that true?

=== Rederive
$
  G'
    &= 2[I+A_l B_m]^(-1) -I\
    &= 2[(I+A_l) + (I+A_l -I) (B_m-I)]^(-1) -I\
$
As $G_l = 2(I+A_l)^{-1} -I$, we have
$
  G'
    &=[I + A_l/(I+A_l) (B_m -I)]^(-1) (I+G_l) -I\
$
$
  A_l/(I+A_l) = I - (I+A_l)^(-1) = I - (I+G_l)/2 = (I - G_l)/2 \
$

thus
$
  G' &= [I + (I-G_l)/2 (B_m -I)]^(-1) (I+G_l) -I\
$

=== Another update
Consider $Delta B = U V^T$

$
  G' = 2(I+B+U V^T)^(-1)
$

Let $X = I+B$, $A = X^(-1)$

#theorem(title: "Woodbury Identity")[
  $
    (X + U V^T)^(-1) = X^(-1) - X^(-1) U (I + V^T X^(-1) U)^(-1) V^T X^(-1)
  $
]

Thus
$
  G' = 2 [A- A U (I + V^T A U)^(-1) V^T A] \
$

