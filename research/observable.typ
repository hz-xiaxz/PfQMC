#import "@preview/physica:0.9.7": *
#import "@preview/theorion:0.4.1": *
#import cosmos.clouds: *
#show: show-theorion
Ref: Fakher 10.6.3.2
We need to compute the ising swap operator

For FTQMC
#definition[
  $
    U_s (tau_2, tau_1) &= product_(m=n_1+1)^n_2 e^(c^dagger V(s_n) c) e^(-Delta_tau c^dagger T c)\
    B_s (tau_2, tau_1) &= product_(...) e^(V(s_n))e^(-Delta_tau T)
  $
]

$
  tr(e^(-beta H) O)/tr(e^(-beta H)) = sum_s P_s expval(O)_s
$
where
$
  P_s         &= "normalize" (det(1+B_s (beta,0)))\
  expval(O)_s &= tr(U_s (beta, tau) O U_s (tau,0))/tr(U_s (beta,0))
$

However, our $O$ can't be written as the bilinear form in Fakher's paper $c^dagger A c$. Need to figure it out.

Pengfei's paper defines
#definition[
  $ S_m^(alpha beta)[sigma_m] ^dagger vec(psi_m^((alpha)), psi_m^((beta))) S_m^(alpha beta ) [sigma_m] = sigma_m vec(psi_m^((beta)), -psi_m^((alpha))) $
]

#problem[
  Can $U_s(beta, tau), U_s(tau,0)$ be extracted in this program?
]

Let's assume now that we have total $U$ with replica structure
$
  U = diag(U_1, U_2, U_3, U_4)
$

here we have
$
  O = S_m^(alpha beta) = e^(pi/2 sigma_m psi_m^((alpha) )psi_m^((beta)))
$

Ah shit, there is a 2 ambiguity in Pengfei's notation, now I denote
$
  gamma = sqrt(2) psi
$
thus we have
$
  {gamma_m,gamma_n} = 2 delta_(m n)
$

Thus
$ 
  O = exp(pi/4 sigma_m gamma^alpha_m gamma^beta_m)= 1/sqrt(2) (1+ sigma_m gamma^alpha_m gamma^beta_m)
 $
  

we follow fakher's approach
$
  expval(O)_s 
  &= (diff )/(diff eta) ln tr[U_s (beta, tau) e^(eta O) U_s (tau,0) ]|_(eta=0)\
  &= (diff )/(diff eta) ln det[1+B_s (beta, tau) e^(eta /sqrt(2)) e^(eta/sqrt(2) sigma_m) U_s (tau,0) ]|_(eta=0)\
  
$