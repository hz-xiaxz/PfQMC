#import "@preview/physica:0.9.7": *
#show link:underline
= ED reproduction of Ising Magic density
I tried to reproduce Yiming et al.'s results on ising magic density using exact diagonalization.

running `XDiag/ising-magic-uniform.jl` gives
```julia
julia> include("XDiag/ising-magic-uniform.jl")
julia> @time ising_magic_finite_T(;N=8, beta=8)
 61.457345 seconds (1.35 M allocations: 123.190 GiB, 5.88% gc time)
(0.2875921768304569, 0.32935638486262475, 0.041764208032167856)
```
which aligns my dot picking result `0.285...` (I used a website to extract their data points from the figure, just look at the data itself is enough).
#image("fig/yiming.png")

Their result is thermal, you can see the last term $S_2 equiv - log (tr (rho^2)) = 0.04$ which is fairly large.

also pushing $beta->oo$ gives
```julia
julia> @time ising_magic_finite_T(;N=8, beta=100)
 55.869051 seconds (1.35 M allocations: 123.190 GiB, 6.81% gc time)
(0.2641815725196151, 0.26418157321637864, 6.967635376987902e-10)
```

A good agreement with my zero temperature result `0.264...`
```julia
julia> @time ising_magic(;N=8)
 16.789665 seconds (9.51 M allocations: 185.189 MiB, 2.50% gc time)
0.2641815725636629
```

Also from the tree tensor network literature see the link #link("https://doi.org/10.1103/PRXQuantum.4.040317")[PRX Quantum 4, 040317 (2023)] 
#image("fig/ttn.png")

Yes they match, so TTN has better ground state approximation. Also calibrating my ED approach is right. The next step is the importance sampling ED, then next change to tV chain.