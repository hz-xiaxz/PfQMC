include("ising-magic-uniform.jl")
include("ising-magic.jl")
using Printf

function compare_methods()
    N = 4
    println("Running Exact Calculation for N=$N...")
    exact_val = ising_magic(N=N)
    @printf "Exact Magic Density: %.6f\n" exact_val

    println("\nRunning MCMC Calculation for N=$N...")
    # Increase samples for better accuracy
    mcmc_val = ising_magic_MCMC(N=N, samples=10000, thermalization=1000)
    @printf "MCMC Magic Density: %.6f\n" mcmc_val
    
    
    diff = abs(exact_val - mcmc_val)
    @printf "\nDifference: %.6f\n" diff
end

compare_methods()
