include("ising-magic-uniform.jl")
include("ising-magic.jl")
using CairoMakie
using Printf

function run_phase_diagram()
    println("--- Phase Diagram Sweep (h=0.8-1.2) ---")
    
    # Parameters for quick run
    N = 8
    samples = 65536*5 
    thermalization = 5000
    
    hs = range(0.8, 1.2, step=0.05)
    magic_vals = Float64[]
    exact_vals = Float64[]
    
    println("System Size N=$N")
    println("Samples = $samples")
    
    for h_val in hs
        println("\nCalculating for h=$h_val...")
        
        # MCMC
        mcmc_val = ising_magic_MCMC(N=N, h=h_val, samples=samples, thermalization=thermalization)
        push!(magic_vals, mcmc_val)
        
        # Exact (for comparison since N=8 is fast enough)
        # exact_val = ising_magic(N=N, h=h_val)
        # push!(exact_vals, exact_val)
        
        # @printf "h=%.2f: MCMC=%.5f, Exact=%.5f, Diff=%.5f\n" h_val mcmc_val exact_val abs(mcmc_val - exact_val)
    end
    exact_vals= [0.18564, 0.21332,0.23897,0.25751,0.26418,0.25813,0.24303,0.22385,0.20412]
    # Plot
    f = Figure()
    ax = Axis(f[1, 1], xlabel="Transverse Field h", ylabel="Magic Density", 
        title="Ising Model Magic Density (N=$N)")
    
    scatterlines!(ax, hs, exact_vals, label="Exact", color=:blue)
    scatterlines!(ax, hs, magic_vals, label="MCMC-10000sweeps", color=:orange, linestyle=:dash)
    
    axislegend(ax)
    save("phase_diagram_MCMC.png", f)
    println("\nSaved phase_diagram.png")
end

run_phase_diagram()
