include("ising-magic-uniform.jl")
include("ising-magic.jl")
using CairoMakie
using BenchmarkTools
using Printf
using Statistics

function run_benchmarks()
    # ==========================================
    # Part 1: Convergence & Bias Check
    # ==========================================
    println("--- Part 1: Convergence Check (N=4) ---")
    N_check = 4
    exact_val = ising_magic(N=N_check)
    println("Exact Value (N=$N_check): $exact_val")
    
    sample_sizes = [1000, 10000, 100000]
    n_runs = 10
    
    means = Float64[]
    stds = Float64[]
    errors = Float64[]
    
    for samples in sample_sizes
        println("Running samples=$samples...")
        vals = Float64[]
        for i in 1:n_runs
            val = ising_magic_MCMC(N=N_check, samples=samples, thermalization=1000)
            push!(vals, val)
        end
        
        mean_val = mean(vals)
        std_val = std(vals)
        err = abs(mean_val - exact_val)
        
        push!(means, mean_val)
        push!(stds, std_val)
        push!(errors, err)
        
        @printf "Samples: %d, Mean: %.6f, Std: %.6f, Error: %.6f\n" samples mean_val std_val err
    end
    
    # Plot Convergence
    f1 = Figure()
    ax1 = Axis(f1[1, 1], xscale=log10, yscale=log10, 
        xlabel="Sample Size", ylabel="Absolute Error", title="MCMC Convergence (N=$N_check)")
    
    scatterlines!(ax1, sample_sizes, errors, label="Mean Error")
    
    # Handle log-scale error bars: clip lower bound to avoid negative values
    # errorbars!(x, y, low, high) takes offsets
    low_errs = [min(s, e * 0.999) for (e, s) in zip(errors, stds)]
    high_errs = stds
    errorbars!(ax1, sample_sizes, errors, low_errs, high_errs, color=:red)
    
    # Reference 1/sqrt(N) line
    ref_x = sample_sizes
    ref_y = errors[1] .* sqrt.(sample_sizes[1] ./ sample_sizes)
    lines!(ax1, ref_x, ref_y, linestyle=:dash, color=:gray, label="1/√N Scaling")
    
    axislegend(ax1)
    save("convergence_plot.png", f1)
    println("Saved convergence_plot.png")

    # ==========================================
    # Part 2: Speed Benchmark
    # ==========================================
    println("\n--- Part 2: Speed Benchmark ---")
    
    # Warm-up to exclude compilation time
    println("Warming up...")
    ising_magic(N=4)
    ising_magic_MCMC(N=4, samples=100, thermalization=10)
    
    Ns = [4, 6, 8, 10]
    times_exact = Float64[]
    times_mcmc = Float64[]
    
    mcmc_samples = 10000 # Fixed samples for speed comparison
    
    for N in Ns
        println("Benchmarking N=$N...")
        
        # Exact
        if N <= 8 # Exact gets very slow
            t_exact = @elapsed ising_magic(N=N)
            push!(times_exact, t_exact)
        else
            push!(times_exact, NaN)
        end
        
        # MCMC
        t_mcmc = @elapsed ising_magic_MCMC(N=N, samples=mcmc_samples, thermalization=1000)
        push!(times_mcmc, t_mcmc)
        
        @printf "N=%d: Exact=%.4fs, MCMC=%.4fs\n" N (N<=10 ? times_exact[end] : NaN) times_mcmc[end]
    end
    
    # Plot Benchmark
    f2 = Figure()
    ax2 = Axis(f2[1, 1], yscale=log10, 
        xlabel="System Size N", ylabel="Time (s)", title="Execution Time: Exact vs MCMC")
    
    valid_exact = .!isnan.(times_exact)
    scatterlines!(ax2, Ns[valid_exact], times_exact[valid_exact], label="Exact", color=:blue)
    scatterlines!(ax2, Ns, times_mcmc, label="MCMC ($mcmc_samples samples)", color=:orange)
    
    axislegend(ax2)
    save("benchmark_plot.png", f2)
    println("Saved benchmark_plot.png")
end

run_benchmarks()
