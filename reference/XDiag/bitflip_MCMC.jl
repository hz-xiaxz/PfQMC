using XDiag
using LinearAlgebra
using Random
using Statistics

# ==============================================================================
# 1. Physics & Math Helpers (Optimized)
# ==============================================================================

function get_ground_state(N::Int, h::Float64)
    # Setup Hamiltonian
    hs = Spinhalf(N)
    ops = OpSum()
    J = 1.0
    for i in 1:N
        s1 = i
        s2 = mod1(i + 1, N)
        ops += -4.0 * J * Op("SzSz", [s1, s2])
        ops += -1.0 * h * Op("S+", [s1])
        ops += -1.0 * h * Op("S-", [s1])
    end
    
    # Solve
    _, psi_state = eig0(ops, hs)
    return ComplexF64.(vector(psi_state))
end

# Fast Bitwise Weight Calculation (O(2^N) per call)
function compute_weight(psi::Vector{ComplexF64}, p_idx::Int, N::Int)::Float64
    dim = length(psi)
    x_mask = 0
    z_mask = 0
    temp = p_idx
    
    for k in 0:(N-1)
        code = temp & 3
        temp >>= 2
        if code == 1; x_mask |= (1 << k); end
        if code == 2; x_mask |= (1 << k); z_mask |= (1 << k); end
        if code == 3; z_mask |= (1 << k); end
    end

    expectation = 0.0 + 0.0im
    for col_i in 0:(dim-1)
        row_i = col_i ⊻ x_mask
        phase_bit = count_ones(col_i & z_mask)
        sign = (phase_bit & 1 == 1) ? -1.0 : 1.0
        expectation += conj(psi[col_i + 1]) * psi[row_i + 1] * sign
    end

    return abs2(expectation) # Returns w(P)
end 

# ==============================================================================
# 2. Standard MCMC (Single Thread, No Tempering)
# ==============================================================================

function run_standard_mcmc(psi::Vector{ComplexF64}, N::Int; 
                           n_samples::Int=100_000, 
                           n_burnin::Int=10_000)
    
    # Initialize Random Walker
    # We must find a starting point with w > 0 to avoid division errors
    current_p = rand(0:(4^N-1))
    current_w = compute_weight(psi, current_p, N)
    
    # If we started at 0, walk randomly until we hit something
    while current_w < 1e-20
        current_p = rand(0:(4^N-1))
        current_w = compute_weight(psi, current_p, N)
    end
    
    samples_w = Float64[]
    sizehint!(samples_w, n_samples)
    accepted_moves = 0
    
    total_steps = n_samples + n_burnin
    
    for step in 1:total_steps
        
        # --- PROPOSAL: Local Flip ---
        # Pick 1 site, change its operator to one of the other 3 options
        # site = rand(0:(N-1))
        # shift = 2 * site
        
        # old_op = (current_p >> shift) & 3
        # # (old + 1..3) % 4 ensures we pick a diff operator
        # new_op = (old_op + rand(1:3)) & 3 
        
        # mask = 3 << shift
        # proposed_p = (current_p & ~mask) | (new_op << shift)
        proposed_p = rand(0:(4^N-1)) 
        # --- CALCULATE WEIGHT ---
        proposed_w = compute_weight(psi, proposed_p, N)
        
        # --- METROPOLIS STEP ---
        # Target Distribution pi(P) ~ w(P)
        # Ratio = w(new) / w(old)
        
        ratio = proposed_w / current_w
        
        if rand() < ratio
            current_p = proposed_p
            current_w = proposed_w
            if step > n_burnin
                accepted_moves += 1
            end
        end
        
        # --- COLLECT SAMPLE ---
        if step > n_burnin
            # We average the weights to estimate Omega
            push!(samples_w, current_w)
        end
    end
    
    avg_w = mean(samples_w)
    m2 = -log(avg_w) / N
    acc_rate = accepted_moves / n_samples
    return m2, acc_rate
end

# ==============================================================================
# 3. Benchmark
# ==============================================================================

function run_benchmark(N::Int)
    println("\n=== Benchmark N=$N (Single Thread, Standard MCMC) ===")
    
    # 1. Physics
    psi = get_ground_state(N, 1.00)
    
    # 2. Exact
    println("\n[Exact Method]")
    t_exact = @elapsed begin
        sum_w2 = 0.0
        for p in 0:(4^N - 1)
            w = compute_weight(psi, p, N)
            sum_w2 += w^2
        end
        m2_exact = -log(sum_w2/2^N) / N
    end
    println("  M2:    $m2_exact")
    println("  Time:  $(round(t_exact, digits=4))s")
    
    # 3. MCMC
    n_samp = 2_000_000
    println("\n[MCMC Method]")
    println("  Samples: $n_samp")
    t_mcmc = @elapsed begin
        m2_mcmc, acc_rate = run_standard_mcmc(psi, N; n_samples=n_samp)
    end
    
    println("  M2:    $m2_mcmc")
    println("  Time:  $(round(t_mcmc, digits=4))s")
    println("  Diff:  $(m2_mcmc - m2_exact)")
    println("  Acceptance Rate: $(round(acc_rate * 100, digits=2))%")
end

# Run
run_benchmark(12)