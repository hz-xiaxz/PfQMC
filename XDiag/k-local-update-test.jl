using XDiag
using LinearAlgebra
using Random
using Printf
using Statistics

# ============================================================
# Pauli application utilities (from ising-magic.jl)
# ============================================================

function apply_pauli_op(psi, site::Int, pauli_idx::Int)
    if pauli_idx == 1 # I
        return copy(psi)
    elseif pauli_idx == 2 # X = S+ + S-
        return apply(Op("S+", site), psi) + apply(Op("S-", site), psi)
    elseif pauli_idx == 3 # Y = -i(S+ - S-)
        v1 = apply(1.0im * Op("S+", site), psi)
        v2 = apply(1.0im * Op("S-", site), psi)
        return v1 - v2
    elseif pauli_idx == 4 # Z = 2Sz
        return apply(2.0 * Op("Sz", site), psi)
    else
        error("Invalid Pauli index")
    end
end

function apply_transition(current_state, site::Int, old_idx::Int, new_idx::Int)
    if old_idx == new_idx
        return current_state
    end

    # Lookup table for Pauli multiplication: T = P_new * P_old
    lookup = Dict{Tuple{Int,Int},Int}()
    for a in 1:4, b in 1:4
        if a == 1
            lookup[(a, b)] = b
        elseif b == 1
            lookup[(a, b)] = a
        elseif a == b
            lookup[(a, b)] = 1
        elseif (a, b) in [(2, 3), (3, 2)]
            lookup[(a, b)] = 4
        elseif (a, b) in [(2, 4), (4, 2)]
            lookup[(a, b)] = 3
        elseif (a, b) in [(3, 4), (4, 3)]
            lookup[(a, b)] = 2
        end
    end

    trans_idx = lookup[(old_idx, new_idx)]
    return apply_pauli_op(current_state, site, trans_idx)
end

function compute_state_from_indices(psi0, indices)
    state = psi0
    for (site, p_idx) in enumerate(indices)
        if p_idx != 1
            state = apply_pauli_op(state, site, p_idx)
        end
    end
    return state
end

# ============================================================
# Hamming distance computation
# ============================================================

"""
Compute Hamming distance between two Pauli index arrays.
"""
function hamming_distance(indices1::Vector{Int}, indices2::Vector{Int})
    @assert length(indices1) == length(indices2)
    return sum(indices1 .!= indices2)
end

# ============================================================
# k-local update function
# ============================================================

"""
Perform a k-site update on the current Pauli string configuration.
Returns: (new_indices, new_state, new_weight, accepted)
"""
function k_local_update!(
    current_P_indices::Vector{Int},
    current_state,
    current_weight::Float64,
    psi0,
    k::Int,
    N::Int
)
    # Select k distinct random sites
    sites = randperm(N)[1:k]

    # Store old indices
    old_indices = [current_P_indices[s] for s in sites]

    # Propose new indices
    new_indices = [rand(1:4) for _ in 1:k]

    # Check if any change is proposed
    any_change = any(old_indices .!= new_indices)

    if !any_change
        # No change proposed, automatically accepted
        return current_P_indices, current_state, current_weight, true, 0
    end

    # Apply incremental updates
    trial_state = current_state
    for (i, site) in enumerate(sites)
        trial_state = apply_transition(trial_state, site, old_indices[i], new_indices[i])
    end

    trial_overlap = dot(psi0, trial_state)
    trial_weight = abs2(trial_overlap)

    ratio = (current_weight < 1e-12) ? 1.0 : (trial_weight / current_weight)

    # Hamming distance of this proposal
    proposed_hamming = sum(old_indices .!= new_indices)

    # Metropolis acceptance
    if rand() < min(1.0, ratio)
        # Accept
        for (i, site) in enumerate(sites)
            current_P_indices[site] = new_indices[i]
        end
        return current_P_indices, trial_state, trial_weight, true, proposed_hamming
    else
        # Reject - return 0 for hamming since we don't count rejected moves
        return current_P_indices, current_state, current_weight, false, 0
    end
end

# ============================================================
# Main k-local update test
# ============================================================

"""
Run k-local update MCMC test for a given k value.
Returns statistics about acceptance rate and Hamming distances.
"""
function run_k_local_test(;
    N::Int=6,
    h::Float64=1.0,
    k::Int=1,
    samples::Int=2000,
    thermalization::Int=500
)
    # Build Hamiltonian and get ground state
    block = Spinhalf(N)
    ops = OpSum()
    J = 1.0

    for i in 1:N
        j = mod1(i + 1, N)
        ops += -4.0 * J * Op("SzSz", [i, j])
        ops += -1.0 * h * Op("S+", [i])
        ops += -1.0 * h * Op("S-", [i])
    end

    _, psi0 = eig0(ops, block)

    # Initialize Pauli string
    current_P_indices = rand(1:4, N)
    current_state = compute_state_from_indices(psi0, current_P_indices)
    current_overlap = dot(psi0, current_state)
    current_weight = abs2(current_overlap)

    # Thermalization
    for _ in 1:thermalization
        current_P_indices, current_state, current_weight, _, _ = k_local_update!(
            current_P_indices, current_state, current_weight, psi0, k, N
        )
    end

    # Sampling with statistics collection
    n_accepted = 0
    n_proposed = 0

    # Track Hamming distances for accepted moves only
    hamming_distances_consecutive = Int[]

    # Store previous accepted state for Hamming distance
    prev_indices = copy(current_P_indices)

    for _ in 1:samples
        n_proposed += 1

        current_P_indices, current_state, current_weight, accepted, accepted_hamming = k_local_update!(
            current_P_indices, current_state, current_weight, psi0, k, N
        )

        if accepted && accepted_hamming > 0
            n_accepted += 1
            # Only record Hamming distance for accepted non-trivial moves
            push!(hamming_distances_consecutive, accepted_hamming)
            prev_indices = copy(current_P_indices)
        elseif accepted && accepted_hamming == 0
            # Trivial acceptance (no change proposed)
            n_accepted += 1
        end
        # Rejected moves: don't count anything
    end

    # Compute statistics
    acceptance_rate = n_accepted / n_proposed

    # Average Hamming distance per accepted move (only counting accepted moves)
    avg_hamming_accepted = isempty(hamming_distances_consecutive) ? 0.0 : mean(hamming_distances_consecutive)
    std_hamming_accepted = isempty(hamming_distances_consecutive) ? 0.0 : std(hamming_distances_consecutive)

    # Number of non-trivial accepted moves
    n_nontrivial_accepted = length(hamming_distances_consecutive)

    return (
        k = k,
        acceptance_rate = acceptance_rate,
        avg_hamming_accepted = avg_hamming_accepted,
        std_hamming_accepted = std_hamming_accepted,
        n_nontrivial_accepted = n_nontrivial_accepted,
        n_proposed = n_proposed
    )
end

"""
Run comprehensive k-local update test from k=1 to k=max_k.
This tests how different locality of updates affects mixing efficiency.
"""
function k_local_update_test(;
    N::Int=6,
    h::Float64=1.0,
    max_k::Int=5,
    samples::Int=2000,
    thermalization::Int=500,
    n_runs::Int=3
)
    println("=" ^ 80)
    println("K-Local Update Test for Ising Model Magic Sampling")
    println("=" ^ 80)
    @printf "System size N = %d, Transverse field h = %.2f\n" N h
    @printf "Samples per run = %d, Thermalization = %d, Runs = %d\n" samples thermalization n_runs
    println("=" ^ 80)
    println()

    # Adjust max_k if larger than N
    max_k = min(max_k, N)

    # Storage for results
    results = []

    for k in 1:max_k
        @printf "Testing k = %d local updates...\n" k

        # Run multiple times and average
        avg_hammings = Float64[]
        n_accepted_list = Int[]
        acc_rates = Float64[]

        for run in 1:n_runs
            res = run_k_local_test(
                N=N, h=h, k=k, samples=samples, thermalization=thermalization
            )
            push!(avg_hammings, res.avg_hamming_accepted)
            push!(n_accepted_list, res.n_nontrivial_accepted)
            push!(acc_rates, res.acceptance_rate)
        end

        push!(results, (
            k = k,
            avg_hamming_mean = mean(avg_hammings),
            avg_hamming_std = std(avg_hammings),
            n_accepted_mean = mean(n_accepted_list),
            n_accepted_std = std(n_accepted_list),
            acceptance_rate_mean = mean(acc_rates),
            acceptance_rate_std = std(acc_rates)
        ))
    end

    # Print summary table
    println()
    println("=" ^ 80)
    println("RESULTS SUMMARY")
    println("=" ^ 80)
    println()

    # Header
    @printf "%-3s | %-12s | %-12s | %-14s | %-18s\n" "k" "Acc Rate" "Normalized" "Avg Hamming" "Num Accepted"
    println("-" ^ 70)

    # For normalization: assume uniform distribution over 4^N Pauli strings
    # A k-local update proposes changing k sites. Under uniform distribution,
    # acceptance rate would be 1 (detailed balance with flat distribution).
    #
    # But a better baseline: if the target distribution has some structure,
    # we expect acceptance ~ (typical weight ratio)^k for independent sites.
    #
    # Simple normalization: divide by (3/4)^k to account for the fact that
    # each site has 3/4 chance of actually changing (not proposing same value).
    # This gives us "acceptance rate per actual change attempted".
    #
    # Alternative: normalize by 1/k to compare "acceptance per site touched"

    for r in results
        # Normalize: acc_rate / (3/4)^k gives acceptance assuming all k sites change
        # This removes the trivial k-dependence from "more sites = harder to accept"
        baseline_change_prob = (3/4)^r.k
        normalized_acc = r.acceptance_rate_mean / baseline_change_prob

        @printf "%-3d | %10.4f | %10.4f | %12.4f | %8.1f ± %-6.1f\n" r.k r.acceptance_rate_mean normalized_acc r.avg_hamming_mean r.n_accepted_mean r.n_accepted_std
    end

    println()
    println("=" ^ 80)
    println("INTERPRETATION")
    println("=" ^ 80)
    println("""
    - Acc Rate: Raw acceptance rate (fraction of proposals accepted)
    - Normalized: Acc Rate / (3/4)^k — removes trivial k-dependence
      If this is constant across k, the distribution factorizes (sites independent).
      If it decreases with k, there are k-body correlations in the distribution.
      If it increases with k, larger moves are "finding better paths" through config space.
    - Avg Hamming: Average sites changed per accepted move (only counts accepted)
    - Num Accepted: Total non-trivial accepted moves out of $samples proposals
    """)

    return results
end

# ============================================================
# Additional diagnostic: Autocorrelation analysis
# ============================================================

"""
Compute autocorrelation of Hamming distances to check mixing.
"""
function analyze_mixing(;
    N::Int=6,
    h::Float64=1.0,
    k::Int=1,
    samples::Int=5000,
    thermalization::Int=1000
)
    # Build system
    block = Spinhalf(N)
    ops = OpSum()
    J = 1.0

    for i in 1:N
        j = mod1(i + 1, N)
        ops += -4.0 * J * Op("SzSz", [i, j])
        ops += -1.0 * h * Op("S+", [i])
        ops += -1.0 * h * Op("S-", [i])
    end

    _, psi0 = eig0(ops, block)

    # Initialize
    current_P_indices = rand(1:4, N)
    current_state = compute_state_from_indices(psi0, current_P_indices)
    current_weight = abs2(dot(psi0, current_state))

    # Thermalization
    for _ in 1:thermalization
        current_P_indices, current_state, current_weight, _, _ = k_local_update!(
            current_P_indices, current_state, current_weight, psi0, k, N
        )
    end

    # Record reference state
    reference_indices = copy(current_P_indices)

    # Track Hamming distance from reference over time
    hamming_from_reference = Int[]

    for _ in 1:samples
        current_P_indices, current_state, current_weight, _, _ = k_local_update!(
            current_P_indices, current_state, current_weight, psi0, k, N
        )
        push!(hamming_from_reference, hamming_distance(reference_indices, current_P_indices))
    end

    # Compute statistics
    mean_hamming = mean(hamming_from_reference)
    std_hamming = std(hamming_from_reference)
    max_hamming = maximum(hamming_from_reference)

    @printf "k=%d: Mean Hamming from reference = %.2f ± %.2f, Max = %d\n" k mean_hamming std_hamming max_hamming

    return hamming_from_reference
end

# ============================================================
# Parameter sweep over h (transverse field)
# ============================================================

"""
Sweep h parameter and measure normalized acceptance rate and Hamming distance for each k.
"""
function sweep_h_parameter(;
    N::Int=6,
    h_values::Vector{Float64}=[0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0],
    max_k::Int=5,
    samples::Int=3000,
    thermalization::Int=500,
    n_runs::Int=3
)
    println("=" ^ 90)
    println("Parameter Sweep: h from $(minimum(h_values)) to $(maximum(h_values))")
    println("=" ^ 90)
    @printf "N = %d, samples = %d, thermalization = %d, runs = %d\n" N samples thermalization n_runs
    println("=" ^ 90)
    println()

    # Store all results
    all_results = Dict{Float64, Vector{NamedTuple}}()

    for h in h_values
        @printf "h = %.2f: " h
        results_for_h = []

        for k in 1:max_k
            avg_hammings = Float64[]
            acc_rates = Float64[]

            for run in 1:n_runs
                res = run_k_local_test(N=N, h=h, k=k, samples=samples, thermalization=thermalization)
                push!(avg_hammings, res.avg_hamming_accepted)
                push!(acc_rates, res.acceptance_rate)
            end

            push!(results_for_h, (
                k = k,
                avg_hamming = mean(avg_hammings),
                acceptance_rate = mean(acc_rates),
                normalized_acc = mean(acc_rates) / (3/4)^k
            ))
            print("k=$k ")
        end
        println("done")

        all_results[h] = results_for_h
    end

    # Print summary table
    println()
    println("=" ^ 90)
    println("NORMALIZED ACCEPTANCE RATE BY h AND k")
    println("=" ^ 90)

    # Header
    print(@sprintf "%-6s |" "h")
    for k in 1:max_k
        print(@sprintf " k=%-6d |" k)
    end
    println()
    println("-" ^ (8 + 10 * max_k))

    for h in h_values
        print(@sprintf "%-6.2f |" h)
        for r in all_results[h]
            print(@sprintf " %7.4f |" r.normalized_acc)
        end
        println()
    end

    println()
    println("=" ^ 90)
    println("AVERAGE HAMMING DISTANCE (ACCEPTED) BY h AND k")
    println("=" ^ 90)

    # Header
    print(@sprintf "%-6s |" "h")
    for k in 1:max_k
        print(@sprintf " k=%-6d |" k)
    end
    println()
    println("-" ^ (8 + 10 * max_k))

    for h in h_values
        print(@sprintf "%-6.2f |" h)
        for r in all_results[h]
            print(@sprintf " %7.4f |" r.avg_hamming)
        end
        println()
    end

    return all_results
end

# ============================================================
# Run the test
# ============================================================

if abspath(PROGRAM_FILE) == @__FILE__
    # Run the main test with more samples to check scaling
    println("PART 1: Scaling test with more samples")
    println("=" ^ 80)
    results = k_local_update_test(N=6, h=1.0, max_k=5, samples=5000, thermalization=1000, n_runs=5)

    println()
    println("Additional mixing analysis (Hamming distance from reference state):")
    println("-" ^ 60)
    for k in 1:5
        analyze_mixing(N=6, h=1.0, k=k, samples=5000, thermalization=1000)
    end

    # Sweep h parameter around critical region (h=0 has zero magic, importance sampling fails)
    println()
    println()
    println("PART 2: Parameter sweep over h (around critical region h≈1)")
    sweep_results = sweep_h_parameter(
        N=6,
        h_values=[0.8, 0.9, 1.0, 1.1, 1.2],
        max_k=5,
        samples=3000,
        thermalization=500,
        n_runs=3
    )
end
