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

    # Metropolis acceptance
    if rand() < min(1.0, ratio)
        # Accept
        for (i, site) in enumerate(sites)
            current_P_indices[site] = new_indices[i]
        end
        proposed_hamming = sum(old_indices .!= new_indices)
        return current_P_indices, trial_state, trial_weight, true, proposed_hamming
    else
        # Reject
        proposed_hamming = sum(old_indices .!= new_indices)
        return current_P_indices, current_state, current_weight, false, proposed_hamming
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

    # Track Hamming distances between consecutive accepted states
    hamming_distances_consecutive = Int[]

    # Track proposed Hamming distances (before acceptance)
    proposed_hamming_sum = 0

    # Store previous accepted state for Hamming distance
    prev_indices = copy(current_P_indices)

    for _ in 1:samples
        n_proposed += 1
        old_indices = copy(current_P_indices)

        current_P_indices, current_state, current_weight, accepted, prop_hamming = k_local_update!(
            current_P_indices, current_state, current_weight, psi0, k, N
        )

        proposed_hamming_sum += prop_hamming

        if accepted && prop_hamming > 0
            n_accepted += 1
            # Compute Hamming distance from previous state
            h_dist = hamming_distance(prev_indices, current_P_indices)
            push!(hamming_distances_consecutive, h_dist)
            prev_indices = copy(current_P_indices)
        elseif accepted && prop_hamming == 0
            # Trivial acceptance (no change proposed)
            n_accepted += 1
        end
    end

    # Compute statistics
    acceptance_rate = n_accepted / n_proposed

    # Average Hamming distance per accepted move
    avg_hamming_consecutive = isempty(hamming_distances_consecutive) ? 0.0 : mean(hamming_distances_consecutive)
    std_hamming_consecutive = isempty(hamming_distances_consecutive) ? 0.0 : std(hamming_distances_consecutive)

    # Average proposed Hamming distance (how many sites would change if accepted)
    avg_proposed_hamming = proposed_hamming_sum / n_proposed

    # Effective Hamming distance per step = acceptance_rate * avg_hamming_per_accepted_move
    effective_hamming_per_step = acceptance_rate * avg_hamming_consecutive

    return (
        k = k,
        acceptance_rate = acceptance_rate,
        avg_hamming_consecutive = avg_hamming_consecutive,
        std_hamming_consecutive = std_hamming_consecutive,
        avg_proposed_hamming = avg_proposed_hamming,
        effective_hamming_per_step = effective_hamming_per_step,
        n_accepted = n_accepted,
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
        acc_rates = Float64[]
        avg_hammings = Float64[]
        effective_hammings = Float64[]

        for run in 1:n_runs
            res = run_k_local_test(
                N=N, h=h, k=k, samples=samples, thermalization=thermalization
            )
            push!(acc_rates, res.acceptance_rate)
            push!(avg_hammings, res.avg_hamming_consecutive)
            push!(effective_hammings, res.effective_hamming_per_step)
        end

        push!(results, (
            k = k,
            acceptance_rate_mean = mean(acc_rates),
            acceptance_rate_std = std(acc_rates),
            avg_hamming_mean = mean(avg_hammings),
            avg_hamming_std = std(avg_hammings),
            effective_hamming_mean = mean(effective_hammings),
            effective_hamming_std = std(effective_hammings)
        ))
    end

    # Print summary table
    println()
    println("=" ^ 80)
    println("RESULTS SUMMARY")
    println("=" ^ 80)
    println()

    # Header
    @printf "%-5s | %-20s | %-20s | %-20s\n" "k" "Acceptance Rate" "Avg Hamming Dist" "Effective Hamming/Step"
    println("-" ^ 80)

    for r in results
        @printf "%-5d | %7.4f ± %-10.4f | %7.4f ± %-10.4f | %7.4f ± %-10.4f\n" r.k r.acceptance_rate_mean r.acceptance_rate_std r.avg_hamming_mean r.avg_hamming_std r.effective_hamming_mean r.effective_hamming_std
    end

    println()
    println("=" ^ 80)
    println("INTERPRETATION")
    println("=" ^ 80)
    println("""
    - Acceptance Rate: Fraction of proposed moves that are accepted
    - Avg Hamming Dist: Average number of sites changed per accepted move
    - Effective Hamming/Step: acceptance_rate × avg_hamming_dist
      This measures the average "distance traveled" in configuration space per MCMC step.
      Higher values indicate better mixing efficiency.

    Key insights:
    - k=1 updates typically have high acceptance but move only 1 site at a time
    - Larger k updates can move more sites but may have lower acceptance
    - The effective Hamming/step balances these tradeoffs
    - If effective Hamming/step increases with k, larger updates improve mixing
    - If it decreases, smaller updates are more efficient
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
# Run the test
# ============================================================

if abspath(PROGRAM_FILE) == @__FILE__
    # Run the main test
    results = k_local_update_test(N=6, h=1.0, max_k=5, samples=2000, thermalization=500, n_runs=3)

    println()
    println("Additional mixing analysis (Hamming distance from reference state):")
    println("-" ^ 60)
    for k in 1:5
        analyze_mixing(N=6, h=1.0, k=k, samples=3000, thermalization=500)
    end
end
