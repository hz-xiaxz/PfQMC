using XDiag
using LinearAlgebra
using Random
using Statistics

# Include the original logic (re-implemented here to capture trace)
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
    
    # Lookup table logic from original file
    lookup = Dict{Tuple{Int,Int},Int}()
    for a in 1:4, b in 1:4
        if a == 1; lookup[(a, b)] = b; continue; end
        if b == 1; lookup[(a, b)] = a; continue; end
        if a == b; lookup[(a, b)] = 1; continue; end
        if (a, b) in [(2, 3), (3, 2)]; lookup[(a, b)] = 4; end
        if (a, b) in [(2, 4), (4, 2)]; lookup[(a, b)] = 3; end
        if (a, b) in [(3, 4), (4, 3)]; lookup[(a, b)] = 2; end
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

function compute_autocorr(data, max_lag)
    n = length(data)
    mean_val = mean(data)
    var_val = var(data)
    
    acf = Float64[]
    for k in 0:max_lag
        s = 0.0
        for t in 1:(n-k)
            s += (data[t] - mean_val) * (data[t+k] - mean_val)
        end
        push!(acf, s / ((n-k) * var_val))
    end
    return acf
end

function run_autocorr_analysis(; N::Int=8, h::Float64=1.0, samples::Int=10000, thermalization::Int=1000)
    println("Running Autocorrelation Analysis for N=$N, h=$h, samples=$samples...")

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

    current_P_indices = rand(1:4, N)
    current_state = compute_state_from_indices(psi0, current_P_indices)
    current_overlap = dot(psi0, current_state)
    current_weight = abs2(current_overlap)

    # Thermalization
    for t in 1:thermalization
        # 2-site update only (as per user's latest code preference in conversation history)
        site1 = rand(1:N)
        site2 = rand(1:N)
        while site2 == site1; site2 = rand(1:N); end
        
        old_idx1 = current_P_indices[site1]
        old_idx2 = current_P_indices[site2]
        new_idx1 = rand(1:4)
        new_idx2 = rand(1:4)
        
        if old_idx1 != new_idx1 || old_idx2 != new_idx2
            temp_state = apply_transition(current_state, site1, old_idx1, new_idx1)
            trial_state = apply_transition(temp_state, site2, old_idx2, new_idx2)
            trial_overlap = dot(psi0, trial_state)
            trial_weight = abs2(trial_overlap)
            ratio = (current_weight < 1e-12) ? 1.0 : (trial_weight / current_weight)
            if rand() < min(1.0, ratio)
                current_P_indices[site1] = new_idx1
                current_P_indices[site2] = new_idx2
                current_state = trial_state
                current_weight = trial_weight
            end
        end
    end

    # Sampling
    trace = Float64[]
    n_acc = 0
    
    for s in 1:samples
        site1 = rand(1:N)
        site2 = rand(1:N)
        while site2 == site1; site2 = rand(1:N); end
        
        old_idx1 = current_P_indices[site1]
        old_idx2 = current_P_indices[site2]
        new_idx1 = rand(1:4)
        new_idx2 = rand(1:4)
        
        if old_idx1 != new_idx1 || old_idx2 != new_idx2
            temp_state = apply_transition(current_state, site1, old_idx1, new_idx1)
            trial_state = apply_transition(temp_state, site2, old_idx2, new_idx2)
            trial_overlap = dot(psi0, trial_state)
            trial_weight = abs2(trial_overlap)
            ratio = (current_weight < 1e-12) ? 1.0 : (trial_weight / current_weight)
            if rand() < min(1.0, ratio)
                current_P_indices[site1] = new_idx1
                current_P_indices[site2] = new_idx2
                current_state = trial_state
                current_weight = trial_weight
                n_acc += 1
            end
        else
            n_acc += 1
        end
        
        # Observable: Magic Density contribution -log(weight)/N
        # Note: This is the instantaneous value. The magic density is -log(<weight>)/N.
        # But for autocorrelation we can look at the weight itself or log(weight).
        # Let's track log(weight) as it's the local energy-like term.
        push!(trace, -log(max(1e-20, current_weight)) / N)
    end
    
    println("Acceptance Rate: $(n_acc / samples)")
    
    # Calculate Autocorrelation
    max_lag = min(2000, samples ÷ 10)
    acf = compute_autocorr(trace, max_lag)
    
    # Estimate Integrated Autocorrelation Time (tau_int)
    # Simple sum until ACF drops below 0 or small threshold
    tau_int = 0.5
    for k in 2:length(acf)
        if acf[k] < 0.05 # Threshold
            break
        end
        tau_int += acf[k]
    end
    
    println("Estimated Autocorrelation Time (tau_int): $(tau_int)")
    println("Effective Sample Size (ESS): $(samples / (2 * tau_int))")
    
    # Print first few ACF values
    println("ACF[0:10]: ", acf[1:11])
end

run_autocorr_analysis(N=8, samples=10000)
