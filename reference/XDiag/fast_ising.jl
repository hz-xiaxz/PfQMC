using XDiag
using LinearAlgebra
using Random
using Printf

# ==========================================
# Optimized Pauli Application (Zero Allocation)
# ==========================================
function apply_pauli_fast!(out::Vector{ComplexF64}, inp::Vector{ComplexF64}, site::Int, pauli_idx::Int, N::Int)
    if pauli_idx == 1 # I
        copyto!(out, inp)
        return
    end

    stride = 1 << (site - 1)
    dim = length(inp)
    step = 2 * stride
    
    @inbounds for base in 1:step:dim
        for offset in 0:(stride-1)
            idx0 = base + offset       # Spin UP (0) at site
            idx1 = idx0 + stride       # Spin DOWN (1) at site
            
            c0 = inp[idx0]
            c1 = inp[idx1]
            
            if pauli_idx == 2 # X
                out[idx0] = c1
                out[idx1] = c0
            elseif pauli_idx == 3 # Y
                out[idx0] = -1.0im * c1
                out[idx1] =  1.0im * c0
            elseif pauli_idx == 4 # Z
                out[idx0] =  c0
                out[idx1] = -c1
            end
        end
    end
end

# Transition lookup table (const)
const TRANSITION_LOOKUP = let
    tbl = Dict{Tuple{Int,Int},Int}()
    for a in 1:4, b in 1:4
        if a == 1; tbl[(a, b)] = b; continue; end
        if b == 1; tbl[(a, b)] = a; continue; end
        if a == b; tbl[(a, b)] = 1; continue; end
        if (a, b) in [(2, 3), (3, 2)]; tbl[(a, b)] = 4; end
        if (a, b) in [(2, 4), (4, 2)]; tbl[(a, b)] = 3; end
        if (a, b) in [(3, 4), (4, 3)]; tbl[(a, b)] = 2; end
    end
    tbl
end

function get_transition_op(old_idx::Int, new_idx::Int)
    return TRANSITION_LOOKUP[(old_idx, new_idx)]
end

function fast_ising_magic(; N::Int=8, h::Float64=1.0, samples::Int=100000, thermalization::Int=1000)
    println("Running FAST Ising Magic MCMC for N=$N, h=$h, samples=$samples...")

    # 1. Setup Hamiltonian & Ground State
    block = Spinhalf(N)
    ops = OpSum()
    J = 1.0
    for i in 1:N
        j = mod1(i + 1, N)
        ops += -4.0 * J * Op("SzSz", [i, j])
        ops += -1.0 * h * Op("S+", [i])
        ops += -1.0 * h * Op("S-", [i])
    end
    
    # Use full diagonalization to get Julia Vector
    # This works for N <= 12 comfortably
    H_mat = matrix(ops, block)
    vals, vecs = eigen(H_mat)
    
    # Ground state is the first eigenvector
    # Convert to ComplexF64 to support Y operations
    psi0 = ComplexF64.(vecs[:, 1])
    
    # 2. Initialize MCMC
    current_P_indices = rand(1:4, N)
    
    # Pre-allocate state vectors
    current_state = copy(psi0)
    tmp = similar(current_state)
    
    # Apply initial Pauli string
    for (site, p_idx) in enumerate(current_P_indices)
        if p_idx != 1
            apply_pauli_fast!(tmp, current_state, site, p_idx, N)
            copyto!(current_state, tmp)
        end
    end
    
    trial_state = similar(current_state)
    temp_state = similar(current_state)
    
    current_overlap = dot(psi0, current_state)
    current_weight = abs2(current_overlap)
    
    # 3. MCMC Loop
    n_acc = 0
    weight_sum = 0.0
    
    total_steps = thermalization + samples
    
    for step in 1:total_steps
        # 2-site update
        site1 = rand(1:N)
        site2 = rand(1:N)
        while site2 == site1; site2 = rand(1:N); end
        
        old_idx1 = current_P_indices[site1]
        old_idx2 = current_P_indices[site2]
        new_idx1 = rand(1:4)
        new_idx2 = rand(1:4)
        
        if old_idx1 != new_idx1 || old_idx2 != new_idx2
            # Apply transition 1: current -> temp
            trans1 = get_transition_op(old_idx1, new_idx1)
            apply_pauli_fast!(temp_state, current_state, site1, trans1, N)
            
            # Apply transition 2: temp -> trial
            trans2 = get_transition_op(old_idx2, new_idx2)
            apply_pauli_fast!(trial_state, temp_state, site2, trans2, N)
            
            trial_overlap = dot(psi0, trial_state)
            trial_weight = abs2(trial_overlap)
            
            ratio = (current_weight < 1e-12) ? 1.0 : (trial_weight / current_weight)
            
            if rand() < min(1.0, ratio)
                # Accept
                current_P_indices[site1] = new_idx1
                current_P_indices[site2] = new_idx2
                
                copyto!(current_state, trial_state)
                
                current_weight = trial_weight
                if step > thermalization
                    n_acc += 1
                end
            else
                # Reject
            end
        else
            if step > thermalization
                n_acc += 1
            end
        end
        
        if step > thermalization
            weight_sum += current_weight
        end
    end
    
    avg_weight = weight_sum / samples
    magic_density = -log(avg_weight) / N
    
    println("Fast MCMC Result: $magic_density")
    println("Acceptance Rate: $(n_acc / samples)")
    
    return magic_density
end

# Run Benchmark
using BenchmarkTools
println("Benchmarking Fast Implementation...")
# Run a small case first to compile
fast_ising_magic(N=4, samples=1000, thermalization=100)

# Run the target case
@time fast_ising_magic(N=8, samples=1000000, thermalization=10000)
