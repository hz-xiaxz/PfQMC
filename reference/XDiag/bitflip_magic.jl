using XDiag
using LinearAlgebra

# --- 1. The Physics (XDiag) ---
function get_ground_state_xdiag(N::Int, h::Float64)
    hs = Spinhalf(N)
    ops = OpSum()
    J = 1.0
    for i in 1:N
        s1 = i
        s2 = mod1(i + 1, N)
        # Transverse Field Ising Model Hamiltonian
        ops += -4.0 * J * Op("SzSz", [s1, s2])
        ops += -1.0 * h * Op("S+", [s1])
        ops += -1.0 * h * Op("S-", [s1])
    end
    
    # Solve Ground State
    # eig0 returns (energy, state)
    _, psi_state = eig0(ops, hs)
    
    # EXTRACT DATA HERE:
    # vector(state) returns a standard Julia Vector
    # We convert to ComplexF64 to ensure type stability in the next step
    return ComplexF64.(vector(psi_state))
end

# --- 2. The Math (Bitwise Magic Calculation) ---
function calc_magic_density(psi::Vector{ComplexF64}, N::Int)
    dim = 2^N
    total_paulis = 4^N
    
    val_sum_ref = 0.0

    # Parallel loop over all 4^N Pauli strings
    for p_idx in 0:(total_paulis-1)
        
        # A. Decode Pauli string to bitmasks
        # x_mask: 1 where Pauli is X or Y (flips bit)
        # z_mask: 1 where Pauli is Z or Y (adds phase)
        x_mask = 0
        z_mask = 0
        temp = p_idx
        
        for k in 0:(N-1)
            code = temp & 3 # Fast modulo 4
            temp >>= 2      # Fast div 4
            
            if code == 1     # X
                x_mask |= (1 << k)
            elseif code == 2 # Y
                x_mask |= (1 << k)
                z_mask |= (1 << k)
            elseif code == 3 # Z
                z_mask |= (1 << k)
            end
        end

        # B. Compute <psi | P | psi>
        expectation = 0.0 + 0.0im
        
        for col_i in 0:(dim-1)
            row_i = col_i ⊻ x_mask
            
            phase_bit = count_ones(col_i & z_mask)
            sign = (phase_bit & 1 == 1) ? -1.0 : 1.0
            
            expectation += conj(psi[col_i + 1]) * psi[row_i + 1] * sign
        end

        # Add |<P>|^4
        val_sum_ref += abs2(abs2(expectation))
    end

    return -log(val_sum_ref / 2^N) / N
end

# --- 3. Main Execution ---

function run_benchmark()
    println("Compiling functions...")
    # Warmup with small system
    v4 = get_ground_state_xdiag(4, 1.0)
    calc_magic_density(v4, 4)

    target_N = 10
    println("\n=== Running N=$target_N ===")
    
    # Step 1: Get State
    t0 = time()
    psi = get_ground_state_xdiag(target_N, 1.0)
    println("Exact Diagonalization Time: $(round(time() - t0, digits=4))s")
    
    # Step 2: Calculate Magic
    t1 = time()
    m2 = calc_magic_density(psi, target_N)
    println("Magic Calculation Time:     $(round(time() - t1, digits=4))s")
    
    println("---------------------------")
    println("Magic Density M2 = $m2")
end

run_benchmark()