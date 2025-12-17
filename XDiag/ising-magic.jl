using XDiag
using LinearAlgebra
using Random

# 在 XDiag 中应用单个 Pauli 算符到态矢量
function apply_pauli_op(psi, site::Int, pauli_idx::Int)
    if pauli_idx == 1 # I
        return copy(psi)
    elseif pauli_idx == 2 # X = S+ + S-
        # apply 返回新向量
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
function multiply_paulis(old_idx::Int, new_idx::Int)
    if old_idx == new_idx
        return (1, 1.0 + 0im)
    end
    if old_idx == 1
        return (new_idx, 1.0 + 0im)
    end
    if new_idx == 1
        return (old_idx, 1.0 + 0im)
    end # Pauli 是厄米的，所以 P*I = P

    # 乘法表: T = new * old
    # X(2), Y(3), Z(4)
    # XY=iZ, YZ=iX, ZX=iY
    # YX=-iZ, ZY=-iX, XZ=-iY
    if (new_idx, old_idx) == (2, 3)
        return (4, -1.0im)
    end  # Y* (-iZ) = X 
    if (new_idx, old_idx) == (3, 2)
        return (4, 1.0im)
    end # X * (iZ) = Y 

    if (new_idx, old_idx) == (3, 4)
        return (2, -1.0im)
    end  # Z *(-iX) =Y
    if (new_idx, old_idx) == (4, 3)
        return (2, 1.0im)
    end #  Y *(iX) = Z

    if (new_idx, old_idx) == (4, 2)
        return (3, -1.0im)
    end  # X *(-iZ) = Y
    if (new_idx, old_idx) == (2, 4)
        return (3, 1.0im)
    end # Z * (iX) = Y 

    error("Unreachable Pauli multiplication")
end

function apply_transition(current_state, site::Int, old_idx::Int, new_idx::Int)
    if old_idx == new_idx
        return current_state
    end

    # 计算 P_new * P_old 对应的 Pauli 算符 (忽略全局相位 ±1, ±i，因为我们只关心概率模长)
    # Multiplication table (indices):
    # I(1) * X(2) = X(2), X(2) * I(1) = X(2), etc.
    # X(2) * X(2) = I(1)
    # X(2) * Y(3) = iZ(4) -> index 4
    # ...

    # 简单查找表 (old, new) -> transition_op_index
    # 这里的逻辑是: 我们需要应用哪个算符 T，使得 T * P_old ~ P_new
    # 即 T ~ P_new * P_old

    lookup = Dict{Tuple{Int,Int},Int}()

    # 填充查找表
    for a in 1:4, b in 1:4
        if a == 1
            lookup[(a, b)] = b
            continue
        end # 1 * b = b
        if b == 1
            lookup[(a, b)] = a
            continue
        end # a * 1 = a
        if a == b
            lookup[(a, b)] = 1
            continue
        end # a * a = I

        # Pauli 乘法表 (1:I, 2:X, 3:Y, 4:Z)
        # (2,3)->4, (2,4)->3, (3,2)->4, (3,4)->2, (4,2)->3, (4,3)->2
        if (a, b) in [(2, 3), (3, 2)]
            lookup[(a, b)] = 4
        end
        if (a, b) in [(2, 4), (4, 2)]
            lookup[(a, b)] = 3
        end
        if (a, b) in [(3, 4), (4, 3)]
            lookup[(a, b)] = 2
        end
    end

    trans_idx = lookup[(old_idx, new_idx)]

    # 应用转移算符
    # 注意：这里忽略了相位(i, -i, -1)，对于 Metropolis 接受率 p ~ |<psi|P|psi>|^2 没有影响
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


function ising_magic_MCMC(; N::Int=4, h::Float64=1.0, samples::Int=1000, thermalization::Int=200)
    println("Running Pauli Sampling with XDiag for N=$N, h=$h...")

    # ==========================================
    # 0. 准备工作: 定义 Hamiltonian 并获取基态 ρ = |ψ⟩⟨ψ|
    # ==========================================
    block = Spinhalf(N)

    # 定义 TFIM Hamiltonian: H = -J ∑ Z_i Z_{i+1} - h ∑ X_i
    ops = OpSum()
    J = 1.0
    # h is passed as argument

    for i in 1:N
        j = mod1(i + 1, N)
        ops += -4.0 * J * Op("SzSz", [i, j])
        ops += -1.0 * h * Op("S+", [i])
        ops += -1.0 * h * Op("S-", [i])
    end

    # 计算基态 |ψ0⟩
    _, psi0 = eig0(ops, block)

    # ==========================================
    # Algorithm 1: Monte Carlo sampling
    # ==========================================

    # --- Step 1: Initialize Pauli string P ---
    # 随机初始化 P
    current_P_indices = rand(1:4, N)

    # 构造初始状态 |ϕ⟩ = P |ψ0⟩
    current_state = compute_state_from_indices(psi0, current_P_indices)

    # --- Step 2: Compute Tr(ρP) and Π_P ---
    current_overlap = dot(psi0, current_state)
    current_weight = abs2(current_overlap)

    # 预热 (Thermalization)
    for t in 1:thermalization
        # Decide move type: 50% 1-site, 50% 2-site
        # if rand() < 0.5
        #     # --- 1-site update ---
        #     site = rand(1:N)
        #     old_idx = current_P_indices[site]
        #     new_idx = rand(1:4)

        #     if old_idx != new_idx
        #         # Incremental update
        #         trial_state = apply_transition(current_state, site, old_idx, new_idx)
        #         trial_overlap = dot(psi0, trial_state)
        #         trial_weight = abs2(trial_overlap)

        #         ratio = (current_weight < 1e-12) ? 1.0 : (trial_weight / current_weight)

        #         if rand() < min(1.0, ratio)
        #             current_P_indices[site] = new_idx
        #             current_state = trial_state
        #             current_weight = trial_weight
        #         end
        #     end
        # else
            # --- 2-site update ---
            site1 = rand(1:N)
            site2 = rand(1:N)
            while site2 == site1
                site2 = rand(1:N)
            end
            
            old_idx1 = current_P_indices[site1]
            old_idx2 = current_P_indices[site2]
            
            new_idx1 = rand(1:4)
            new_idx2 = rand(1:4)
            
            if old_idx1 != new_idx1 || old_idx2 != new_idx2
                # Incremental update (apply twice)
                # Note: Order doesn't matter for distinct sites
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
        # end
    end

    # 采样 (Sampling)
    weight_sum = 0.0 # 用于计算 entropy 的 estimator
    
    # n_prop_1 = 0
    # n_acc_1 = 0
    n_prop_2 = 0
    n_acc_2 = 0

    # --- Step 3: Loop ---
    for s in 1:samples
        # Decide move type: 50% 1-site, 50% 2-site
        # if rand() < 0.5
        #     # --- 1-site update ---
        #     n_prop_1 += 1
        #     site = rand(1:N)
        #     old_idx = current_P_indices[site]
        #     new_idx = rand(1:4)

        #     if old_idx != new_idx
        #         # Incremental update
        #         trial_state = apply_transition(current_state, site, old_idx, new_idx)
        #         trial_overlap = dot(psi0, trial_state)
        #         trial_weight = abs2(trial_overlap)

        #         ratio = (current_weight < 1e-12) ? 1.0 : (trial_weight / current_weight)

        #         if rand() < min(1.0, ratio)
        #             current_P_indices[site] = new_idx
        #             current_state = trial_state
        #             current_weight = trial_weight
        #             n_acc_1 += 1
        #         end
        #     else
        #         # Proposing same state is accepted (no change)
        #         n_acc_1 += 1
        #     end
        # else
            # --- 2-site update ---
            n_prop_2 += 1
            site1 = rand(1:N)
            site2 = rand(1:N)
            while site2 == site1
                site2 = rand(1:N)
            end
            
            old_idx1 = current_P_indices[site1]
            old_idx2 = current_P_indices[site2]
            
            new_idx1 = rand(1:4)
            new_idx2 = rand(1:4)
            
            if old_idx1 != new_idx1 || old_idx2 != new_idx2
                # Incremental update (apply twice)
                # Note: Order doesn't matter for distinct sites
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
                    n_acc_2 += 1
                end
            else
                n_acc_2 += 1
            end
        # end

        # --- Step 7: Measure ---
        weight_sum += current_weight
    end

    # 计算结果
    avg_weight = weight_sum / samples
    
    m2 = -log(avg_weight) / N
    
    # println("Acceptance Rate 1-site: $(n_acc_1 / max(1, n_prop_1))")
    println("Acceptance Rate 2-site: $(n_acc_2 / max(1, n_prop_2))")

    return m2
end
