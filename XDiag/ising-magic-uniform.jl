using XDiag
using LinearAlgebra

function ising_magic(; N::Int=8, h::Float64=1.0)
    # 1. 构造哈密顿量 (不变)
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
    _, psi0 = eig0(ops, hs)

    total_paulis = 4^N
    val_sum = 0.0

    # 2. 遍历 Pauli String
    for idx in 0:(total_paulis-1)
        temp_idx = idx

        # 复制一份基态用于演化 (或者使用不可变更新)
        # apply(op, psi) 返回新向量，不会修改 psi0
        current_psi = psi0
        coeff = 1.0 + 0.0im # 追踪系数 (比如 Pauli Y 的 -i)

        # 逐个格点应用算符
        for i in 1:N
            code = temp_idx % 4
            temp_idx = div(temp_idx, 4)

            if code == 0
                # Identity: 没有任何操作
                continue

            elseif code == 1
                # Pauli X = S+ + S-
                # 逻辑：psi = apply(S+, psi) + apply(S-, psi)
                # 这种写法虽然会有内存分配，但在 XDiag 中是标准写法
                v_plus = apply(Op("S+", [i]), current_psi)
                v_minus = apply(Op("S-", [i]), current_psi)
                current_psi = v_plus + v_minus

            elseif code == 2
                # Pauli Y = -i (S+ - S-)
                v_plus = apply(Op("S+", [i]), current_psi)
                v_minus = apply(Op("S-", [i]), current_psi)
                # 更新波函数，并把 -i 乘进去
                current_psi = v_plus - v_minus
                coeff *= -1.0im

            elseif code == 3
                # Pauli Z = 2 * Sz
                current_psi = apply(Op("Sz", [i]), current_psi)
                coeff *= 2.0
            end
        end

        # 计算重叠
        overlap = dot(psi0, current_psi) * coeff
        val_sum += abs2(abs2(overlap))
    end

    return -log(val_sum / 2^N) / N
end
using Printf

function ising_magic_finite_T(; N::Int=8, beta::Int=8)
    # ---------------------------------------------------------
    # 1. 构造哈密顿量 (PBC)
    # ---------------------------------------------------------
    hs = Spinhalf(N)
    ops = OpSum()
    J = 1.0; h = 1.0
    for i in 1:N
        s1 = i; s2 = mod1(i+1, N)
        ops += -4.0 * J * Op("SzSz", [s1, s2])
        ops += -1.0 * h * Op("S+", [s1])
        ops += -1.0 * h * Op("S-", [s1])
    end

    # 注意：我们需要所有本征值来构建热密度矩阵，所以必须转为稠密矩阵
    H_mat = matrix(ops, hs)
    
    # 全对角化
    eig_vals, eig_vecs = eigen(H_mat)

    # ---------------------------------------------------------
    # 2. 构建热密度矩阵 ρ (β = N)
    # ---------------------------------------------------------
    
    # 玻尔兹曼权重 P_n = exp(-β E_n)
    # 为了数值稳定性，先减去最小能量: exp(-β (E_n - E_0))
    E0 = minimum(eig_vals)
    weights = exp.(-beta .* (eig_vals .- E0))
    Z = sum(weights) # 配分函数
    weights ./= Z    # 归一化
    
    # 构建 ρ = Σ p_n |n><n|
    # 这里不需要显式构建巨大的 ρ 矩阵，我们只需要知道它在能基下的形式
    # 但为了计算 Tr(ρ P) 方便，在 N=8 时直接构建矩阵是可以的 (256x256)
    rho = eig_vecs * Diagonal(weights) * eig_vecs'

    # ---------------------------------------------------------
    # 3. 计算 Renyi Entropy S2(ρ)
    # ---------------------------------------------------------
    # S2 = -log(Tr(ρ^2))
    # Tr(ρ^2) = Σ p_n^2
    tr_rho2 = sum(abs2, weights)
    S2 = -log(tr_rho2)

    # ---------------------------------------------------------
    # 4. 计算 Pauli Magic M2(ρ)
    # ---------------------------------------------------------
    total_paulis = 4^N
    sum_tr_rho_P_4 = 0.0

    # 遍历 Pauli String
    # 优化：为了速度，我们可以直接把 Pauli 矩阵构建出来做 Trace
    # N=8, Matrix size 256. Trace很快。
    
    # 预定义 Pauli 矩阵
    Id = [1.0 0.0; 0.0 1.0]
    Sx = [0.0 1.0; 1.0 0.0]
    Sy = [0.0 -1.0im; 1.0im 0.0]
    Sz = [1.0 0.0; 0.0 -1.0]
    bases = [Id, Sx, Sy, Sz]

    # 使用 index 递归或者 base-4 循环
    # 这里用简单的循环构造 kron (虽然稍慢但逻辑最清晰)
    for idx in 0:(total_paulis - 1)
        temp_idx = idx
        
        # 我们可以不用 kron 整个大矩阵，利用 trace 技巧:
        # Tr(ρ P) = Σ <n|ρ|m> <m|P|n> ... 比较麻烦
        # 直接 kron 出 P 矩阵做 trace 吧，256x256 很快
        
        # 快速构建 P 的 Kronecker 积
        # 注意：Julia 的 kron 是反序的? 需要确认 XDiag 基底顺序
        # XDiag 通常是 site 1 在低位还是高位？
        # 一般来说，我们之前计算 pure state 没遇到顺序问题是因为是对称的 Ising
        # 这里为了保险，我们可以不用 matrix trace，而是用 sum_n p_n <n|P|n>
        
        expectation_val = 0.0 + 0.0im
        
        # 计算 Tr(ρ P) = Σ_k p_k <k|P|k>
        # 这比矩阵乘法 P * rho 快，因为 P 是稀疏的 (Product state)
        # 我们只需要对每个本征态 |k> 算 <k|P|k> 然后加权
        
        # 这里为了避免 N*2^N*4^N 的复杂度，我们还是退回到: 
        # N=8 很小，直接 matrix trace 比较容易写
        
        mats = Vector{Matrix{ComplexF64}}(undef, N)
        y_count = 0
        
        for i in 1:N
            code = temp_idx % 4
            temp_idx = div(temp_idx, 4)
            mats[i] = bases[code + 1]
            if code == 2; y_count += 1; end
        end
        
        # 优化：Ising 模型哈密顿量是实的，本征向量是实的
        # 如果 Pauli 串含有奇数个 Y，则 <n|P|n> 为纯虚数
        # 且由于对称性 Tr(rho P) 应该为 0
        if isodd(y_count)
            continue
        end

        # 构建 Pauli 矩阵 P
        # 注意 Kronecker 顺序：Matrix(ops) 默认是 1 ⊗ 2 ... 还是 N ⊗ ... 1?
        # XDiag/Julia 通常 site 1 是最内层 (stride 1) 或者最外层
        # XDiag standard: site 1 is least significant bit (stride 1).
        # 这意味着 kron 应该是 kron(PN, ..., P1) 因为 Julia 是 column-major
        # 或者直接用 reduce(kron, reverse(mats))
        P_mat = reduce(kron, reverse(mats)) 
        
        # Tr(ρ P) = dot(P_mat', rho) ? No, trace(P * rho)
        # trace(A*B) = sum(A .* B')
        val = tr(rho * P_mat)
        
        sum_tr_rho_P_4 += abs2(abs2(val))
    end

    M2 = -log(sum_tr_rho_P_4 / 2^N)
    
    # ---------------------------------------------------------
    # 5. 最终结果 M~2 = M2 - S2
    # ---------------------------------------------------------
    magic_density = (M2 - S2) / N
    
    return magic_density, M2/N, S2/N
end

# 运行
# res, m2, s2 = ising_magic_finite_T(N=8)
# @printf "N=8, β=8 (Finite T):\n"
# @printf "Full SRE Density (M2-S2)/N : %.5f\n" res
# @printf "Raw Magic Density M2/N      : %.5f\n" m2
# @printf "Renyi Entropy Density S2/N  : %.5f\n" s2