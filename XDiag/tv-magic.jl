using XDiag
using LinearAlgebra
function tv_magic(; N::Int=4, samples::Int=500)
    hs = Spinhalf(N)
    ops = OpSum()
    for i in 1:N
        s1 = i
        s2 = mod1(i + 1, N)
        ops += "t" * Op("Exchange", [s1, s2])
        ops += "V" * Op("SzSz", [s1, s2])
        ops += "V" * Op("Sz", [s1])
    end
    t = 1.0
    V = 1.0
    ops["t"] = -2t
    ops["V"] = V
    _, psi0 = eig0(ops, hs)
    val_sum = 0.0

    for _ in 1:samples
        pauli_ops = OpSum()
        for i in 1:N
            p = rand(1:4)
            p == 1 && continue
            if p == 2
                pauli_ops += 0.5 * Op("S+", [i])
                pauli_ops += 0.5 * Op("S-", [i])
            elseif p == 3
                pauli_ops += -0.5im * Op("S+", [i])
                pauli_ops -= -0.5im * Op("S-", [i])
            elseif p == 4
                pauli_ops += Op("Sz", [i])
            end
        end

        pauli_ops == OpSum() && continue


        # inner(A, B) 计算 <A|B>
        # 因为 P 是厄米的，期望值一定是实数，取 real 防止浮点误差
        exp_val = real(inner(pauli_ops, psi0))

        val_sum += exp_val^4
    end

    # 5. 归一化与计算熵
    # sum ≈ 4^N * (val_sum / samples)
    # Omega = 1/2^N * Sum
    # M2 = -ln(Omega) = -ln( 2^N * mean )

    mean_val = val_sum / samples
    m2 = -log(mean_val * (2^N))

    return m2
end
