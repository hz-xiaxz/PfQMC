from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian
import numpy as np
from scipy import sparse  # 必须导入 scipy.sparse

def magic_quspin(N=2, samples=50):
    J_list = [[1.0, i, (i+1)%N] for i in range(N)]
    D_list = [[1.0, i, (i+1)%N] for i in range(N)]
    static = [["xx", J_list], ["yy", J_list], ["zz", D_list]]
    
    # 这里的 basis 不需要 symmetries，因为我们要算随机泡利的期望值
    basis = spin_basis_1d(L=N)
    H = hamiltonian(static, [], basis=basis, dtype=np.float64, check_symm=False, check_pcon=False)
    
    # 2. 求解基态
    print("Solving Ground State...")
    # eigsh 返回的是列向量，psi[:, 0]
    E, psi = H.eigsh(k=1, which='SA')
    psi = psi[:, 0]

    # 3. 计算 Magic (极速版)
    print(f"Sampling Magic (samples={samples})...")
    pauli_map = {0: 'I', 1: 'x', 2: 'y', 3: 'z'}
    val_sum = 0.0
    
    rand_cfgs = np.random.randint(0, 4, size=(samples, N))
    
    for k, cfg in enumerate(rand_cfgs):
        # 提取非 I 的位点和字符
        # 例如: indices=[0, 2], ops=['x', 'z']
        indices = [i for i, c in enumerate(cfg) if c != 0]
        ops = [pauli_map[cfg[i]] for i in indices]
        
        if not indices: # Identity
            exp_val = 1.0
        else:
            # 拼接算符字符串: "xz" (注意 QuSpin basis.Op 不需要 "|")
            op_str = "".join(ops)
            
            # 使用 basis.Op 直接生成稀疏矩阵
            # 签名: Op(opstr, indx, J, dtype)
            # 这里 J=1.0, indx=indices (列表)
            # 注意: indices 必须和 op_str 的字符一一对应
            P_raw= basis.Op(op_str, indices, 1.0, np.complex128)
            if isinstance(P_raw, tuple):
                # 如果是元组，说明是 (ME, row, col)，手动构造矩阵
                # 格式: (data, (row, col))
                P_mat = sparse.csr_matrix((P_raw[0], (P_raw[1], P_raw[2])), shape=(N, N))
            else:
                # 否则它已经是矩阵了
                P_mat = P_raw
            # 计算 <psi | P | psi>
            # P_mat 是稀疏矩阵，P_mat.dot(psi) 很快
            # 结果可能是复数，但 Pauli 期望值应为实数
            exp_val = np.vdot(psi, P_mat.dot(psi)).real
            
        val_sum += exp_val**4

        if (k+1) % 1000 == 0:
            print(f"Progress: {k+1}/{samples}")

    m2 = -np.log(val_sum / samples)
    return m2

print(f"SRE: {magic_quspin(N=10)}")