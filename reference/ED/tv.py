#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from quspin.operators import hamiltonian
from quspin.basis import spinless_fermion_basis_1d

def calculate_structure_factor(basis, psi, L):
    # S(pi) = (1/L) * < (Sum (-1)^i n_i)^2 >
    # 展开后其实是关联函数 <n_i n_j> 的傅里叶变换
    staggered_density_list = [[(-1.0)**i, i] for i in range(L)]
    static_op = [["n", staggered_density_list]]
    O_op = hamiltonian(static_op, [], basis=basis, check_herm=True, dtype=np.float64)
    return (O_op * O_op).expt_value(psi).real / L

def run_corrected_tV_model():
    # === 参数 ===
    L = 12
    N_f = L // 2
    t = 1.0
    # 细化扫描点，特别是在相变点 V=2 附近
    V_list = np.linspace(0, 5.0, 51)
    
    basis = spinless_fermion_basis_1d(L=L, Nf=N_f)
    print(f"System: L={L}, Nf={N_f}, APBC Boundary")
    
    # === 关键修正：APBC (反周期边界条件) ===
    # 使得 hopping 在跨越边界时变号: c^dag_L c_1 -> - c^dag_L c_1
    # 这样可以消除 L=4n 系统的费米面简并
    hop_list = []
    for i in range(L - 1):
        hop_list.append([-t, i, i+1])     # 内部跳跃: -t
        hop_list.append([+t, i, i+1])     # h.c. (+t because -t* = -t, wait... logic below)
    
    # 边界跳跃 (L-1 -> 0): 符号反转
    # 正常是 -t，APBC 变成 +t
    hop_list.append([+t, L-1, 0]) 
    hop_list.append([-t, L-1, 0]) # h.c.

    # QuSpin 的 hopping 是 hermitian conjugate 自动处理吗？
    # 最好显式写出 +- 和 -+。
    # Hamiltonian: -t \sum (c^dag_i c_{i+1} + h.c.)
    # APBC: bond (L-1, 0) has hopping +t
    
    hop_pm = [] # c^dag_i c_j
    hop_mp = [] # c^dag_j c_i
    
    for i in range(L-1):
        hop_pm.append([-t, i, i+1]) # -t c^dag_i c_{i+1}
        hop_mp.append([-t, i, i+1]) # -t c^dag_{i+1} c_i (Hermitian conj of above? No, real t)
        # 解释：(-t c^dag_i c_{i+1}) 的厄米共轭是 (-t c^dag_{i+1} c_i)
        
    # 边界项 (APBC): 系数变号，由 -t 变为 +t
    hop_pm.append([+t, L-1, 0]) 
    hop_mp.append([+t, L-1, 0]) 

    hopping_static = [
        ["+-", hop_pm], 
        ["-+", hop_mp]
    ]

    S_pi_list = []
    gap01_list = [] # E1 - E0 (Symmetry breaking gap)
    gap02_list = [] # E2 - E0 (Excitation gap)

    for V in V_list:
        # 相互作用项 (不变)
        interaction_list = [[V, i, (i+1)%L] for i in range(L)]
        interaction_static = [["nn", interaction_list]]
        
        static = hopping_static + interaction_static
        dynamic = []
        
        H = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_herm=False)
        
        # 计算前3个能级
        E, vectors = H.eigh()
        
        E0, E1, E2 = E[0], E[1], E[2]
        psi0 = vectors[:, 0]
        
        S_pi_list.append(calculate_structure_factor(basis, psi0, L))
        gap01_list.append(E1 - E0)
        gap02_list.append(E2 - E0)

    # === 绘图 ===
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # 1. 结构因子
    color = 'tab:red'
    ax1.set_xlabel('Interaction Strength V / t', fontsize=14)
    ax1.set_ylabel('Structure Factor $S(\pi)$', color=color, fontsize=14)
    ax1.plot(V_list, S_pi_list, 'o-', color=color, markersize=4, label='CDW Order $S(\pi)$')
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_ylim(0, max(S_pi_list)*1.1)

    # 2. 能隙 (双Y轴)
    ax2 = ax1.twinx()
    color_g1 = 'tab:blue'
    color_g2 = 'tab:green'
    
    ax2.set_ylabel('Energy Gaps', color='k', fontsize=14)
    
    # 画出 E2-E0 (Excitation Gap)
    line2, = ax2.plot(V_list, gap02_list, 's--', color=color_g2, alpha=0.8, label='Excitation Gap $\Delta_{02}$')
    
    # 画出 E1-E0 (Splitting)
    line3, = ax2.plot(V_list, gap01_list, '^:', color=color_g1, alpha=0.6, label='Splitting $\Delta_{01}$')
    
    # 理论相变线
    ax1.axvline(x=2.0, color='gray', linestyle='--', alpha=0.5)
    ax1.text(2.05, 0.5, 'Theoretical $V_c=2t$', rotation=90, color='gray')

    # Legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2 = [line2, line3]
    labels2 = [l.get_label() for l in lines2]
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

    plt.title(f'Spinless t-V Model (L={L}, APBC)\nLL to CDW Transition', fontsize=16)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run_corrected_tV_model()