#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from quspin.operators import hamiltonian
from quspin.basis import spinless_fermion_basis_1d

def run_fidelity_L10():
    # === 参数 ===
    L = 26  # 使用刚才验证过的 L=10
    N_f = L // 2
    t = 1.0
    
    # 细致扫描 V，特别是在 2.0 - 3.0 区域
    V_list = np.linspace(0.1, 5.0, 20) 
    dV = V_list[1] - V_list[0]
    
    print(f"=== Fidelity Susceptibility Scan (L={L}, PBC) ===")
    print("Detecting how fast the Ground State Wavefunction changes...")
    
    # 既然刚才验证了 L=10 时 k=0 永远是基态
    # 我们就锁定 k=0 扇区计算波函数
    basis = spinless_fermion_basis_1d(L=L, Nf=N_f, kblock=0)
    
    hop_pm = [[-t, i, (i+1)%L] for i in range(L)]
    hop_mp = [[+t, i, (i+1)%L] for i in range(L)]
    hopping_static = [["+-", hop_pm], ["-+", hop_mp]]

    fidelity_list = []
    psi_prev = None
    
    for V in V_list:
        interaction_list = [[V, i, (i+1)%L] for i in range(L)]
        interaction_static = [["nn", interaction_list]]
        static = hopping_static + interaction_static
        dynamic = []
        
        # 必须用 complex128
        H = hamiltonian(static, dynamic, basis=basis, dtype=np.complex128, check_herm=False)
        
        # 计算基态
        E, vectors = H.eigsh(k=1, which='SA')
        psi_current = vectors[:, 0]
        
        if psi_prev is not None:
            # 计算重叠 <psi(V) | psi(V+dV)>
            overlap = np.abs(np.vdot(psi_prev, psi_current))
            
            # Fidelity Susceptibility
            chi = 2 * (1.0 - overlap) / (dV**2)
            fidelity_list.append(chi)
        else:
            fidelity_list.append(np.nan)
            
        psi_prev = psi_current

    # === 绘图 ===
    V_plot = V_list[1:]
    chi_plot = fidelity_list[1:]

    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 画出 Fidelity
    ax.plot(V_plot, chi_plot, 'o-', color='crimson', markersize=4, label='Fidelity Susceptibility $\chi_F$')
    
    # 寻找峰值位置
    peak_idx = np.argmax(chi_plot)
    peak_V = V_plot[peak_idx]
    
    ax.axvline(x=peak_V, color='black', linestyle='--', label=f'Peak @ V={peak_V:.2f}')
    ax.axvline(x=2.0, color='green', linestyle='--', alpha=0.5, label='Thermodynamic $V_c=2t$')
    
    ax.set_xlabel('Interaction $V/t$', fontsize=14)
    ax.set_ylabel('Wavefunction Change Rate $\chi_F$', fontsize=14)
    ax.set_title(f'Phase Transition Detection (L={L}, PBC)\nLooking at Wavefunction instead of Energy', fontsize=16)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run_fidelity_L10()