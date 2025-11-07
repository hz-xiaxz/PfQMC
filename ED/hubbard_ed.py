#!/usr/bin/env python3
"""
Exact Diagonalization of 1D Hubbard Model using QuSpin

Hamiltonian: H = -t Σ_{i,σ} (c†_{i,σ} c_{i+1,σ} + h.c.) + U Σ_i n_{i,↑} n_{i,↓} - μ Σ_{i,σ} n_{i,σ}

Usage:
    python hubbard_ed.py --L 4 --U 4.0 --mu 0.0 --boundary OBC
"""

import numpy as np
from quspin.operators import hamiltonian
from quspin.basis import spinful_fermion_basis_1d
import argparse
import sys

def build_hubbard_hamiltonian(L, t=1.0, U=0.0, mu=0.0, boundary='OBC'):
    """
    Build 1D Hubbard model Hamiltonian

    Args:
        L: Number of lattice sites
        t: Hopping amplitude (default: 1.0)
        U: On-site interaction strength
        mu: Chemical potential
        boundary: 'OBC' (open) or 'PBC' (periodic)

    Returns:
        H: Hamiltonian operator
        basis: Fermion basis
    """

    # Define basis with conserved particle number
    # At half-filling (mu=0): N_up = L/2, N_down = L/2
    if mu == 0.0:
        if L % 2 != 0:
            raise ValueError(f"For half-filling with equal spins, L must be even. Got L={L}")
        basis = spinful_fermion_basis_1d(L=L, Nf=(L//2, L//2))
        print(f"Using particle number conservation: N_up={L//2}, N_down={L//2}")
    else:
        # Away from half-filling, use full Hilbert space
        basis = spinful_fermion_basis_1d(L)

    print(f"\n=== System Parameters ===")
    print(f"Sites: L = {L}")
    print(f"Hopping: t = {t}")
    print(f"Interaction: U = {U}")
    print(f"Chemical potential: μ = {mu}")
    print(f"Boundary: {boundary}")
    print(f"Hilbert space dimension: {basis.Ns}")

    # Hopping term: -t Σ_{i,σ} (c†_{i,σ} c_{i+1,σ} + h.c.)
    # QuSpin convention: use J = t (positive)
    # hop_left and hop_right with opposite signs
    J = t  # Hopping amplitude

    if boundary == 'OBC':
        hop_right = [[-J, i, i+1] for i in range(L-1)]  # OBC
        hop_left = [[J, i, i+1] for i in range(L-1)]    # OBC
    else:  # PBC
        hop_right = [[-J, i, (i+1)%L] for i in range(L)]  # PBC
        hop_left = [[J, i, (i+1)%L] for i in range(L)]    # PBC

    # Interaction term: U Σ_i n_{i,↑} n_{i,↓}
    # In QuSpin: use 'n|n' operator (density up times density down)
    interaction = [[U, i, i] for i in range(L)]

    # Chemical potential: -μ Σ_{i,σ} n_{i,σ}
    # In QuSpin: 'n' operator counts particles
    chemical_pot = [[-mu, i] for i in range(L)]

    # Define operator strings for QuSpin (CORRECTED)
    # Format: [operator_string, coupling_list]
    static = [
        ["+-|", hop_left],       # up hop
        ["-+|", hop_right],      # up hop (hermitian conjugate)
        ["|+-", hop_left],       # down hop
        ["|-+", hop_right],      # down hop (hermitian conjugate)
        ["n|n", interaction],    # n_{i,↑} n_{i,↓}
        ["n|", chemical_pot],    # n_{i,↑}
        ["|n", chemical_pot],    # n_{i,↓}
    ]

    # No time-dependent terms
    dynamic = []

    # Build Hamiltonian (check_herm=False to suppress hermiticity warnings)
    H = hamiltonian(static, dynamic, basis=basis, check_herm=False, dtype=np.float64)

    return H, basis


def calculate_observables(H, basis, L):
    """
    Calculate ground state observables

    Args:
        H: Hamiltonian
        basis: Fermion basis
        L: Number of sites

    Returns:
        dict: Dictionary of observables
    """

    # Diagonalize to get ground state
    print("\n=== Diagonalizing Hamiltonian ===")
    E, V = H.eigh()
    E0 = E[0]
    psi0 = V[:, 0]

    print(f"Ground state energy: E0 = {E0:.10f}")
    print(f"First excited state: E1 = {E[1]:.10f}")
    print(f"Energy gap: ΔE = {E[1] - E[0]:.10f}")

    # Particle number operators
    n_up_list = [[1.0, i] for i in range(L)]
    n_down_list = [[1.0, i] for i in range(L)]

    static_n_up = [["n|", n_up_list]]
    static_n_down = [["|n", n_down_list]]

    N_up_op = hamiltonian(static_n_up, [], basis=basis, check_herm=False, dtype=np.float64)
    N_down_op = hamiltonian(static_n_down, [], basis=basis, check_herm=False, dtype=np.float64)

    N_up = N_up_op.expt_value(psi0).real
    N_down = N_down_op.expt_value(psi0).real
    N_total = N_up + N_down

    print(f"\n=== Particle Numbers ===")
    print(f"<N_up> = {N_up:.10f}")
    print(f"<N_down> = {N_down:.10f}")
    print(f"<N_total> = {N_total:.10f}")

    # Spin structure factor S(q) = (1/L) Σ_{i,j} e^{iq(i-j)} <S^z_i S^z_j>
    # S^z_i = (n_{i,↑} - n_{i,↓}) / 2
    spin_structure_factor = {}

    # Calculate for q = π
    q = np.pi
    S_q = 0.0

    for i in range(L):
        for j in range(L):
            # Build S^z_i operator: (n_{i,↑} - n_{i,↓}) / 2
            static_szi = [["n|", [[0.5, i]]], ["|n", [[-0.5, i]]]]
            static_szj = [["n|", [[0.5, j]]], ["|n", [[-0.5, j]]]]

            Szi_op = hamiltonian(static_szi, [], basis=basis, check_herm=False, dtype=np.float64)
            Szj_op = hamiltonian(static_szj, [], basis=basis, check_herm=False, dtype=np.float64)

            # <S^z_i S^z_j> = <psi| S^z_i S^z_j |psi>
            SziSzj = Szi_op.matrix_ele(psi0, Szj_op.dot(psi0)).real

            S_q += np.exp(1j * q * (i - j)) * SziSzj

    S_q = S_q / L

    print(f"\n=== Structure Factors ===")
    print(f"Spin structure factor S(π) = {S_q.real:.10f}")

    # Charge structure factor N(q) = (1/L) Σ_{i,j} e^{iq(i-j)} <n_i n_j>
    # n_i = n_{i,↑} + n_{i,↓}
    N_q = 0.0

    for i in range(L):
        for j in range(L):
            # Build n_i operator: n_{i,↑} + n_{i,↓}
            static_ni = [["n|", [[1.0, i]]], ["|n", [[1.0, i]]]]
            static_nj = [["n|", [[1.0, j]]], ["|n", [[1.0, j]]]]

            ni_op = hamiltonian(static_ni, [], basis=basis, check_herm=False, dtype=np.float64)
            nj_op = hamiltonian(static_nj, [], basis=basis, check_herm=False, dtype=np.float64)

            # <n_i n_j>
            ninj = ni_op.matrix_ele(psi0, nj_op.dot(psi0)).real

            N_q += np.exp(1j * q * (i - j)) * ninj

    N_q = N_q / L

    print(f"Charge structure factor N(π) = {N_q.real:.10f}")

    # Double occupancy: <D> = Σ_i <n_{i,↑} n_{i,↓}> / L
    double_occ = 0.0
    for i in range(L):
        static_double = [["n|n", [[1.0, i, i]]]]
        double_op = hamiltonian(static_double, [], basis=basis, check_herm=False, dtype=np.float64)
        double_occ += double_op.expt_value(psi0).real

    double_occ /= L
    print(f"Double occupancy <D> = {double_occ:.10f}")

    # Nearest-neighbor spin correlation
    # <S_i · S_{i+1}> = <S^x_i S^x_{i+1}> + <S^y_i S^y_{i+1}> + <S^z_i S^z_{i+1}>
    # For spin-1/2: S^+ = c†_↑ c_↓, S^- = c†_↓ c_↑, S^z = (n_↑ - n_↓)/2
    # <S_i · S_{i+1}> = 1/2 <S^+_i S^-_{i+1} + S^-_i S^+_{i+1}> + <S^z_i S^z_{i+1}>

    nn_spin_corr = 0.0
    for i in range(L-1):
        # S^z_i S^z_{i+1}
        static_szi = [["n|", [[0.5, i]]], ["|n", [[-0.5, i]]]]
        static_szj = [["n|", [[0.5, i+1]]], ["|n", [[-0.5, i+1]]]]
        Szi_op = hamiltonian(static_szi, [], basis=basis, check_herm=False, dtype=np.float64)
        Szj_op = hamiltonian(static_szj, [], basis=basis, check_herm=False, dtype=np.float64)
        nn_spin_corr += Szi_op.matrix_ele(psi0, Szj_op.dot(psi0)).real

        # S^+_i S^-_{i+1} = (c†_{i,↑} c_{i,↓}) (c†_{i+1,↓} c_{i+1,↑})
        static_splus_i = [["+|-", [[1.0, i, i]]]]  # c†_↑ c_↓
        static_sminus_j = [["-|+", [[1.0, i+1, i+1]]]]  # c†_↓ c_↑
        Splus_op = hamiltonian(static_splus_i, [], basis=basis, check_herm=False, dtype=np.complex128)
        Sminus_op = hamiltonian(static_sminus_j, [], basis=basis, check_herm=False, dtype=np.complex128)
        nn_spin_corr += 0.5 * Splus_op.matrix_ele(psi0, Sminus_op.dot(psi0)).real

        # S^-_i S^+_{i+1}
        static_sminus_i = [["-|+", [[1.0, i, i]]]]
        static_splus_j = [["+|-", [[1.0, i+1, i+1]]]]
        Sminus_op_i = hamiltonian(static_sminus_i, [], basis=basis, check_herm=False, dtype=np.complex128)
        Splus_op_j = hamiltonian(static_splus_j, [], basis=basis, check_herm=False, dtype=np.complex128)
        nn_spin_corr += 0.5 * Sminus_op_i.matrix_ele(psi0, Splus_op_j.dot(psi0)).real

    nn_spin_corr /= (L - 1)
    print(f"NN spin correlation <S_i · S_{i+1}> = {nn_spin_corr:.10f}")

    results = {
        'E0': E0,
        'E1': E[1],
        'gap': E[1] - E[0],
        'N_up': N_up,
        'N_down': N_down,
        'N_total': N_total,
        'S_pi': S_q.real,
        'N_pi': N_q.real,
        'double_occupancy': double_occ,
        'nn_spin_corr': nn_spin_corr,
    }

    return results, E, psi0


def main():
    parser = argparse.ArgumentParser(description='Exact Diagonalization of 1D Hubbard Model')
    parser.add_argument('--L', type=int, default=4, help='Number of sites (default: 4)')
    parser.add_argument('--t', type=float, default=1.0, help='Hopping amplitude (default: 1.0)')
    parser.add_argument('--U', type=float, default=4.0, help='On-site interaction (default: 4.0)')
    parser.add_argument('--mu', type=float, default=0.0, help='Chemical potential (default: 0.0)')
    parser.add_argument('--boundary', type=str, default='OBC', choices=['OBC', 'PBC'],
                        help='Boundary condition (default: OBC)')
    parser.add_argument('--spectrum', action='store_true', help='Print energy spectrum')
    parser.add_argument('--save', type=str, default=None, help='Save results to file')

    args = parser.parse_args()

    print("="*60)
    print("1D Hubbard Model - Exact Diagonalization")
    print("="*60)

    # Build Hamiltonian
    H, basis = build_hubbard_hamiltonian(args.L, args.t, args.U, args.mu, args.boundary)

    # Calculate observables
    results, E, psi0 = calculate_observables(H, basis, args.L)

    # Print energy spectrum if requested
    if args.spectrum:
        print(f"\n=== Energy Spectrum (first 10 states) ===")
        for i in range(min(10, len(E))):
            print(f"E[{i}] = {E[i]:.10f}")

    # Save results to file
    if args.save:
        with open(args.save, 'w') as f:
            f.write("# 1D Hubbard Model - Exact Diagonalization Results\n")
            f.write(f"# L={args.L}, t={args.t}, U={args.U}, mu={args.mu}, boundary={args.boundary}\n")
            f.write(f"# Hilbert space dimension: {basis.Ns}\n")
            f.write("#\n")
            f.write(f"E0 = {results['E0']:.15f}\n")
            f.write(f"gap = {results['gap']:.15f}\n")
            f.write(f"N_total = {results['N_total']:.15f}\n")
            f.write(f"N_up = {results['N_up']:.15f}\n")
            f.write(f"N_down = {results['N_down']:.15f}\n")
            f.write(f"S_pi = {results['S_pi']:.15f}\n")
            f.write(f"N_pi = {results['N_pi']:.15f}\n")
            f.write(f"double_occupancy = {results['double_occupancy']:.15f}\n")
            f.write(f"nn_spin_corr = {results['nn_spin_corr']:.15f}\n")
        print(f"\nResults saved to {args.save}")

    print("\n" + "="*60)
    print("Calculation complete!")
    print("="*60)

    return results


if __name__ == "__main__":
    main()
