#!/usr/bin/env python3
"""
exact diagonalization of 1d hubbard model using quspin with translation symmetry

this version uses momentum (k) symmetry to reduce the hilbert space dimension.
the hamiltonian block-diagonalizes into l separate momentum sectors.

hamiltonian: h = -t σ_{i,σ} (c†_{i,σ} c_{i+1,σ} + h.c.) + u σ_i n_{i,↑} n_{i,↓} - μ σ_{i,σ} n_{i,σ}

usage:
    python hubbard_ed_ksymmetry.py --l 4 --u 4.0 --mu 0.0 --boundary pbc
"""

import numpy as np
from quspin.operators import hamiltonian
from quspin.basis import spinful_fermion_basis_1d
import argparse


def build_hubbard_hamiltonian_k(l, k, t=1.0, u=0.0, mu=0.0, boundary='pbc'):
    """
    build 1d hubbard model hamiltonian in momentum sector k

    args:
        l: number of lattice sites
        k: momentum quantum number (0, 1, ..., l-1)
        t: hopping amplitude (default: 1.0)
        u: on-site interaction strength
        mu: chemical potential
        boundary: 'pbc' (periodic) only - translation symmetry requires pbc

    returns:
        h: hamiltonian operator in k-sector
        basis: fermion basis with momentum k
    """

    if boundary != 'pbc':
        raise valueerror("translation symmetry requires pbc boundary conditions")

    # define basis with conserved particle number and momentum
    # at half-filling (mu=0): n_up = l/2, n_down = l/2
    if mu == 0.0:
        if l % 2 != 0:
            raise valueerror(f"for half-filling with equal spins, l must be even. got l={l}")
        basis = spinful_fermion_basis_1d(l=l, nf=(l//2, l//2), kblock=k)
    else:
        # away from half-filling, use full hilbert space with k symmetry
        basis = spinful_fermion_basis_1d(l=l, kblock=k)

    # hopping term: -t σ_{i,σ} (c†_{i,σ} c_{i+1,σ} + h.c.)
    j = t  # hopping amplitude
    hop_right = [[-j, i, (i+1)%l] for i in range(l)]  # pbc
    hop_left = [[j, i, (i+1)%l] for i in range(l)]    # pbc

    # interaction term: u σ_i n_{i,↑} n_{i,↓}
    interaction = [[u, i, i] for i in range(l)]

    # chemical potential: -μ σ_{i,σ} n_{i,σ}
    chemical_pot = [[-mu, i] for i in range(l)]

    # define operator strings for quspin
    static = [
        ["+-|", hop_left],       # up hop
        ["-+|", hop_right],      # up hop (hermitian conjugate)
        ["|+-", hop_left],       # down hop
        ["|-+", hop_right],      # down hop (hermitian conjugate)
        ["n|n", interaction],    # n_{i,↑} n_{i,↓}
        ["n|", chemical_pot],    # n_{i,↑}
        ["|n", chemical_pot],    # n_{i,↓}
    ]

    dynamic = []

    # build hamiltonian
    h = hamiltonian(static, dynamic, basis=basis, check_herm=false, dtype=np.float64)

    return h, basis


def calculate_observables_k(h, basis, l, k):
    """
    calculate ground state observables in momentum sector k

    args:
        h: hamiltonian in k-sector
        basis: fermion basis
        l: number of sites
        k: momentum quantum number

    returns:
        dict: dictionary of observables
    """

    # diagonalize to get ground state in this k-sector
    e, v = h.eigh()
    e0 = e[0]
    psi0 = v[:, 0]

    # particle number operators
    n_up_list = [[1.0, i] for i in range(l)]
    n_down_list = [[1.0, i] for i in range(l)]

    static_n_up = [["n|", n_up_list]]
    static_n_down = [["|n", n_down_list]]

    n_up_op = hamiltonian(static_n_up, [], basis=basis, check_herm=false, dtype=np.float64)
    n_down_op = hamiltonian(static_n_down, [], basis=basis, check_herm=false, dtype=np.float64)

    n_up = n_up_op.expt_value(psi0).real
    n_down = n_down_op.expt_value(psi0).real
    n_total = n_up + n_down

    # double occupancy: <d> = σ_i <n_{i,↑} n_{i,↓}> / l
    double_occ = 0.0
    for i in range(l):
        static_double = [["n|n", [[1.0, i, i]]]]
        double_op = hamiltonian(static_double, [], basis=basis, check_herm=false, dtype=np.float64)
        double_occ += double_op.expt_value(psi0).real

    double_occ /= l

    results = {
        'k': k,
        'e0': e0,
        'e1': e[1] if len(e) > 1 else np.nan,
        'gap': e[1] - e[0] if len(e) > 1 else np.nan,
        'n_up': n_up,
        'n_down': n_down,
        'n_total': n_total,
        'double_occupancy': double_occ,
        'hilbert_dim': basis.ns,
    }

    return results, e, psi0


def hubbard_ed_k0(l, t=1.0, u=0.0, mu=0.0, verbose=true):
    """
    quick calculation for k=0 momentum sector only (most commonly used).

    args:
        l: number of lattice sites
        t: hopping amplitude (default: 1.0)
        u: on-site interaction strength
        mu: chemical potential (default: 0.0)
        verbose: print results (default: true)

    returns:
        dict: dictionary with observables:
            - e0: ground state energy
            - gap: energy gap to first excited state
            - n_up, n_down, n_total: particle numbers
            - double_occupancy: <σ_i n_{i,↑} n_{i,↓}> / l
            - hilbert_dim: hilbert space dimension

    example:
        >>> results = hubbard_ed_k0(l=4, u=4.0, mu=0.0)
        >>> print(f"ground state energy: {results['e0']:.6f}")
    """

    k = 0  # k=0 momentum sector

    if verbose:
        print("="*60)
        print(f"1d hubbard model ed: k=0 sector")
        print("="*60)
        print(f"l={l}, t={t}, u={u}, μ={mu}")

    # build hamiltonian for k=0
    h, basis = build_hubbard_hamiltonian_k(l, k, t, u, mu, boundary='pbc')

    if verbose:
        print(f"hilbert space dimension: {basis.ns} (reduced from full space)")

    # calculate observables
    results, e, psi0 = calculate_observables_k(h, basis, l, k)

    if verbose:
        print(f"\n=== results ===")
        print(f"e0 = {results['e0']:.10f}")
        if not np.isnan(results['gap']):
            print(f"gap = {results['gap']:.10f}")
        print(f"<n_up> = {results['n_up']:.6f}")
        print(f"<n_down> = {results['n_down']:.6f}")
        print(f"<n_total> = {results['n_total']:.6f}")
        print(f"double occupancy = {results['double_occupancy']:.6f}")
        print("="*60)

    return results


def main():
    parser = argparse.argumentparser(
        description='exact diagonalization of 1d hubbard model with translation symmetry'
    )
    parser.add_argument('--l', type=int, default=4, help='number of sites (default: 4)')
    parser.add_argument('--t', type=float, default=1.0, help='hopping amplitude (default: 1.0)')
    parser.add_argument('--u', type=float, default=4.0, help='on-site interaction (default: 4.0)')
    parser.add_argument('--mu', type=float, default=0.0, help='chemical potential (default: 0.0)')
    parser.add_argument('--boundary', type=str, default='pbc', choices=['pbc'],
                        help='boundary condition (must be pbc for k-symmetry)')
    parser.add_argument('--spectrum', action='store_true', help='print energy spectrum for each k')
    parser.add_argument('--save', type=str, default=none, help='save results to file')
    parser.add_argument('--k', type=int, default=none,
                        help='calculate only specific k-sector (default: all sectors)')

    args = parser.parse_args()

    print("="*70)
    print("1d hubbard model - exact diagonalization with translation symmetry")
    print("="*70)

    print(f"\n=== system parameters ===")
    print(f"sites: l = {args.l}")
    print(f"hopping: t = {args.t}")
    print(f"interaction: u = {args.u}")
    print(f"chemical potential: μ = {args.mu}")
    print(f"boundary: {args.boundary}")

    # determine which k-sectors to calculate
    if args.k is not none:
        k_sectors = [args.k]
        print(f"calculating only k = {args.k} sector")
    else:
        k_sectors = range(args.l)
        print(f"calculating all k-sectors: 0, 1, ..., {args.l-1}")

    # store results for all k-sectors
    all_results = []
    all_energies = []
    all_wavefunctions = []

    # loop over momentum sectors
    for k in k_sectors:
        print(f"\n{'='*70}")
        print(f"momentum sector k = {k} (momentum = 2πk/l = {2*np.pi*k/args.l:.4f})")
        print(f"{'='*70}")

        # build hamiltonian for this k-sector
        h, basis = build_hubbard_hamiltonian_k(args.l, k, args.t, args.u, args.mu, args.boundary)

        print(f"hilbert space dimension in k={k} sector: {basis.ns}")

        # calculate observables
        results, e, psi0 = calculate_observables_k(h, basis, args.l, k)

        print(f"\n=== results for k = {k} ===")
        print(f"ground state energy: e0(k={k}) = {results['e0']:.10f}")
        if not np.isnan(results['e1']):
            print(f"first excited state: e1(k={k}) = {results['e1']:.10f}")
            print(f"energy gap: δe(k={k}) = {results['gap']:.10f}")
        print(f"<n_up> = {results['n_up']:.10f}")
        print(f"<n_down> = {results['n_down']:.10f}")
        print(f"<n_total> = {results['n_total']:.10f}")
        print(f"double occupancy <d> = {results['double_occupancy']:.10f}")

        # print energy spectrum if requested
        if args.spectrum:
            print(f"\n=== energy spectrum for k={k} (first 10 states) ===")
            for i in range(min(10, len(e))):
                print(f"e[{i}](k={k}) = {e[i]:.10f}")

        all_results.append(results)
        all_energies.append(e)
        all_wavefunctions.append(psi0)

    # find global ground state across all k-sectors
    print(f"\n{'='*70}")
    print("global ground state")
    print(f"{'='*70}")

    ground_state_k = min(range(len(all_results)), key=lambda i: all_results[i]['e0'])
    gs_results = all_results[ground_state_k]

    print(f"\nglobal ground state is in momentum sector k = {gs_results['k']}")
    print(f"ground state energy: e0 = {gs_results['e0']:.10f}")
    print(f"ground state momentum: p = 2πk/l = {2*np.pi*gs_results['k']/args.l:.4f}")
    print(f"<n_up> = {gs_results['n_up']:.10f}")
    print(f"<n_down> = {gs_results['n_down']:.10f}")
    print(f"<n_total> = {gs_results['n_total']:.10f}")
    print(f"double occupancy <d> = {gs_results['double_occupancy']:.10f}")

    # print energy vs momentum
    print(f"\n=== ground state energy vs momentum ===")
    print(f"{'k':<5} {'2πk/l':<12} {'e0(k)':<18} {'hilbert dim':<12}")
    print("-"*50)
    for res in all_results:
        momentum = 2*np.pi*res['k']/args.l
        marker = " <-- ground state" if res['k'] == gs_results['k'] else ""
        print(f"{res['k']:<5} {momentum:<12.6f} {res['e0']:<18.10f} {res['hilbert_dim']:<12}{marker}")

    # save results to file
    if args.save:
        with open(args.save, 'w') as f:
            f.write("# 1d hubbard model - exact diagonalization with translation symmetry\n")
            f.write(f"# l={args.l}, t={args.t}, u={args.u}, mu={args.mu}, boundary={args.boundary}\n")
            f.write("#\n")
            f.write(f"# global ground state in k={gs_results['k']} sector\n")
            f.write(f"e0 = {gs_results['e0']:.15f}\n")
            f.write(f"ground_state_k = {gs_results['k']}\n")
            f.write(f"ground_state_momentum = {2*np.pi*gs_results['k']/args.l:.15f}\n")
            f.write(f"gap = {gs_results['gap']:.15f}\n")
            f.write(f"n_total = {gs_results['n_total']:.15f}\n")
            f.write(f"n_up = {gs_results['n_up']:.15f}\n")
            f.write(f"n_down = {gs_results['n_down']:.15f}\n")
            f.write(f"double_occupancy = {gs_results['double_occupancy']:.15f}\n")
            f.write("#\n")
            f.write("# energy vs momentum for all k-sectors:\n")
            f.write("# k, momentum, e0(k), hilbert_dim\n")
            for res in all_results:
                momentum = 2*np.pi*res['k']/args.l
                f.write(f"{res['k']}, {momentum:.15f}, {res['e0']:.15f}, {res['hilbert_dim']}\n")
        print(f"\nresults saved to {args.save}")

    print("\n" + "="*70)
    print("calculation complete!")
    print("="*70)

    return all_results, gs_results


if __name__ == "__main__":
    main()
