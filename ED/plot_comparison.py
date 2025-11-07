#!/usr/bin/env python3
"""
Plot QMC vs ED energy comparison for 1D Hubbard model

This script generates comparison plots between QMC and ED results.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import subprocess
import re
import os


def run_comparison(L, U_values, mu, boundary, LTau, dt, hsScheme, nthreads, nseed, evaluationLength, qmc_exec, ed_script):
    """
    Run QMC and ED for multiple U values and collect results

    Returns:
        dict: Results for each U value
    """
    results = {}

    for U in U_values:
        print(f"\n{'='*60}")
        print(f"Running L={L}, U={U}")
        print(f"{'='*60}")

        # Run ED
        boundary_flag = boundary
        ed_cmd = f"source {os.path.dirname(os.path.abspath(ed_script))}/.venv/bin/activate && python3 {ed_script} --L {L} --U {U} --mu {mu} --boundary {boundary_flag}"
        ed_result = subprocess.run(ed_cmd, capture_output=True, text=True, shell=True, executable='/bin/bash')

        ed_energy = None
        for line in ed_result.stdout.split('\n'):
            if 'Ground state energy: E0 =' in line:
                ed_energy = float(line.split('=')[1].strip())
                break

        # Run QMC
        filepath = "./output/"
        os.makedirs(filepath, exist_ok=True)
        boundary_int = 1 if boundary == 'OBC' else 0

        qmc_cmd = f"source /opt/intel/oneapi/setvars.sh --force > /dev/null 2>&1 && {qmc_exec} --hubbard {L} {LTau} {dt} {U} {nthreads} {nseed} {evaluationLength} {filepath} {mu} {hsScheme} {boundary_int}"
        qmc_result = subprocess.run(qmc_cmd, capture_output=True, text=True, shell=True, executable='/bin/bash')

        qmc_energy = None
        qmc_sign = None
        for line in qmc_result.stdout.split('\n'):
            if 'AveEnergy' in line:
                energy_match = re.search(r'AveEnergy = \(([-\d.e]+)', line)
                sign_match = re.search(r'AveSign = \(([-\d.e]+)', line)
                if energy_match:
                    qmc_energy = float(energy_match.group(1))
                if sign_match:
                    qmc_sign = float(sign_match.group(1))
                break

        results[U] = {
            'ed_energy': ed_energy,
            'qmc_energy': qmc_energy,
            'qmc_sign': qmc_sign
        }

        print(f"ED Energy: {ed_energy:.6f}")
        print(f"QMC Energy: {qmc_energy:.6f}")
        print(f"QMC Sign: {qmc_sign:.6f}")

    return results


def plot_energy_comparison(results, L, mu, boundary, output_file='energy_comparison.png'):
    """
    Create comparison plots
    """
    U_values = sorted(results.keys())
    ed_energies = [results[U]['ed_energy'] for U in U_values]
    qmc_energies = [results[U]['qmc_energy'] for U in U_values]
    qmc_signs = [results[U]['qmc_sign'] for U in U_values]

    # Calculate errors
    errors = [qmc - ed for qmc, ed in zip(qmc_energies, ed_energies)]
    rel_errors = [abs(err/ed)*100 if ed != 0 else 0 for err, ed in zip(errors, ed_energies)]

    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'QMC vs ED: 1D Hubbard Model (L={L}, μ={mu}, {boundary})', fontsize=14, fontweight='bold')

    # Plot 1: Energy comparison
    ax1 = axes[0, 0]
    ax1.plot(U_values, ed_energies, 'o-', label='ED (Exact)', markersize=8, linewidth=2)
    ax1.plot(U_values, qmc_energies, 's--', label='QMC', markersize=8, linewidth=2)
    ax1.set_xlabel('U/t', fontsize=12)
    ax1.set_ylabel('Ground State Energy', fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_title('Energy Comparison')

    # Plot 2: Absolute error
    ax2 = axes[0, 1]
    ax2.plot(U_values, errors, 'o-', color='red', markersize=8, linewidth=2)
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax2.set_xlabel('U/t', fontsize=12)
    ax2.set_ylabel('Error (QMC - ED)', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.set_title('Absolute Error')

    # Plot 3: Relative error
    ax3 = axes[1, 0]
    ax3.plot(U_values, rel_errors, 'o-', color='orange', markersize=8, linewidth=2)
    ax3.set_xlabel('U/t', fontsize=12)
    ax3.set_ylabel('Relative Error (%)', fontsize=12)
    ax3.grid(True, alpha=0.3)
    ax3.set_title('Relative Error')
    ax3.set_yscale('log')

    # Plot 4: QMC Sign
    ax4 = axes[1, 1]
    ax4.plot(U_values, qmc_signs, 'o-', color='green', markersize=8, linewidth=2)
    ax4.axhline(y=1.0, color='k', linestyle='--', alpha=0.3)
    ax4.axhline(y=0.9, color='orange', linestyle='--', alpha=0.3, label='Warning threshold')
    ax4.set_xlabel('U/t', fontsize=12)
    ax4.set_ylabel('Average Sign', fontsize=12)
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)
    ax4.set_title('QMC Sign Problem')
    ax4.set_ylim([0, 1.1])

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Plot saved to {output_file}")

    # Print summary table
    print(f"\n{'='*80}")
    print(f"{'U/t':<8} {'ED Energy':<15} {'QMC Energy':<15} {'Abs Error':<12} {'Rel Error (%)':<15} {'Sign':<8}")
    print(f"{'='*80}")
    for U, ed_e, qmc_e, err, rel_err, sign in zip(U_values, ed_energies, qmc_energies, errors, rel_errors, qmc_signs):
        print(f"{U:<8.2f} {ed_e:<15.6f} {qmc_e:<15.6f} {err:<12.2e} {rel_err:<15.4f} {sign:<8.6f}")
    print(f"{'='*80}")


def main():
    parser = argparse.ArgumentParser(description='Plot QMC vs ED energy comparison')
    parser.add_argument('--L', type=int, default=4, help='Number of sites (default: 4)')
    parser.add_argument('--U-min', type=float, default=0.0, help='Minimum U value (default: 0.0)')
    parser.add_argument('--U-max', type=float, default=8.0, help='Maximum U value (default: 8.0)')
    parser.add_argument('--U-step', type=float, default=2.0, help='U step size (default: 2.0)')
    parser.add_argument('--mu', type=float, default=0.0, help='Chemical potential (default: 0.0)')
    parser.add_argument('--boundary', type=str, default='OBC', choices=['OBC', 'PBC'],
                        help='Boundary condition (default: OBC)')

    # QMC parameters
    parser.add_argument('--LTau', type=int, default=200, help='QMC: Time slices (default: 200)')
    parser.add_argument('--dt', type=float, default=0.1, help='QMC: Time step (default: 0.1)')
    parser.add_argument('--hsScheme', type=int, default=0, help='QMC: HS scheme (default: 0)')
    parser.add_argument('--nthreads', type=int, default=8, help='QMC: Threads (default: 8)')
    parser.add_argument('--nseed', type=int, default=42, help='QMC: Random seed (default: 42)')
    parser.add_argument('--evaluationLength', type=int, default=20000, help='QMC: Measurements (default: 20000)')

    # Executable paths
    parser.add_argument('--qmc-exec', type=str, default='../build/main',
                        help='Path to QMC executable (default: ../build/main)')
    parser.add_argument('--ed-script', type=str, default='./hubbard_ed.py',
                        help='Path to ED script (default: ./hubbard_ed.py)')

    # Output
    parser.add_argument('--output', type=str, default='energy_comparison.png',
                        help='Output plot filename (default: energy_comparison.png)')

    args = parser.parse_args()

    # Generate U values
    U_values = np.arange(args.U_min, args.U_max + args.U_step/2, args.U_step)

    print("="*60)
    print("QMC vs ED Energy Comparison Plot Generator")
    print("="*60)
    print(f"Parameters: L={args.L}, μ={args.mu}, boundary={args.boundary}")
    print(f"U values: {U_values}")
    print(f"QMC: β={args.LTau * args.dt}, dt={args.dt}, measurements={args.evaluationLength}")

    # Run comparisons
    results = run_comparison(
        args.L, U_values, args.mu, args.boundary,
        args.LTau, args.dt, args.hsScheme, args.nthreads,
        args.nseed, args.evaluationLength, args.qmc_exec, args.ed_script
    )

    # Create plots
    plot_energy_comparison(results, args.L, args.mu, args.boundary, args.output)


if __name__ == "__main__":
    main()
