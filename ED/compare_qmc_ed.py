#!/usr/bin/env python3
"""
Compare QMC and ED results for 1D Hubbard model

This script runs both QMC and ED for the same parameters and compares results.
"""

import numpy as np
import subprocess
import argparse
import os
import re
import sys

def run_qmc(L, LTau, dt, U, mu, hsScheme, boundary, nthreads, nseed, evaluationLength, qmc_executable):
    """
    Run QMC simulation

    Returns:
        dict: Parsed QMC results
    """
    filepath = "./output/"
    os.makedirs(filepath, exist_ok=True)

    boundary_flag = 1 if boundary == 'OBC' else 0

    # Source Intel oneAPI environment and run QMC
    qmc_cmd = [
        qmc_executable, "--hubbard",
        str(L), str(LTau), str(dt), str(U),
        str(nthreads), str(nseed), str(evaluationLength),
        filepath, str(mu), str(hsScheme), str(boundary_flag)
    ]

    # Wrap command with Intel oneAPI environment sourcing
    cmd = f"source /opt/intel/oneapi/setvars.sh --force > /dev/null 2>&1 && {' '.join(qmc_cmd)}"

    print("\n" + "="*60)
    print("Running QMC simulation...")
    print("Command:", cmd)
    print("="*60)

    result = subprocess.run(cmd, capture_output=True, text=True, shell=True, executable='/bin/bash')

    # Print QMC output for debugging
    print("\n--- QMC stdout ---")
    print(result.stdout)
    if result.stderr:
        print("\n--- QMC stderr ---")
        print(result.stderr)
    print("--- End QMC output ---\n")

    # Parse output
    qmc_results = {}
    for line in result.stdout.split('\n'):
        if 'AveEnergy' in line:
            # Parse: AveEnergy = (1.95236,1.34772e-12) AveSign = (1,8.74162e-10) ...
            energy_match = re.search(r'AveEnergy = \(([-\d.e]+)', line)
            sign_match = re.search(r'AveSign = \(([-\d.e]+)', line)
            n_match = re.search(r'AveN = \(([-\d.e]+)', line)
            s_pi_match = re.search(r'AveS\(π\) = \(([-\d.e]+)', line)
            n_pi_match = re.search(r'AveN\(π\) = \(([-\d.e]+)', line)

            if energy_match:
                qmc_results['energy'] = float(energy_match.group(1))
            if sign_match:
                qmc_results['sign'] = float(sign_match.group(1))
            if n_match:
                qmc_results['N_total'] = float(n_match.group(1))
            if s_pi_match:
                qmc_results['S_pi'] = float(s_pi_match.group(1))
            if n_pi_match:
                qmc_results['N_pi'] = float(n_pi_match.group(1))

    print("\n--- Parsed QMC results ---")
    for key, val in qmc_results.items():
        print(f"  {key}: {val}")
    print("--- End parsed QMC results ---\n")

    return qmc_results


def run_ed(L, U, mu, boundary, ed_script):
    """
    Run ED calculation

    Returns:
        dict: Parsed ED results
    """
    boundary_flag = boundary

    # Activate venv and run ED script
    # Get absolute path to ED script directory
    ed_script_abs = os.path.abspath(ed_script)
    ed_dir = os.path.dirname(ed_script_abs)
    venv_activate = os.path.join(ed_dir, ".venv/bin/activate")

    cmd = f"source {venv_activate} && python3 {ed_script_abs} --L {L} --U {U} --mu {mu} --boundary {boundary_flag}"

    print("\n" + "="*60)
    print("Running ED calculation...")
    print("Command:", cmd)
    print("="*60)

    result = subprocess.run(cmd, capture_output=True, text=True, shell=True, executable='/bin/bash')
    print(result.stdout)
    if result.stderr:
        print("--- ED stderr ---")
        print(result.stderr)

    # Parse output
    ed_results = {}
    for line in result.stdout.split('\n'):
        if 'Ground state energy: E0 =' in line:
            ed_results['energy'] = float(line.split('=')[1].strip())
        elif '<N_total> =' in line:
            ed_results['N_total'] = float(line.split('=')[1].strip())
        elif 'Spin structure factor S(π) =' in line:
            ed_results['S_pi'] = float(line.split('=')[1].strip())
        elif 'Charge structure factor N(π) =' in line:
            ed_results['N_pi'] = float(line.split('=')[1].strip())
        elif 'Double occupancy <D> =' in line:
            ed_results['double_occ'] = float(line.split('=')[1].strip())
        elif 'NN spin correlation' in line:
            ed_results['nn_spin_corr'] = float(line.split('=')[1].strip())

    print("\n--- Parsed ED results ---")
    for key, val in ed_results.items():
        print(f"  {key}: {val}")
    print("--- End parsed ED results ---\n")

    return ed_results


def compare_results(qmc_results, ed_results):
    """
    Compare QMC and ED results
    """
    print("\n" + "="*60)
    print("COMPARISON: QMC vs ED")
    print("="*60)

    print(f"\n{'Observable':<20} {'QMC':<20} {'ED (Exact)':<20} {'Diff':<15} {'Rel. Error'}")
    print("-"*90)

    for key in ['energy', 'N_total', 'S_pi', 'N_pi']:
        if key in qmc_results and key in ed_results:
            qmc_val = qmc_results[key]
            ed_val = ed_results[key]
            diff = qmc_val - ed_val
            rel_err = abs(diff / ed_val) * 100 if ed_val != 0 else 0

            print(f"{key:<20} {qmc_val:<20.10f} {ed_val:<20.10f} {diff:<15.2e} {rel_err:.4f}%")

    if 'sign' in qmc_results:
        print(f"\nQMC average sign: {qmc_results['sign']:.10f}")
        if qmc_results['sign'] < 0.5:
            print("⚠️  WARNING: Severe sign problem! QMC results may be unreliable.")
        elif qmc_results['sign'] < 0.9:
            print("⚠️  WARNING: Moderate sign problem. Consider longer runs.")
        else:
            print("✓ Good sign - QMC is reliable")

    print("\n" + "="*60)


def main():
    parser = argparse.ArgumentParser(description='Compare QMC and ED for 1D Hubbard model')
    parser.add_argument('--L', type=int, default=4, help='Number of sites (default: 4, max ~10 for ED)')
    parser.add_argument('--U', type=float, default=4.0, help='Interaction strength (default: 4.0)')
    parser.add_argument('--mu', type=float, default=0.0, help='Chemical potential (default: 0.0)')
    parser.add_argument('--boundary', type=str, default='OBC', choices=['OBC', 'PBC'],
                        help='Boundary condition (default: OBC)')

    # QMC parameters
    parser.add_argument('--LTau', type=int, default=200, help='QMC: Number of time slices (default: 200)')
    parser.add_argument('--dt', type=float, default=0.1, help='QMC: Time step (default: 0.1)')
    parser.add_argument('--hsScheme', type=int, default=0, help='QMC: HS scheme (default: 0)')
    parser.add_argument('--nthreads', type=int, default=8, help='QMC: OpenMP threads (default: 8)')
    parser.add_argument('--nseed', type=int, default=42, help='QMC: Random seed (default: 42)')
    parser.add_argument('--evaluationLength', type=int, default=5000, help='QMC: Measurements (default: 5000)')

    # Executable paths
    parser.add_argument('--qmc-exec', type=str, default='../build/main',
                        help='Path to QMC executable (default: ../build/main)')
    parser.add_argument('--ed-script', type=str, default='./hubbard_ed.py',
                        help='Path to ED script (default: ./hubbard_ed.py)')

    args = parser.parse_args()

    if args.L > 10:
        print("WARNING: ED is exponentially expensive. L > 10 may take very long time!")
        response = input("Continue? (y/n): ")
        if response.lower() != 'y':
            sys.exit(0)

    print("="*60)
    print("QMC vs ED Comparison for 1D Hubbard Model")
    print("="*60)
    print(f"Parameters: L={args.L}, U={args.U}, μ={args.mu}, boundary={args.boundary}")
    print(f"QMC: β={args.LTau * args.dt}, dt={args.dt}, measurements={args.evaluationLength}")

    # Run ED (exact, fast for small systems)
    ed_results = run_ed(args.L, args.U, args.mu, args.boundary, args.ed_script)

    # Run QMC
    qmc_results = run_qmc(
        args.L, args.LTau, args.dt, args.U, args.mu,
        args.hsScheme, args.boundary, args.nthreads,
        args.nseed, args.evaluationLength, args.qmc_exec
    )

    # Compare
    compare_results(qmc_results, ed_results)


if __name__ == "__main__":
    main()
