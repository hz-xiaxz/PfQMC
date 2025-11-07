# Exact Diagonalization for 1D Hubbard Model

This directory contains exact diagonalization (ED) code for the 1D Hubbard model using the QuSpin library. ED provides exact results for small systems and is useful for benchmarking QMC.

## Prerequisites

Install QuSpin:

```bash
pip install quspin
```

Or with conda:

```bash
conda install -c weinbe58 quspin
```

## Usage

### Basic ED Calculation

```bash
cd ED
python hubbard_ed.py --L 4 --U 4.0 --mu 0.0 --boundary OBC
```

### Parameters

- `--L`: Number of lattice sites (default: 4, max ~10 for reasonable runtime)
- `--t`: Hopping amplitude (default: 1.0)
- `--U`: On-site interaction strength (default: 4.0)
- `--mu`: Chemical potential (default: 0.0)
- `--boundary`: Boundary condition, `OBC` or `PBC` (default: OBC)
- `--spectrum`: Print energy spectrum
- `--save`: Save results to file

### Example

```bash
# Half-filled 4-site chain with U=4
python hubbard_ed.py --L 4 --U 4.0 --mu 0.0 --boundary OBC

# With energy spectrum
python hubbard_ed.py --L 4 --U 4.0 --spectrum

# Save results
python hubbard_ed.py --L 6 --U 2.0 --save results_L6_U2.txt
```

## Comparing QMC and ED

Use the comparison script to benchmark QMC against exact ED:

```bash
cd ED
python compare_qmc_ed.py --L 4 --U 4.0 --mu 0.0
```

This will:
1. Run ED to get exact ground state
2. Run QMC simulation
3. Compare all observables (energy, particle number, structure factors)
4. Report agreement and check for sign problem

## Quick Test

Test your QMC results against ED:

```bash
# Make scripts executable
chmod +x hubbard_ed.py compare_qmc_ed.py

# Run comparison for L=4, U=4.0
python compare_qmc_ed.py --L 4 --U 4.0 --evaluationLength 5000
```

Expected output shows QMC vs ED agreement with relative error < 0.1%.

## Size Limitations

ED scales exponentially with system size:

| L (sites) | Hilbert space | Memory    | Time       |
|-----------|---------------|-----------|------------|
| 2         | 16            | Instant   | < 1s       |
| 4         | 256           | Instant   | < 1s       |
| 6         | 4,096         | ~130 KB   | ~1s        |
| 8         | 65,536        | ~33 MB    | ~10s       |
| 10        | 1,048,576     | ~8 GB     | ~5 min     |

**Recommendation**: Use ED for L d 8. Use QMC for larger systems.

## Observables

1. **Energy**: Ground state energy E€
2. **Particle numbers**: èN‘é, èN“é, èNé
3. **Spin structure factor S(À)**: Antiferromagnetic correlations
4. **Charge structure factor N(À)**: Charge density wave
5. **Double occupancy**: èn‘ n“é
6. **NN spin correlation**: èS·Sé

## References

- **QuSpin**: Weinberg & Bukov, SciPost Phys. 2, 003 (2017)
- **Hubbard Model**: Essler et al., "The One-Dimensional Hubbard Model" (2005)
