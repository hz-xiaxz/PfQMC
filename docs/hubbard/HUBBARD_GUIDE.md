# 1D Hubbard Model Implementation - Usage Guide

## Overview

This implementation adds the **1D Hubbard model** in Majorana fermion representation to the PfQMC codebase.

### Hamiltonian

```
H = -t Σ_{i,σ} (c†_{i,σ} c_{i+1,σ} + h.c.) + U Σ_i n_{i,↑} n_{i,↓} - μ Σ_{i,σ} n_{i,σ}
```

where:
- `t = 1` (hopping amplitude, set to 1 as energy unit)
- `U` = on-site interaction strength (repulsive if U > 0, attractive if U < 0)
- `μ` = chemical potential
- `σ ∈ {↑, ↓}` = spin index

### Majorana Representation

Each lattice site has **4 Majorana modes**:
- `γ^1_{i,↑}, γ^2_{i,↑}` for spin-up
- `γ^1_{i,↓}, γ^2_{i,↓}` for spin-down

System dimension: `nDim = 4L` (where L is the number of sites)

## Building

Make sure you've already built the main PfQMC code:

```bash
# From the PfQMC directory
cd build
cmake .. -DEIGEN3_INCLUDE_DIR=/path/to/eigen3
make
```

## Running Simulations

### Command Line Syntax

```bash
mpirun -np <num_processes> ./build/main --hubbard Lx LTau dt U nthreads nseed evaluationLength filepath mu hsScheme boundary
```

### Parameters

- `Lx`: Number of lattice sites (e.g., 4, 6, 8, 10)
- `LTau`: Number of imaginary time slices (e.g., 100, 200)
- `dt`: Imaginary time step (e.g., 0.1, 0.05)
- `U`: On-site interaction strength (e.g., 1.0, 2.0, 4.0)
- `nthreads`: Number of OpenMP threads per MPI process (e.g., 8)
- `nseed`: Random seed for reproducibility
- `evaluationLength`: Number of measurement sweeps (e.g., 1000, 5000)
- `filepath`: Output directory path (must end with `/`)
- `mu`: Chemical potential (e.g., 0.0 for half-filling)
- `hsScheme`: Hubbard-Stratonovich scheme (0 or 1, try both)
- `boundary`: Boundary condition (0 = PBC, 1 = OBC)

### Example Commands

#### 1. Half-filled system (L=8, U=4)

```bash
mpirun -np 4 ./build/main --hubbard 8 200 0.1 4.0 8 42 5000 ./output/ 0.0 0 1
```

This runs:
- 8-site chain
- 200 time slices with dt=0.1 (β = 20)
- U = 4.0 (moderate interaction)
- Half-filling (μ = 0.0)
- OBC (boundary = 1)
- 4 MPI processes with different random seeds

#### 2. Attractive Hubbard model (U < 0)

```bash
mpirun -np 4 ./build/main --hubbard 6 150 0.1 -2.0 8 123 3000 ./output/ 0.0 0 0
```

This runs:
- 6-site chain
- U = -2.0 (attractive interaction, favors pairing)
- PBC (boundary = 0)
- Better sign behavior expected

#### 3. Doped system (away from half-filling)

```bash
mpirun -np 4 ./build/main --hubbard 8 200 0.1 4.0 8 42 5000 ./output/ 1.0 0 1
```

This adds chemical potential μ = 1.0 to shift away from half-filling.

## Output

### Output Files

Each MPI process generates a file:
```
<filepath>hubbard-<Lx>-<LTau>-<dt>-<U>-<process_id>-<seed>-<mu>-BC=<boundary>.out
```

Example: `output/hubbard-8-200-0.10-4.00-0-12345-0.00-BC=1.out`

### Output Format

Each line during measurement contains:
```
iter = <i> sign = <sign> energy = <E> N = <particle_number> S(π) = <spin_SF> N(π) = <charge_SF>
```

At the end:
```
AveEnergy = <E_avg> AveSign = <sign_avg> AveN = <N_avg> AveS(π) = <S(π)_avg> AveN(π) = <N(π)_avg>
```

### Observables

1. **Energy**: Total energy `<H>`
2. **Particle Number**: Total number of particles `<N> = Σ_{i,σ} <n_{i,σ}>`
3. **Spin Structure Factor S(π)**: Measures antiferromagnetic correlations
   - `S(q) = (1/L) Σ_{i,j} e^{iq(i-j)} <S^z_i S^z_j>`
   - Large S(π) indicates AFM order
4. **Charge Structure Factor N(π)**: Measures charge density wave
   - `N(q) = (1/L) Σ_{i,j} e^{iq(i-j)} <n_i n_j>`
   - Large N(π) indicates CDW order

## Physical Regimes

### 1. Weak Coupling (U/t ≪ 1)
- Nearly free fermions
- Small sign problem
- Good convergence

**Example:**
```bash
mpirun -np 4 ./build/main --hubbard 8 150 0.1 0.5 8 42 3000 ./output/ 0.0 0 1
```

### 2. Moderate Coupling (U/t ~ 1-4)
- Interesting physics: AFM vs metallic
- Moderate sign problem
- May need longer runs

**Example:**
```bash
mpirun -np 4 ./build/main --hubbard 8 200 0.1 2.0 8 42 5000 ./output/ 0.0 0 1
```

### 3. Strong Coupling (U/t ≫ 4)
- Mott insulating behavior at half-filling
- Strong AFM correlations
- **Warning**: May have severe sign problem!
- Try hsScheme=1 if sign problem occurs

**Example:**
```bash
mpirun -np 8 ./build/main --hubbard 6 250 0.05 8.0 8 42 10000 ./output/ 0.0 1 1
```

### 4. Attractive Interaction (U < 0)
- Favors Cooper pairing
- Better sign behavior
- Superconducting correlations

**Example:**
```bash
mpirun -np 4 ./build/main --hubbard 8 200 0.1 -2.0 8 42 5000 ./output/ 0.0 0 0
```

## Sign Problem

The fermion sign problem is a major challenge in QMC. Monitor the average sign:
- **Sign ≈ 1**: Good, results are reliable
- **Sign < 0.5**: Significant sign problem, need more statistics
- **Sign ≈ 0**: Severe sign problem, results may be unreliable

### Strategies for Sign Problem

1. **Try both HS schemes** (hsScheme = 0 or 1):
   - Different schemes can have different sign structures
   - Try both and use whichever gives better sign

2. **Reduce temperature** (increase β = LTau * dt):
   - Lower T can sometimes improve sign
   - But may need smaller dt

3. **Attractive U** (U < 0):
   - Generally has better sign behavior
   - Related to repulsive U by particle-hole transformation

4. **Smaller system size**:
   - Sign problem often grows exponentially with size
   - Test on small systems first

## Convergence Checks

### 1. Trotter Error
The time discretization has error O(dt²). Check convergence by reducing dt:

```bash
# dt = 0.1
mpirun -np 4 ./build/main --hubbard 8 200 0.1 4.0 8 42 3000 ./output/ 0.0 0 1

# dt = 0.05 (β kept constant, so LTau doubled)
mpirun -np 4 ./build/main --hubbard 8 400 0.05 4.0 8 42 3000 ./output/ 0.0 0 1
```

Results should agree within error bars.

### 2. Finite Size Effects
Check dependence on system size:

```bash
for L in 4 6 8 10; do
    mpirun -np 4 ./build/main --hubbard $L 200 0.1 4.0 8 42 3000 ./output/ 0.0 0 1
done
```

### 3. Statistical Errors
Run with different random seeds and average:

```bash
for seed in 100 200 300 400 500; do
    mpirun -np 4 ./build/main --hubbard 8 200 0.1 4.0 8 $seed 3000 ./output/ 0.0 0 1
done
```

## Debugging

### Common Issues

1. **Sign drops rapidly during thermalization**
   - Check if U is too large
   - Try different hsScheme
   - Check that β (= LTau * dt) is not too large

2. **Energy seems wrong**
   - Check with exact diagonalization for small systems
   - Verify half-filling: should have N ≈ L at μ = 0

3. **NaN or Inf in output**
   - Reduce dt (stabilization issue)
   - Check stabilizationTime parameter (currently 10)
   - May need to adjust UDT threshold

### Testing on Small Systems

For L=2, you can compare with exact diagonalization:

```bash
# Small system for testing
mpirun -np 2 ./build/main --hubbard 2 100 0.1 2.0 8 42 1000 ./output/ 0.0 0 1
```

Expected at half-filling (N=2):
- Energy should be smooth function of U
- Sign should be close to 1 for U ≤ 4

## Performance Tips

1. **MPI parallelization**: Run multiple independent simulations
   ```bash
   mpirun -np 16 ./build/main --hubbard ...
   ```
   Each process uses different random seed.

2. **OpenMP threads**: Use CPU cores efficiently
   ```bash
   # On a 16-core node, try:
   # 4 MPI × 4 threads, or 8 MPI × 2 threads
   ```

3. **Batch jobs**: For production runs, use job scheduler:
   ```bash
   #!/bin/bash
   #SBATCH -n 32
   #SBATCH -t 24:00:00
   mpirun ./build/main --hubbard 10 300 0.05 4.0 8 42 10000 ./output/ 0.0 0 1
   ```

## References

1. **Hubbard Model**: Hubbard, J. (1963). Proc. R. Soc. Lond. A 276, 238.
2. **PfQMC Algorithm**: [arXiv:2408.10311](https://arxiv.org/abs/2408.10311)
3. **Majorana Fermions**: See CLAUDE.md for detailed explanation
4. **Sign Problem**: Loh et al., Phys. Rev. B 41, 9301 (1990)

## Contact

For issues or questions about the Hubbard model implementation, refer to the main PfQMC repository.
