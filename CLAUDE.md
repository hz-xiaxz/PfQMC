# Pfaffian Quantum Monte Carlo (PfQMC) Codebase Documentation

This document provides a comprehensive overview of the PfQMC codebase structure and implementation details.

## Overview

This is a prototype implementation of the **Pfaffian Quantum Monte Carlo (PfQMC)** algorithm for simulating fermionic quantum many-body systems in Majorana representation. The code implements several models:

- **Spinless t-V model** on square and honeycomb lattices
- **Kitaev chain** (p-wave superconductor chain model)
- **Single Majorana honeycomb model**
- **1D Hubbard model** with spin-1/2 fermions (NEW)

Reference: [arXiv:2408.10311](https://arxiv.org/abs/2408.10311)

## Core Algorithm

### Pfaffian QMC Method

The algorithm uses:
1. **Majorana fermion representation**: Decomposes complex fermions into pairs of Majorana fermions
2. **Imaginary time path integral**: Discretized with Trotter decomposition (time step `dt`, `LTau` slices)
3. **Hubbard-Stratonovich (HS) transformation**: Decouples quartic interaction terms using auxiliary Ising fields
4. **Metropolis sampling**: Updates auxiliary fields to sample configuration space
5. **UDT stabilization**: Numerically stable matrix decomposition for Green's function updates
6. **Pfaffian computation**: Uses PFAPACK library for antisymmetric matrix Pfaffians

### Sign Problem

The code monitors the "sign" which indicates whether the sign-problem-free condition is maintained. Several observables check this via `pfqmc.sign` and compare with `pfqmc.getSignRaw()`.

### Sign Computation Details

The sign represents the fermion sign problem and should be ±1 for sign-problem-free simulations. The code uses two complementary tracking methods:

**1. Accumulated Sign** (pfqmc.cpp:29,40-41,90-91)
- Initialized at construction: `sign = getSignRaw()`
- Updated multiplicatively during each sweep:
  ```cpp
  signCur = op_array[l]->update(g);  // Get sign from operator update
  this->sign *= signCur;              // Accumulate into total sign
  ```
- Fast but accumulates numerical errors over many sweeps

**2. Raw Sign Calculation** (pfqmc.cpp:125-192)
- `getSignRaw()` recomputes sign from scratch using Pfaffians
- More expensive but numerically exact
- Used for initialization and periodic verification

**Sign Computation Algorithm:**
```cpp
DataType PfQMC::getSignRaw() {
    const DataType extraSign = ((nDim / 2) % 2 == 0) ? 1.0 : -1.0;

    signCur = op_array[0]->getSignOfWeight();  // Initialize with first operator

    for (int i = 1; i < op_length; i++) {
        signNext = op_array[i]->getSignOfWeight();

        // Compute sign from Pfaffian of 4n×4n block matrix
        signPfaf = pfaffianForSignOfProduct(gNext, gCur);

        // Accumulate all signs with parity factor
        signCur = signCur * signNext * signPfaf * extraSign;
    }
    return signCur;
}
```

**Pfaffian Block Matrix** (skewMatUtils.cpp:223-264)

The key function `pfaffianForSignOfProduct(G1, G2)` computes the sign of a Pfaffian for a 4n×4n antisymmetric block matrix:

```cpp
// Build antisymmetric matrix:
// A = [ G1    -I  ]
//     [ I     G2  ]
MatType A = MatType::Zero(4*n, 4*n);
A.block(0, 0, 2*n, 2*n) = G1;
A.block(2*n, 2*n, 2*n, 2*n) = G2;
A.block(0, 2*n, 2*n, 2*n) = -Identity;
A.block(2*n, 0, 2*n, 2*n) = Identity;

// Compute sign using PFAPACK's sktrf decomposition
DataType r = signOfPfaf(A);
return r;
```

**Sign Verification Protocol** (main.cpp:416-424)

Every 20 iterations during measurement phase:
```cpp
if (i % 20 == 0) {
    signRaw = pfqmc.getSignRaw();   // Recompute from scratch
    if (std::abs(sign - signRaw) > threshold) {
        fixSign(sign, signRaw, threshold);  // Correct accumulated sign
        pfqmc.sign = sign;
    }
}
```

**Key Properties:**
- Sign computation is **model-independent** - handled entirely by PfQMC class
- For Hubbard model: `nDim = 4L` (4 Majorana modes per site)
- Block structure distinguishes Majorana QMC from conventional DQMC
- Periodic verification catches accumulated numerical errors
- Sign should remain ±1 for sign-problem-free cases

## Directory Structure

```
PfQMC/
├── main.cpp                    # Entry point with different model runners
├── inc/                        # Header files
│   ├── types.h                 # Type definitions and basic utilities
│   ├── operator.h              # Base Operator class hierarchy
│   ├── pfqmc.h                 # Main PfQMC algorithm class
│   ├── qr_udt.h                # UDT decomposition for stabilization
│   ├── skewMatUtils.h          # Pfaffian and antisymmetric matrix utilities
│   ├── spinless_tV.h           # Base classes for spinless t-V models
│   ├── square.h                # Square lattice implementation
│   ├── honeycomb.h             # Honeycomb lattice implementation
│   ├── singleMajoranaHoneycomb.h  # Single Majorana honeycomb
│   ├── kitaevChain.h           # 1D Kitaev chain (p-wave SC chain)
│   ├── hubbardChain.h          # 1D Hubbard model (NEW)
│   └── pfapack/                # PFAPACK library for Pfaffian calculations
├── src/                        # Implementation files
│   ├── pfqmc.cpp               # PfQMC algorithm implementation
│   └── skewMatUtils.cpp        # Pfaffian utilities implementation
└── test/                       # Test files
```

## Key Classes and Architecture

### 1. Type Definitions (types.h:1)

```cpp
typedef std::complex<double> DataType;
typedef Eigen::MatrixXcd MatType;
typedef Eigen::VectorXcd cVecType;
typedef Eigen::VectorXi iVecType;
```

- Uses Eigen library with Intel MKL backend for efficient linear algebra
- Complex double precision for all quantum amplitudes
- `rdGenerator`: Random number generator for auxiliary fields and Metropolis updates

### 2. Operator Hierarchy (operator.h:1)

**Base Class: `Operator`** (operator.h:8)
- Abstract interface for operators in the path integral
- Virtual methods:
  - `left_multiply/right_multiply`: Matrix multiplication
  - `left_propagate/right_propagate`: Update Green's function
  - `update(MatType &g)`: Metropolis update of auxiliary fields
  - `stabilizedLeftMultiply(UDT &F)`: Stable multiplication with UDT

**Derived Classes:**

1. **DenseOperator** (operator.h:38)
   - Represents kinetic term: `exp(-dt * H_kinetic)`
   - Stores matrix, inverse, and associated Green's function
   - No auxiliary fields to update

2. **SpinlessVOperator** (spinless_tV.h:140)
   - Represents interaction term via HS transformation
   - Stores auxiliary Ising fields `s` (±1 on each bond)
   - Implements Metropolis single-flip updates via `singleFlip()`
   - Generates interaction matrices using `InteractionBGenerator()`

### 3. Model Configuration Classes

**Base Class: `SpinlessTvUtils`** (spinless_tV.h:7)
- Stores model parameters: `Lx, Ly, dt, V, l, nDim`
- HS parameters: `lambdaV, chlV, shlV, thlV, etaM`
- Two HS schemes (hsScheme):
  - **Scheme 0** (hopping channel): Decouples `n_i n_j` as `(γ^1_i γ^1_j + γ^2_i γ^2_j)`
  - **Scheme 1** (density channel): Decouples as `(n_i - n_j)^2`
- Virtual method `aux2MajoranaIdx()`: Maps auxiliary field index to Majorana indices
- Methods to generate interaction matrices: `InteractionTanhGenerator()`, `InteractionBGenerator()`

**Derived Classes:**

1. **SpinlessTvChainUtils** (kitaevChain.h:10)
   - 1D chain with hopping, p-wave pairing (`delta`), and chemical potential (`mu`)
   - `boundaryType`: 0 = PBC, 1 = OBC
   - Implements `KineticGenerator()`: Builds kinetic Hamiltonian matrix
   - Observable calculations: `EdgeCorrelator()`, `Z2FermionParity()`, `StructureFactorCDW()`

2. **SpinlessTvSquareUtils** (square.h)
   - 2D square lattice implementation

3. **SpinlessTvHoneycombUtils** (honeycomb.h)
   - 2D honeycomb lattice implementation

### 4. Walker Classes

**Base Class: `Spinless_tV`** (spinless_tV.h:406)
- Container for operator array `op_array`
- Destructor handles cleanup

**Derived Classes:**

1. **Chain_tV** (kitaevChain.h:404)
   - Constructor builds operator sequence:
     - `exp(-dt/2 * H_K)` (half kinetic step)
     - For each time slice: `exp(-dt * H_K)` → `V_bond0` → `V_bond1`
     - `exp(-dt/2 * H_K)` (half kinetic step)
   - Splits bonds into two types for checkerboard decomposition

2. **Square_tV**, **Honeycomb_tV**: Similar structure for 2D lattices

### 5. PfQMC Class (pfqmc.h:7)

**Main algorithm driver:**

```cpp
class PfQMC {
    int stb;                              // Stabilization interval
    int nDim;                             // System dimension (# Majorana modes)
    MatType g;                            // Equal-time Green's function
    std::vector<Operator*> op_array;      // Operator sequence
    std::vector<bool> need_stabilization; // Stabilization checkpoints
    std::vector<UDT> udtL, udtR;          // Left/right UDT decompositions
    DataType sign;                        // Current sign

    void rightInit();                     // Initialize right sweep
    void leftInit();                      // Initialize left sweep
    void rightSweep();                    // Sweep left-to-right
    void leftSweep();                     // Sweep right-to-left
    DataType getSignRaw();                // Compute sign from Pfaffian
};
```

**Workflow:**
1. **Initialization**: `rightInit()` or `leftInit()` computes initial Green's function
2. **Sweeps**: Alternates `rightSweep()` and `leftSweep()`
   - For each operator: attempt Metropolis updates of auxiliary fields
   - Propagate Green's function forward/backward in time
   - Apply UDT stabilization at checkpoints
3. **Measurements**: Compute observables from Green's function `g` after each sweep

### 6. UDT Decomposition (qr_udt.h:6)

**Numerical stabilization technique:**

```cpp
class UDT {
    MatType U;   // Unitary matrix
    dVecType D;  // Diagonal positive values
    MatType T;   // Upper triangular matrix

    void onePlusInv(MatType &g);  // Compute g = 2(1 + UDT)^{-1}
};
```

- Decomposes matrix as `A = U * D * T`
- Prevents overflow/underflow in product of many matrices
- Used at stabilization checkpoints to reconstruct Green's function

## Current Models Implementation

### 1. Kitaev Chain (kitaevChain.h:10)

**Hamiltonian:**
```
H = -t Σ_i (c†_i c_{i+1} + h.c.) + Δ Σ_i (c†_i c†_{i+1} + h.c.) - μ Σ_i n_i + V Σ_i n_i n_{i+1}
```

**Majorana representation:**
- Each site has 2 Majorana modes: `γ^1_i, γ^2_i`
- Fermion: `c_i = (γ^1_i + i γ^2_i)/2`
- Kinetic term: `KineticGenerator()` (kitaevChain.h:47)
- Interaction: HS transformation with auxiliary fields on bonds

**Observables:**
- Energy: `energyFromGreensFunc()` (kitaevChain.h:309)
- Edge correlator: `EdgeCorrelator()` (kitaevChain.h:129) - Majorana edge mode correlation
- Z2 fermion parity: `Z2FermionParity()` (kitaevChain.h:238)
- CDW structure factor: `StructureFactorCDW()` (kitaevChain.h:136)

### 2. Square Lattice (square.h)

Spinless t-V model on 2D square lattice with optional p+ip pairing.

### 3. Honeycomb Lattice (honeycomb.h, singleMajoranaHoneycomb.h)

Spinless t-V model on 2D honeycomb lattice, with variant using single Majorana per site.

### 4. 1D Hubbard Model (hubbardChain.h)

**Hamiltonian:**
```
H = -t Σ_{i,σ} (c†_{i,σ} c_{i+1,σ} + h.c.) + U Σ_i n_{i,↑} n_{i,↓} - μ Σ_{i,σ} n_{i,σ}
```

**Majorana representation:**
- Each site has 4 Majorana modes: `γ^1_{i,↑}, γ^2_{i,↑}, γ^1_{i,↓}, γ^2_{i,↓}`
- System dimension: `nDim = 4L` (vs 2L for spinless models)
- Fermion: `c_{i,σ} = (γ^1_{i,σ} + i γ^2_{i,σ})/2`
- Density: `n_{i,σ} = (1 - i γ^1_{i,σ} γ^2_{i,σ})/2`

**Kinetic term** (spin-conserving hopping):
```
-t (c†_{i,σ} c_{i+1,σ} + h.c.) = i·t/2 (γ^1_{i,σ} γ^2_{i+1,σ} - γ^2_{i,σ} γ^1_{i+1,σ})
```

**Interaction term** (on-site Hubbard U):
```
U n_{i,↑} n_{i,↓} = U/4 (1 - i γ^1_{i,↑} γ^2_{i,↑})(1 - i γ^1_{i,↓} γ^2_{i,↓})
```
- Decoupled via HS transformation with auxiliary Ising fields on each site
- Two HS schemes available (hsScheme = 0 or 1)

**Observables:**
- Energy: `energyFromGreensFunc()` (hubbardChain.h:116)
- Particle number: `particleNumber()` (hubbardChain.h:176)
- Spin structure factor: `spinStructureFactor(g, q)` (hubbardChain.h:189)
- Charge structure factor: `chargeStructureFactor(g, q)` (hubbardChain.h:236)

**Usage:** See `HUBBARD_GUIDE.md` for detailed usage instructions and examples.

## Running Simulations

### Command Line Interface (main.cpp:372)

```bash
# Square lattice
./build/main--square

# Honeycomb lattice
./build/main --honeycomb

# Single Majorana honeycomb (with MPI)
mpirun ./build/main --MRhoneycomb Lx Ly LTau dt V nthreads nseed evaluationLength filepath

# Kitaev chain
mpirun ./build/main --chain Lx LTau dt V delta nthreads nseed evaluationLength filepath mu hsScheme boundary

# 1D Hubbard model (NEW)
mpirun ./build/main --hubbard Lx LTau dt U nthreads nseed evaluationLength filepath mu hsScheme boundary
```

### Key Parameters

- `Lx, Ly`: System size (number of sites)
- `LTau`: Number of imaginary time slices
- `dt`: Time step (Trotter error ~ dt^2)
- `V`: Interaction strength
- `delta`: Pairing amplitude (for Kitaev chain)
- `U`: On-site interaction strength (for Hubbard model)
- `mu`: Chemical potential
- `hsScheme`: Hubbard-Stratonovich scheme (0 or 1)
- `boundary`: Boundary condition (0 = PBC, 1 = OBC)
- `stabilizationTime`: UDT stabilization interval (typically 10)
- `thermalLength`: Thermalization sweeps (typically 200-1000)
- `evaluationLength`: Measurement sweeps (typically 1000-10000)

### Typical Workflow (main.cpp:214)

```cpp
// 1. Create configuration
SpinlessTvChainUtils config(Lx, dt, V, LTau, boundary, delta, mu, hsScheme);

// 2. Create random generator and walker
rdGenerator rd(nseed);
Chain_tV walker(&config, &rd);

// 3. Create PfQMC solver
PfQMC pfqmc(&walker, stabilizationTime);

// 4. Thermalization
for (int i = 0; i < thermalLength; i++) {
    pfqmc.rightSweep();
    pfqmc.leftSweep();
}

// 5. Measurement
for (int i = 0; i < evaluationLength; i++) {
    pfqmc.rightSweep();
    pfqmc.leftSweep();

    // Check/fix sign
    if (i % 20 == 0) {
        signRaw = pfqmc.getSignRaw();
        if (abs(sign - signRaw) > threshold) {
            fixSign(sign, signRaw, threshold);
        }
    }

    // Measure observables
    energy = config.energyFromGreensFunc(pfqmc.g);
    // ... other measurements
}
```

## Adding New Models

To implement a new model (e.g., 1-D Hubbard model):

1. **Create Utils class** (inherits `SpinlessTvUtils`):
   - Define lattice geometry
   - Implement `aux2MajoranaIdx()` mapping
   - Implement `KineticGenerator()` for kinetic Hamiltonian
   - Implement observable calculations

2. **Create Walker class** (inherits `Spinless_tV`):
   - Constructor builds `op_array` with operator sequence
   - Alternate kinetic and interaction operators
   - Use checkerboard decomposition for interaction

3. **Add main function** in `main.cpp`:
   - Parse command line arguments
   - Create configuration and walker
   - Run PfQMC simulation
   - Output results

4. **Add header files** in `inc/` directory

## Implementation Notes for 1-D Hubbard Model in Majorana Representation

### Standard Hubbard Model
```
H = -t Σ_{i,σ} (c†_{i,σ} c_{i+1,σ} + h.c.) + U Σ_i n_{i,↑} n_{i,↓}
```

### Challenges for Majorana Representation

The Hubbard model has **two spin species** (↑, ↓), which requires:

1. **Four Majorana modes per site**: `γ^1_{i,↑}, γ^2_{i,↑}, γ^1_{i,↓}, γ^2_{i,↓}`
2. **Spin-dependent interaction**: `U n_{i,↑} n_{i,↓}` couples different spin sectors
3. **System dimension**: `nDim = 4L` (vs `2L` for spinless model)

### Possible Approaches

#### Option 1: Direct Majorana Hubbard Model
- Use 4 Majorana modes per site
- Implement HS transformation for on-site interaction
- Similar structure to existing `Chain_tV` but with doubled Majorana sector
- Sign problem may appear depending on parameters

#### Option 2: Attractive Hubbard Model (U < 0)
- Can be mapped to spinless model with pairing
- Better sign structure
- Already partially implemented via `kitaevChain.h` with `delta` parameter

#### Option 3: Use Existing Kitaev Chain as Template
- The current `Chain_tV` class (kitaevChain.h:404) is very close to what's needed
- Has hopping (t), pairing (delta), chemical potential (mu), and NN interaction (V)
- To get Hubbard model: need to modify interaction term and add spin degree of freedom

### Recommended Implementation Strategy

**Use the existing Kitaev chain as starting point**, extend to include:
1. Spin index in Majorana indexing
2. Modified interaction term for on-site U interaction
3. Two HS schemes for U term (similar to current V term)
4. Adjust observable calculations for spin-resolved quantities

The architecture is already set up to handle this - main changes needed:
- New `HubbardChainUtils` class extending `SpinlessTvUtils`
- Implement 4-Majorana indexing scheme
- New interaction operator for on-site U term
- New `Hubbard_tV` walker class

## Dependencies

- **Eigen3**: Matrix library
- **Intel MKL**: BLAS/LAPACK backend
- **PFAPACK**: Fortran library for Pfaffian calculations
- **MPI**: Parallel execution across multiple seeds
- **OpenMP**: Thread-level parallelism

## Build Instructions

```bash
# Setup Intel oneAPI environment
source /opt/intel/oneapi/setvars.sh

# Build PFAPACK
cd inc/pfapack/fortran && make mFC=ifx
cd ../c_interface && make

# Build main code
cd /home/hzxiaxz/PfQMC/build
cmake .. -DEIGEN3_INCLUDE_DIR=/home/hzxiaxz/eigen-5.0.0
make -j 8
```

## References

- Paper: [arXiv:2408.10311](https://arxiv.org/abs/2408.10311) - PfQMC algorithm description
- PFAPACK: [arXiv:1102.3440](https://arxiv.org/abs/1102.3440) - Pfaffian computation library
- dont use mpirun when testing

# Repository Guidelines

## Project Structure & Module Organization
Core simulation logic lives in `src/` (`pfqmc.cpp`, `skewMatUtils.cpp`), with shared headers and the PFAPACK Fortran/C bridge under `inc/` (build `inc/pfapack/{fortran,c_interface}` before linking). `main.cpp` drives CLI experiments, while executable and object artifacts land in `bin/` and `obj/`. GoogleTest fixtures are in `test/`, and reference derivations plus Hubbard study notes sit in `_research/` and the `HUBBARD_*.md` sheets—treat them as canonical physics context when adding operators or observables.

## Build, Test, and Development Commands
- `cd inc/pfapack/fortran && make mFC=ifx` then `cd ../c_interface && make` to refresh Pfaffian libraries whenever you change compilers.
- `./local-build.sh [--debug|--release|--clean-only]` sources Intel oneAPI, wipes `build/`, and runs CMake with the Eigen include path; defaults to Debug, so pass `--release` for benchmarks or `--` to forward extra CMake flags.
- `cmake -S . -B build -DEIGEN3_INCLUDE_DIR=/path/to/eigen3` is the manual alternative if you want a separate build tree; follow with `cmake --build build --target main` (or `main_test`).
- `ctest --test-dir build --output-on-failure` or `./build/main_test` executes the full GoogleTest suite after a configure/build.
- `make clean` only applies when using the legacy Makefile rather than CMake; the script already removes stale artifacts each run.

## Coding Style & Naming Conventions
Target C++17 with 4-space indentation, LLVM brace placement, and descriptive camel-case for types (`PfQMC`, `SpinlessTvUtils`) plus lower camelCase for methods (`leftSweep`). Keep Eigen aliases (`MatType`, `DataType`) in headers, prefer `const` references for matrices, and document stabilization logic with brief comments rather than block prose. Run `clang-format -style="{BasedOnStyle: llvm, IndentWidth: 4}" -i <files>` before committing and ensure headers under `inc/` stay self-contained. Clangd users should keep the repo-root `.clangd` to pick up Intel oneAPI include paths so diagnostics mirror the `mpiicpx`/MKL toolchain.

## Testing Guidelines
Unit tests rely on GoogleTest (see `test/main_test.cpp`, `test/squareLatticeTest.cpp`). Mirror production namespaces in test fixture names (`SquareLatticeTest.CheckPropagation`) and cover both numerical accuracy and sign-handling paths. Any change to propagator math must add or extend a deterministic test; stochastic pieces should expose a seeded path. Execute `ctest -R <pattern>` for focused runs and avoid skipping tests in CI—coverage is thin, so every new feature should introduce at least one assertion that protects it.

## Commit & Pull Request Guidelines
History shows concise, imperative commits ("fix sign in Number", "remove two point term when mu=0"); follow that style, referencing issues as `Fix: ... (#42)` when relevant. Each PR should summarize the physics change, list required environment variables (e.g., lattice size, `beta`), and include `ctest` output plus any generated `.out` files that prove convergence. When touching numerics, attach small tables or plots comparing prior runs to justify deltas before requesting review.

## Security & Configuration Tips
Keep MKL and PFAPACK builds in sync with the compiler used by CMake/Make to avoid ABI mismatches. Do not commit Intel oneAPI environment dumps or raw lattice datasets larger than ~5 MB; store them under a reproducible script in `_research/` instead. Validate new environment variables in `README.md` and scrub secrets from shell snippets.
- use ./local_build.sh to compile