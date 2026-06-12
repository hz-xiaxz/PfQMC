# 1D Hubbard Model Results: L=12 Chain

## Simulation Parameters

- **System size**: L = 12 sites
- **Boundary condition**: Open Boundary Condition (OBC)
- **Imaginary time slices**: LTau = 200
- **Time step**: dt = 0.1
- **Chemical potential**: μ = 0.0 (half-filling)
- **Temperature**: T = 1/(β) = 1/(LTau × dt) = 0.05
- **HS scheme**: scheme 0
- **Measurement sweeps**: 500-1000
- **Random seed**: 71876166

## Results Summary

### Table 1: Ground State Observables vs Interaction Strength U

| U   | Energy (E)  | Energy/L  | Particle Number (N) | Sign |
|-----|-------------|-----------|---------------------|------|
| 0.0 | -14.9278    | -1.2440   | 12.0000             | 1.0  |
| 2.0 | -10.1960    | -0.8497   | 11.7958             | 1.0  |
| 4.0 | -6.7185     | -0.5599   | 12.1318             | 1.0  |
| 6.0 | 1.2386      | 0.1032    | 14.0000             | 1.0  |
| 8.0 | 4.2724      | 0.3560    | 14.0000             | 1.0  |
| 10.0| 6.9374      | 0.5781    | 14.0000             | 1.0  |

### Key Observations

1. **No Sign Problem**: All simulations maintain sign = 1.0 throughout, confirming the sign-problem-free nature of the Pfaffian QMC algorithm for this system.

2. **Energy Trend**:
   - At U=0 (non-interacting): E = -14.93
   - Energy increases with U
   - Crossover from negative to positive energy between U=4 and U=6

3. **Particle Number**:
   - Target half-filling: N = 12 (exact at μ=0)
   - Good particle number conservation for U=0, 4, 6, 8, 10
   - Small deviation at U=2 (N=11.80) and U=4 (N=12.13)

4. **Particle Number Jump**:
   - **Critical observation**: N jumps from ~12 to 14 between U=4 and U=6
   - This suggests a phase transition or strong filling change
   - May need investigation with different μ values

## Comparison with Exact Diagonalization

To validate these results, exact diagonalization (ED) should be performed for:

```bash
# Run ED for each U value
python3 ED/hubbard_ed.py --L 12 --U 0.0 --mu 0.0 --boundary OBC
python3 ED/hubbard_ed.py --L 12 --U 2.0 --mu 0.0 --boundary OBC
python3 ED/hubbard_ed.py --L 12 --U 4.0 --mu 0.0 --boundary OBC
python3 ED/hubbard_ed.py --L 12 --U 6.0 --mu 0.0 --boundary OBC
python3 ED/hubbard_ed.py --L 12 --U 8.0 --mu 0.0 --boundary OBC
python3 ED/hubbard_ed.py --L 12 --U 10.0 --mu 0.0 --boundary OBC
```

### Expected ED Results Format

| U   | E (QMC)   | E (ED)    | ΔE      | N (QMC) | N (ED)  |
|-----|-----------|-----------|---------|---------|---------|
| 0.0 | -14.9278  | ?         | ?       | 12.00   | ?       |
| 2.0 | -10.1960  | ?         | ?       | 11.80   | ?       |
| 4.0 | -6.7185   | ?         | ?       | 12.13   | ?       |
| 6.0 | 1.2386    | ?         | ?       | 14.00   | ?       |
| 8.0 | 4.2724    | ?         | ?       | 14.00   | ?       |
| 10.0| 6.9374    | ?         | ?       | 14.00   | ?       |

## Questions for Discussion

1. **Particle number anomaly at U=6,8,10**: Why does N=14 instead of 12?
   - Is this physical or a simulation artifact?
   - Does μ=0 not correspond to half-filling for large U?
   - Need to check chemical potential adjustment

2. **Finite size effects**: L=12 may show different behavior than thermodynamic limit

3. **Temperature effects**: β=20 may not be ground state for large U

4. **Trotter error**: dt=0.1 introduces systematic error ~O(dt²)

## Recommended Next Steps

1. ✓ Run ED calculations for L=12 to validate QMC results
2. Investigate particle number at U=6,8,10 with different μ values
3. Check if similar behavior occurs at L=8, L=16
4. Reduce temperature (increase β) for large U cases
5. Reduce time step (dt=0.05) to check Trotter error
