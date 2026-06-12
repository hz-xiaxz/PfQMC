# Honeycomb Replica Verification Results

**Date:** 2026-01-14  
**Test:** 4-Replica Honeycomb with TrivialSwapOperator vs Standard Single Replica

## Test Configuration

| Parameter | Value |
|-----------|-------|
| Lattice   | 2×2 Honeycomb |
| LTau      | 100 |
| dt        | 0.1 |
| Threads   | 7 |
| Seed      | 42 |
| evalLen   | 10 |
| Tolerance | 0.01 |

## Results Summary

| V   | Single Replica Energy | 4-Replica Energy | Diff      | Status  |
|-----|----------------------|------------------|-----------|---------|
| 0.1 | -6.047620            | -24.190500       | 2.0e-05   | ✅ PASS |
| 0.2 | -6.084440            | -24.337800       | 4.0e-05   | ✅ PASS |
| 0.3 | -6.123570            | -24.494300       | 2.0e-05   | ✅ PASS |
| 0.4 | -6.172350            | -24.689400       | 0.0       | ✅ PASS |
| 0.5 | -6.236030            | -24.944100       | 2.0e-05   | ✅ PASS |

## Analysis

The 4-replica energy is exactly 4× the single replica energy (within floating-point precision), confirming that the `TrivialSwapOperator` correctly implements an identity mixing operator.

**Observed discrepancies** (order of 10⁻⁵) are due to:
- Floating-point truncation in output formatting
- Machine precision differences between single and multi-replica code paths

All differences are well below the 0.01 tolerance threshold.

## Conclusion

✅ **All V values passed verification.** The multi-replica implementation with `TrivialSwapOperator` produces energies consistent with 4 independent single-replica simulations.
