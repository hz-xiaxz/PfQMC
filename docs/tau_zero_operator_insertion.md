# Operator Insertion at $\tau=0$ in PfQMC

This document outlines the considerations and necessary changes for inserting a general operator at $\tau=0$ in the Auxiliary Field Quantum Monte Carlo (AFQMC) simulation.

## Current Architecture

In the current `PfQMC` implementation, $\tau=0$ acts as the boundary between the backward and forward propagation sweeps.

*   **Right Sweep ($0 \to \beta$):** The `mixingOp` (if present) is applied **first**, before any time slices.
*   **Left Sweep ($\beta \to 0$):** The `mixingOp` (if present) is applied **last**, after all time slices.

## Implementation Considerations

### 1. Propagation Order
When inserting a new operator $O_{new}$ at $\tau=0$, the order of application relative to the existing `mixingOp` is critical if they do not commute.

*   **Right Sweep:** If $O_{new}$ acts on the state before `mixingOp`, it must be applied before `mixingOp` in the code.
    *   *Code Location:* `PfQMC::rightSweep` (start of function).
*   **Left Sweep:** The application order must be the exact reverse of the Right Sweep.
    *   *Code Location:* `PfQMC::leftSweep` (end of function).

### 2. Stabilization (Numerical Stability)
Operators at the boundary can affect the condition number of the propagators.

*   **Immediate Stabilization:** If $O_{new}$ is singular (e.g., a projector) or significantly alters scales, a UDT decomposition (stabilization) should be performed immediately after its application.
*   **Current Behavior:** The code currently forces a stabilization (`mixUDT`) immediately after `mixingOp` in `rightSweep`.
*   **Recommendation:** Include $O_{new}$ in the `Aseg` accumulation and ensure a UDT decomposition happens right after.

### 3. Sampling and Updates
If $O_{new}$ contains auxiliary fields that require sampling (dynamic operator):

*   **Update Logic:** You must call `O_{new}->update(g)` during the **Left Sweep**.
*   **Location:** In `PfQMC::leftSweep`, at the point corresponding to its position in the operator chain (likely at $l=0$ or adjacent to `mixingOp`).
*   **Consequence:** Failing to update will leave the operator frozen at its initial configuration.

### 4. Initialization (`getSignRaw`)
The `getSignRaw` function constructs the initial Green's function and weight from scratch.

*   **Requirement:** Any operator added to the sweeps must also be added to `getSignRaw`.
*   **Code Location:** `PfQMC::getSignRaw`, specifically around where `mixingOp->stabilizedLeftMultiply(A)` is called.
*   **Risk:** Omitting this will result in an incorrect initial sign/phase and weight, invalidating the simulation start.

## Summary Checklist

To safely insert an operator at $\tau=0$:

1.  **Class Member:** Add the operator pointer to `PfQMC`.
2.  **Right Sweep:** Apply `left_propagate` at the start (respecting order).
3.  **Left Sweep:** Apply `right_propagate` at the end (respecting reverse order) and `update`.
4.  **Initialization:** Apply `stabilizedLeftMultiply` in `getSignRaw`.
5.  **Stabilization:** Ensure `UDT` is recomputed after these operations.
