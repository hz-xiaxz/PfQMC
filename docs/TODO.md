# Spinful Hubbard TODO

The Monte Carlo flow implemented in `inc/spinless_tV.h` (Majorana ladder, spinless) cannot be reused verbatim by `inc/hubbardChain.h`. We need a dedicated spinful implementation that mirrors the PfQMC walker semantics but honors the doubled site basis. Track the work here.

1. **Branch for investigation** – `git checkout -b feature/spinful-hubbard` *(DONE)* so exploratory work stays isolated until validated.
2. **Port core logic** – Clone the stabilizer/update code in `inc/spinless_tV.h` into a new spinful header/source pair (e.g., `inc/hubbardChainSpinful.h`), adjusting constructor signatures to accept spin-resolved operators. *(DONE via `inc/spinful_tV.h` and `HubbardChainUtils` now inheriting from it).* 
3. **Implement `majoranaCoord2Idx` (spinful)** – Encode `(site, spin, Majorana)` into contiguous indices. Proposed mapping: `idx = 4 * site + 2 * spin + majoranaComponent`, guaranteeing the even/odd structure required by PFAPACK.
4. **Review HS factors** – Confirm whether the spinful interaction should still use `lambdaV = acosh(exp(0.25 * U * dt))` or a modified prefactor (likely `0.5 * U * dt`). Derive the correct expression for up/down sectors before coding.
5. **Refactor HS schemes** – Reduce to a single transformation of the form `gamma_{u1} gamma_{u2} + gamma_{d1} gamma_{d2} - 1`, removing unused auxiliary-field mixes to keep stabilization predictable.
6. **Update propagation scheme** – Re-derive the Green’s function update with spin blocks; expect extra sign bookkeeping. When the formal calculation is ready, integrate it here along with new regression tests.

Add references to `HUBBARD_HS_CORRECT_DERIVATION.md` as steps land so future contributors have context.
