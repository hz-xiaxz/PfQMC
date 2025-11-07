#ifndef HUBBARDCHAIN_H
#define HUBBARDCHAIN_H

#include "operator.h"
#include "skewMatUtils.h"
#include "spinless_tV.h"
#include "types.h"
#include <iostream>

/**
 * 1D Hubbard Model in Majorana Representation
 *
 * Hamiltonian: H = -t Σ_{i,σ} (c†_{i,σ} c_{i+1,σ} + h.c.) + U Σ_i n_{i,↑} n_{i,↓} - μ Σ_{i,σ} n_{i,σ}
 *
 * Majorana representation:
 * - Each site i has 4 Majorana modes: γ^1_{i,↑}, γ^2_{i,↑}, γ^1_{i,↓}, γ^2_{i,↓}
 * - Fermion operator: c_{i,σ} = (γ^1_{i,σ} + i γ^2_{i,σ}) / 2
 * - Density: n_{i,σ} = c†_{i,σ} c_{i,σ} = (1 - i γ^1_{i,σ} γ^2_{i,σ}) / 2
 *
 * System dimension: nDim = 4 * Lx (4 Majorana modes per site)
 */

class HubbardChainUtils : public SpinlessTvUtils {
  public:
    int nsites;
    double U;         // On-site interaction strength
    double mu;        // Chemical potential
    int boundaryType; // 0: PBC, 1: OBC

    HubbardChainUtils(int _L, double _dt, double _U, int _l, int _boundary, double _mu = 0.0, int _hsScheme = 0)
        : SpinlessTvUtils(_L, 1, _dt, _U, _l, _L * 4, false, _hsScheme) {
        // Note: nDim = _L * 4 (4 Majorana modes per site: 2 spins × 2 Majoranas)
        boundaryType = _boundary;
        nsites = _L;
        U = _U;
        mu = _mu;

        // Recompute HS parameters for on-site interaction
        // For on-site U interaction, we use same HS transformation structure
        // Testing factor = 0.25 based on analysis
        lambdaV = acosh(exp(0.25 * U * dt));
        chlV = cosh(lambdaV);
        shlV = sinh(lambdaV);
        thlV = tanh(lambdaV);
        etaM = chlV * chlV;
    }

    /**
     * Majorana indexing: map (site, spin, majorana_species) to linear index
     * @param ix: site index (0 to nsites-1)
     * @param ispin: spin index (0 = up, 1 = down)
     * @param imaj: Majorana species (0 or 1)
     * @return linear index in [0, 4*nsites)
     */
    inline int majoranaCoord2Idx(int ix, int ispin, int imaj) const {
        return ix * 4 + ispin * 2 + imaj;
    }

    /**
     * Map auxiliary field index to Majorana pair for on-site interaction
     * For on-site Hubbard interaction U n_↑ n_↓
     *
     * hsScheme = 0: Decouple as exp(λ s [i γ^1_↑ γ^2_↑ + i γ^1_↓ γ^2_↓])
     * hsScheme = 1: Decouple as exp(λ s [i γ^1_↑ γ^2_↑ - i γ^1_↓ γ^2_↓])
     *
     * @param idAux: site index (on-site interaction)
     * @param imaj: majorana species index (not used for on-site, kept for interface)
     * @param bType: bond type (not used for on-site)
     * @param idx1, idx2: output Majorana indices
     */
    inline void aux2MajoranaIdx(int idAux, int imaj, int bType, int &idx1, int &idx2) const override {
        // For on-site interaction, idAux is just the site index
        // We need to couple the two spins on the same site
        if (hsScheme == 0) {
            // Scheme 0: couple same-spin Majoranas
            // This gives: i γ^1_{i,σ} γ^2_{i,σ} terms
            idx1 = majoranaCoord2Idx(idAux, imaj, 0); // γ^1_{i,σ}
            idx2 = majoranaCoord2Idx(idAux, imaj, 1); // γ^2_{i,σ}
        } else {
            // Scheme 1: alternative coupling
            idx1 = majoranaCoord2Idx(idAux, imaj, 0);
            idx2 = majoranaCoord2Idx(idAux, imaj, 1);
        }
    }

    /**
     * Generate kinetic Hamiltonian matrix H_kinetic in Majorana basis
     * H_kinetic = -t Σ_{i,σ} (c†_{i,σ} c_{i+1,σ} + h.c.) - μ Σ_{i,σ} n_{i,σ}
     *
     * In Majorana basis:
     * Hopping: -t (c†_{i,σ} c_{i+1,σ} + h.c.) = i·t/2 (γ^1_{i,σ} γ^2_{i+1,σ} - γ^2_{i,σ} γ^1_{i+1,σ})
     * Chemical potential: -μ n_{i,σ} = -μ/2 (1 - i γ^1_{i,σ} γ^2_{i,σ})
     */
    inline void KineticGenerator(MatType &H) const {
        H.setZero();
        int idx1, idx2;
        DataType tmp = (0.5i);    // Factor for hopping
        DataType tmpMu = (0.5i) * mu; // Factor for chemical potential

        // Hopping term: -t Σ_{i,σ} (c†_{i,σ} c_{i+1,σ} + h.c.)
        // In Majorana: i·t/2 (γ^1_{i,σ} γ^2_{i+1,σ} - γ^2_{i,σ} γ^1_{i+1,σ})
        int L_bonds = (boundaryType == 0) ? Lx : (Lx - 1); // PBC: Lx bonds, OBC: Lx-1 bonds

        for (int i = 0; i < L_bonds; i++) {
            for (int ispin = 0; ispin < 2; ispin++) {
                // γ^1_{i,σ} γ^2_{i+1,σ} term
                idx1 = majoranaCoord2Idx(i, ispin, 0);
                idx2 = majoranaCoord2Idx((i + 1) % Lx, ispin, 1);
                H(idx1, idx2) += tmp;
                H(idx2, idx1) += -tmp;

                // -γ^2_{i,σ} γ^1_{i+1,σ} term
                idx1 = majoranaCoord2Idx(i, ispin, 1);
                idx2 = majoranaCoord2Idx((i + 1) % Lx, ispin, 0);
                H(idx1, idx2) += -tmp;
                H(idx2, idx1) += +tmp;
            }
        }

        // Chemical potential term: -μ Σ_{i,σ} n_{i,σ} = -μ/2 Σ_{i,σ} (1 - i γ^1_{i,σ} γ^2_{i,σ})
        // Only the i γ^1_{i,σ} γ^2_{i,σ} part contributes to H (constant part drops out)
        for (int i = 0; i < Lx; i++) {
            for (int ispin = 0; ispin < 2; ispin++) {
                idx1 = majoranaCoord2Idx(i, ispin, 0);
                idx2 = majoranaCoord2Idx(i, ispin, 1);
                H(idx1, idx2) += -tmpMu;
                H(idx2, idx1) += +tmpMu;
            }
        }
    }

    /**
     * Calculate energy from Green's function
     * E = <H> = -t Σ_{i,σ} <c†_{i,σ} c_{i+1,σ} + h.c.> + U Σ_i <n_{i,↑} n_{i,↓}> - μ Σ_{i,σ} <n_{i,σ}>
     */
    inline DataType energyFromGreensFunc(const MatType &g) {
        DataType r = 0.0;
        DataType tmp = (0.5i);
        DataType tmpMu = (0.5i) * mu;
        int idx1, idx2;

        // Kinetic energy: hopping term
        int L_bonds = (boundaryType == 0) ? Lx : (Lx - 1);
        for (int i = 0; i < L_bonds; i++) {
            for (int ispin = 0; ispin < 2; ispin++) {
                // <γ^1_{i,σ} γ^2_{i+1,σ}>
                idx1 = majoranaCoord2Idx(i, ispin, 0);
                idx2 = majoranaCoord2Idx((i + 1) % Lx, ispin, 1);
                r += tmp * g(idx1, idx2);

                // <-γ^2_{i,σ} γ^1_{i+1,σ}>
                idx1 = majoranaCoord2Idx(i, ispin, 1);
                idx2 = majoranaCoord2Idx((i + 1) % Lx, ispin, 0);
                r -= tmp * g(idx1, idx2);
            }
        }

        // Chemical potential: -μ Σ_{i,σ} <n_{i,σ}>
        for (int i = 0; i < Lx; i++) {
            for (int ispin = 0; ispin < 2; ispin++) {
                idx1 = majoranaCoord2Idx(i, ispin, 0);
                idx2 = majoranaCoord2Idx(i, ispin, 1);
                r -= tmpMu * g(idx1, idx2);
            }
        }

        // Interaction energy: U Σ_i <n_{i,↑} n_{i,↓}>
        // n_{i,σ} = (1 - i γ^1_{i,σ} γ^2_{i,σ}) / 2
        // n_{i,↑} n_{i,↓} = [1 - i(γ^1_{i,↑} γ^2_{i,↑} + γ^1_{i,↓} γ^2_{i,↓})
        //                     - γ^1_{i,↑} γ^2_{i,↑} γ^1_{i,↓} γ^2_{i,↓}] / 4
        int idxu1, idxu2, idxd1, idxd2;
        DataType tmpU = 0.25 * U;

        // Constant term: U/4 per site
        r += tmpU * double(Lx);

        for (int i = 0; i < Lx; i++) {
            idxu1 = majoranaCoord2Idx(i, 0, 0); // γ^1_{i,↑}
            idxu2 = majoranaCoord2Idx(i, 0, 1); // γ^2_{i,↑}
            idxd1 = majoranaCoord2Idx(i, 1, 0); // γ^1_{i,↓}
            idxd2 = majoranaCoord2Idx(i, 1, 1); // γ^2_{i,↓}

            // Two-point terms: -U/4 (<i γ^1_↑ γ^2_↑> + <i γ^1_↓ γ^2_↓>)
            r -= tmpU * (1.0i) * (g(idxu1, idxu2) + g(idxd1, idxd2));

            // <n_{i,↑} n_{i,↓}> four-point term with ALL Wick contractions
            // Following Kitaev chain pattern (kitaevChain.h:386-388): signs are +, +, -

            DataType nn_interaction = tmpU * (
                + g(idxu1, idxd1) * g(idxu2, idxd2)  // +<γ^1_↑ γ^1_↓><γ^2_↑ γ^2_↓>
                + g(idxu1, idxd2) * g(idxu2, idxd1)  // +<γ^1_↑ γ^2_↓><γ^2_↑ γ^1_↓>
                - g(idxu1, idxu2) * g(idxd1, idxd2)  // -<γ^1_↑ γ^2_↑><γ^1_↓ γ^2_↓> (disconnected)
            );
            r -= nn_interaction;  // Note: minus sign for four-point term
        }

        return r;
    }

    /**
     * Calculate total particle number <N> = Σ_{i,σ} <n_{i,σ}>
     */
    inline DataType particleNumber(const MatType &g) const {
        DataType n_total = 0.0;
        int idx1, idx2;

        for (int i = 0; i < Lx; i++) {
            for (int ispin = 0; ispin < 2; ispin++) {
                idx1 = majoranaCoord2Idx(i, ispin, 0);
                idx2 = majoranaCoord2Idx(i, ispin, 1);
                // <n_{i,σ}> = 1/2 - i/2 <γ^1_{i,σ} γ^2_{i,σ}>
                n_total += 0.5 * (1.0 - (1.0i) * g(idx1, idx2));
            }
        }

        return n_total;
    }
};

    /**
     * Spin structure factor S(q) = (1/L) Σ_{i,j} e^{iq(i-j)} <S^z_i S^z_j>
     * where S^z_i = (n_{i,↑} - n_{i,↓}) / 2
     */
    // inline DataType spinStructureFactor(const MatType &g, double q) const {
    //     DataType Sq = 0.0;

    //     for (int i = 0; i < Lx; i++) {
    //         for (int j = 0; j < Lx; j++) {
    //             int idxu1_i = majoranaCoord2Idx(i, 0, 0);
    //             int idxu2_i = majoranaCoord2Idx(i, 0, 1);
    //             int idxd1_i = majoranaCoord2Idx(i, 1, 0);
    //             int idxd2_i = majoranaCoord2Idx(i, 1, 1);

    //             int idxu1_j = majoranaCoord2Idx(j, 0, 0);
    //             int idxu2_j = majoranaCoord2Idx(j, 0, 1);
    //             int idxd1_j = majoranaCoord2Idx(j, 1, 0);
    //             int idxd2_j = majoranaCoord2Idx(j, 1, 1);

    //             // S^z_i = (n_{i,↑} - n_{i,↓}) / 2
    //             // <S^z_i S^z_j> = (1/4)[<n↑n↑> - <n↑n↓> - <n↓n↑> + <n↓n↓>]
    //             // Only connected Wick contractions (exclude disconnected diagrams)

    //             DataType block_upup = -g(idxu1_i, idxu1_j) * g(idxu2_i, idxu2_j)
    //                                  + g(idxu1_i, idxu2_j) * g(idxu2_i, idxu1_j);
    //             DataType block_updn = -g(idxu1_i, idxd1_j) * g(idxu2_i, idxd2_j)
    //                                  + g(idxu1_i, idxd2_j) * g(idxu2_i, idxd1_j);
    //             DataType block_dnup = -g(idxd1_i, idxu1_j) * g(idxd2_i, idxu2_j)
    //                                  + g(idxd1_i, idxu2_j) * g(idxd2_i, idxu1_j);
    //             DataType block_dndn = -g(idxd1_i, idxd1_j) * g(idxd2_i, idxd2_j)
    //                                  + g(idxd1_i, idxd2_j) * g(idxd2_i, idxd1_j);

    //             // Spin: +upup -updn -dnup +dndn
    //             DataType corr = 0.25 * (block_upup - block_updn - block_dnup + block_dndn);

    //             Sq += corr * exp(1.0i * q * double(i - j));
    //         }
    //     }

    //     return Sq / double(Lx);
    // }

    /**
     * Charge structure factor N(q) = (1/L) Σ_{i,j} e^{iq(i-j)} <n_i n_j>
     * where n_i = n_{i,↑} + n_{i,↓}
     */
//     inline DataType chargeStructureFactor(const MatType &g, double q) const {
//         DataType Nq = 0.0;

//         for (int i = 0; i < Lx; i++) {
//             for (int j = 0; j < Lx; j++) {
//                 // Similar calculation as spin structure factor but with n_i = n_{i,↑} + n_{i,↓}
//                 int idxu1_i = majoranaCoord2Idx(i, 0, 0);
//                 int idxu2_i = majoranaCoord2Idx(i, 0, 1);
//                 int idxd1_i = majoranaCoord2Idx(i, 1, 0);
//                 int idxd2_i = majoranaCoord2Idx(i, 1, 1);

//                 int idxu1_j = majoranaCoord2Idx(j, 0, 0);
//                 int idxu2_j = majoranaCoord2Idx(j, 0, 1);
//                 int idxd1_j = majoranaCoord2Idx(j, 1, 0);
//                 int idxd2_j = majoranaCoord2Idx(j, 1, 1);

//                 // <n_i n_j> = <n↑n↑> + <n↑n↓> + <n↓n↑> + <n↓n↓>
//                 // Only connected Wick contractions (exclude disconnected diagrams)

//                 DataType block_upup = -g(idxu1_i, idxu1_j) * g(idxu2_i, idxu2_j)
//                                      + g(idxu1_i, idxu2_j) * g(idxu2_i, idxu1_j);
//                 DataType block_updn = -g(idxu1_i, idxd1_j) * g(idxu2_i, idxd2_j)
//                                      + g(idxu1_i, idxd2_j) * g(idxu2_i, idxd1_j);
//                 DataType block_dnup = -g(idxd1_i, idxu1_j) * g(idxd2_i, idxu2_j)
//                                      + g(idxd1_i, idxu2_j) * g(idxd2_i, idxu1_j);
//                 DataType block_dndn = -g(idxd1_i, idxd1_j) * g(idxd2_i, idxd2_j)
//                                      + g(idxd1_i, idxd2_j) * g(idxd2_i, idxd1_j);

//                 // Charge: +upup +updn +dnup +dndn (ALL added)
//                 DataType corr = 0.25 * (block_upup + block_updn + block_dnup + block_dndn);

//                 Nq += corr * exp(1.0i * q * double(i - j));
//             }
//         }

//         return Nq / double(Lx);
//     }
// };

/**
 * HubbardChain_tU: Walker class for 1D Hubbard model
 *
 * Constructs the operator array for QMC sampling:
 * op_array = [exp(-dt/2 * H_K), V_0, V_1, ..., V_{L-1}, exp(-dt * H_K), ..., exp(-dt/2 * H_K)]
 *
 * where:
 * - H_K is the kinetic Hamiltonian (hopping + chemical potential)
 * - V_i is the on-site interaction operator at site i
 */
class HubbardChain_tU : public Spinless_tV {
  public:
    const HubbardChainUtils *modelConfig;
    double dt;
    int l, nSites;
    rdGenerator *rd;

    HubbardChain_tU(HubbardChainUtils *_config, rdGenerator *_rd) {
        modelConfig = _config;
        dt = _config->dt;
        l = _config->l;
        nSites = _config->nsites;
        rd = _rd;

        nDim = 4 * nSites; // 4 Majorana modes per site

        // Generate kinetic Hamiltonian and its exponentials
        MatType Ht(nDim, nDim);
        const MatType identity = MatType::Identity(nDim, nDim);
        modelConfig->KineticGenerator(Ht);

        MatType Hcopy(Ht);
        MatType expK = expm(Hcopy, -dt);
        Hcopy = Ht;
        MatType expKhalf = expm(Hcopy, -dt / 2.0);

        // Compute signs associated with K and Khalf
        Hcopy = dt * Ht;
        DataType signK = signOfHamiltonian(Hcopy);
        Hcopy = (0.5 * dt) * Ht;
        DataType signKHalf = signOfHamiltonian(Hcopy);

        // Build operator array
        // Structure: [K/2, U_0, U_1, ..., U_{L-1}, K, U_0, ..., K/2]
        // For on-site interaction, all sites can be updated in one operator
        op_array = std::vector<Operator *>(2 * l + 1);
        iVecType *s;

        for (int i = 0; i < l; i++) {
            if (i == 0) {
                op_array[0] = new DenseOperator(expKhalf, signKHalf);
            } else {
                op_array[2 * i] = new DenseOperator(expK, signK);
            }

            // On-site interaction: one operator for all sites
            s = new iVecType(nSites);
            for (int k = 0; k < nSites; k++) {
                (*s)(k) = rd->rdZ2();
            }
            op_array[2 * i + 1] = new SpinlessVOperator(modelConfig, s, 0, rd);
        }
        op_array[2 * l] = new DenseOperator(expKhalf, signKHalf);
    }
};

#endif
